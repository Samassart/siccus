import xarray as xr
from shapely import geometry
from equi7grid.equi7grid import Equi7Grid
from tqdm import tqdm

import dask
from dask.diagnostics import ProgressBar

import numpy as np
from datetime import datetime
import os

# Set the environment variables to limit resources


wgs_84_attrs = 'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AXIS["Latitude",NORTH],AXIS["Longitude",EAST],AUTHORITY["EPSG","4326"]]'



def read_bbox(
        scene:xr.DataArray,
        geom:geometry.multipolygon.MultiPolygon | geometry.polygon.Polygon,
    ):
    # in this case resample if needed
    try:
        crs_of_file = scene.crs.attrs["spatial_ref"]
    except:
        crs_of_file = scene.spatial_ref.attrs["spatial_ref"]
        
    if not crs_of_file == wgs_84_attrs:
        # make sure the data is reprojected to EPSG:4326
        # Reproject once and keep the information somewhere for quicker processing
        scene = scene.rio.reproject("EPSG:4326")

    scene.rio.set_spatial_dims(x_dim="x", y_dim="y", inplace=True)
    scene.rio.write_crs("epsg:4326", inplace=True)
    # crop with the (multi)polygon
    clipped_scene = scene.rio.clip(geom.values)
    return clipped_scene


def from_e7key_to_lonlat(tile, sampling):
    sampling = 500
    e7 = Equi7Grid(sampling)
    e7tile = e7.create_tile("EU500M_"+tile)
    txmin, tymin, txmax, tymax = e7tile._limits_m()
    xx = np.arange(txmin+sampling/2., txmax+sampling/2., sampling)
    yy = np.arange(tymax-sampling/2., tymin-sampling/2., -sampling)
    x, y = np.meshgrid(xx, yy)
    lons, lats = e7.EU.xy2lonlat(x,y)
    return lons, lats


# Define the function to compute percentile using np.percentile
def compute_percentile(arr, percentile):
    return np.percentile(arr, percentile, axis=0)

def process_tiff_files(
        tiff_files:list[str], 
        chunk_size:int, 
        datetime_list:list, 
        percentiles:list[int]
        ) -> np.array:
    # Create an empty list to store the DataArrays for each file
    da = None
    result_dict =  {}
    date_format = "%Y%m%d%H%M" 

    # Open all TIFF files and concatenate into a single DataArray along the time dimension
    print("Load the information for the quantile processing.")
    for i, tiff_file in enumerate(tiff_files):
        # get datetime information from tif file
        date_to_export = datetime_list[i]
        encoded_time = np.array([date_to_export], dtype='int')

        if da is None:
            da = xr.open_rasterio(tiff_file, chunks={'band': -1, 'x': chunk_size, 'y': chunk_size})
            da = xr.DataArray(da, coords={'time': encoded_time}, dims=('time', 'y', 'x'))

        # Open the TIFF file as a DataArray
        else:
            da_array = xr.open_rasterio(tiff_file, chunks={'band': -1, 'x': chunk_size, 'y': chunk_size})
            da_array = xr.DataArray(da_array, coords={'time': encoded_time}, dims=('time', 'y', 'x'))

            # use numpy 
            da = xr.concat([da, da_array], dim="time")
        
    # Try and rechunk the dask xarray dataarray
    rechunked_da = da.chunk(chunks=(-1, chunk_size, chunk_size))

    # Initialize array 
    for perc in percentiles:
        percentile_array_full = np.zeros((rechunked_da.x.sizes["x"], rechunked_da.y.sizes["y"]))
        result_dict[perc] = percentile_array_full
    
    # prepare shape from chunk information
    X, Y = percentile_array_full.shape

    # First load chunkwise
    for i in range(0, X, chunk_size):
        print(f"Processing for chunk {i}/{int(np.floor(X))}")
        for j in tqdm(range(0, Y, chunk_size)):
            rechunked_array = rechunked_da[:, i:i+chunk_size, j:j+chunk_size].compute()
            # apply on each element of our 3D chunked array
            for row_idx, row in enumerate(range(rechunked_array.values.shape[1])):
                for col_idx, element in enumerate(range(rechunked_array.values.shape[2])):
                    # solve index for global map
                    global_row_idx = i + row_idx
                    global_col_idx = j + col_idx

                    # extract individual timeseries
                    individual_timeseries = rechunked_array.values[:, row_idx, col_idx]
                    for perc in percentiles:
                        percentile_val = np.percentile(individual_timeseries, q=perc) # get an int
                        
                        # store your value in the right percentile index, at the right x,y value
                        result_dict[perc][global_row_idx,global_col_idx] = percentile_val


    return result_dict



import os
from datetime import datetime
import pandas as pd
from .temporal_function import (
    from_day_to_dekad,
    from_dekad_to_day,
    filter_dataframe,
    datetime_to_year_month_dekad
    )
from .spatial_function import read_bbox
from tqdm import tqdm
from typing import List
import rioxarray as rio
import xarray as xr
from shapely import geometry
import copy
import numpy as np
import matplotlib.pyplot as plt
wgs_84_attrs = 'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AXIS["Latitude",NORTH],AXIS["Longitude",EAST],AUTHORITY["EPSG","4326"]]'


class environmental_dataset:
    """
    Parent class for common reading functions for plotting
    mapping or analyzing datasets
    """
    def __init__(
        self,
        directory:os.PathLike
    ):
        self.directory = directory

    def read_bbox_at_datetime(
            self,
            datetime_to_match:datetime,
            geom:geometry.multipolygon.MultiPolygon | geometry.polygon.Polygon,
    ) -> xr.DataArray:
        """
        Take a date and bounding box as input
        Return a DataArray of the dekad for the 
        bounding box
        """
        # find closest time
        datetime_index_str = np.array(self.dataframe_values.index)
        datetime_index_dt = []

        for values in datetime_index_str:
            day_values = from_dekad_to_day(int(str(values)[6]))
            complete_str_value = f"{int(str(values)[0:6])}{day_values:02d}"
            complete_dt_value = datetime.strptime(
                complete_str_value,
                "%Y%m%d"
            )
            datetime_index_dt.append(complete_dt_value)
        dataframe_values_dt = copy.copy(self.dataframe_values) # copy to not override :)
        dataframe_values_dt.index = datetime_index_dt

        #
        dataframe_values_dt = dataframe_values_dt.sort_index()
        idx = dataframe_values_dt.index[
            dataframe_values_dt.index.get_indexer(
                [datetime_to_match], method='nearest'
                    )
                ]

        # read the .tif/.nc file
        file_to_read = dataframe_values_dt.loc[idx]
        scene = rio.open_rasterio(file_to_read.file_name.values[0])

        # apply the proper clip function
        clipped_scene = read_bbox(scene=scene, geom=geom)
        return clipped_scene

    def read_ts(
            self,
            lat:int,
            lon:int,
            start_time:datetime = None,
            end_time:datetime = None,
            nodata:int | float = None,
    ) -> pd.DataFrame:
        """
        Read and store a timeseries of the chirps dataset at a given lat/lon
        for given years/months/dekads period
        """
        # first convert the day into dekad*
        try:
            thresh_start = datetime_to_year_month_dekad(start_time)
            thresh_end = datetime_to_year_month_dekad(end_time)

            # filter the dataset 
            filtered_dataset = filter_dataframe(
                self.dataframe_values,
                thresh_start,
                thresh_end
                    )

        except:
            # no dates given just take everything
            filtered_dataset = self.dataframe_values
        timeseries_to_extract = {}
        preprocessed_data = None

        for file in filtered_dataset.iterrows():
            # get time values
            file_to_extract_date = str(file[0])
            day_of_file = from_dekad_to_day(int(file_to_extract_date[6]))
            year_and_month_file = file_to_extract_date[0:6]
            file_datetime_date = datetime.strptime(
                f"{year_and_month_file}{day_of_file:02d}",
                "%Y%m%d"
                )
            # get file data
            file_to_extract = rio.open_rasterio(file[1].values[0])
            try:
                crs_of_file = file_to_extract.crs.attrs["spatial_ref"]
            except:
                crs_of_file = file_to_extract.spatial_ref.attrs["spatial_ref"]
                
            if not crs_of_file == wgs_84_attrs:
                # make sure the data is reprojected to EPSG:4326
                # Reproject once and keep the information somewhere for quicker processing
                if preprocessed_data is None:
                    file_to_extract_reproc = file_to_extract.rio.reproject("EPSG:4326")

                    # for this case extract from the reprojected
                    file_to_extract_value = file_to_extract_reproc.sel(
                            x=lon,
                            y=lat,
                            method="nearest"
                            )
                    lon_for_index= file_to_extract_value.x.values
                    lat_for_index= file_to_extract_value.y.values
                    # get x and y values
                    idx_x = int(np.where(file_to_extract_reproc.x.values == lon_for_index)[0])
                    idx_y = int(np.where(file_to_extract_reproc.y.values == lat_for_index)[0])

                    # change state of reprocessed date for next loop
                    preprocessed_data = True
                    specific_pixel_xr = file_to_extract_value

                else:
                    # extract back x,y values from positional index at first loop
                    lon_in_x_term = file_to_extract.x.values[idx_x]
                    lat_in_y_term = file_to_extract.y.values[idx_y]

                    # extract with the x and y values
                    specific_pixel_xr = file_to_extract.sel(
                            x=lon_in_x_term,
                            y=lat_in_y_term,
                            method="nearest"
                            )
            else:
                # simply load via lat/lon
                specific_pixel_xr = file_to_extract.sel(
                        x=lon,
                        y=lat,
                        method="nearest"
                        )

            # nodata management
            if not nodata is None:
                if specific_pixel_xr.values[0] == nodata:
                    timeseries_to_extract[file_datetime_date] = np.nan
                else:
                    timeseries_to_extract[file_datetime_date] = specific_pixel_xr.values[0]
            else:
                # simply read
                timeseries_to_extract[file_datetime_date] = specific_pixel_xr.values[0]
        # from dict to pandas        
        timeseries_to_extract_pd = pd.DataFrame(
            {f"{self.dataset}":timeseries_to_extract.values()},
            index=timeseries_to_extract.keys()
            )


        return timeseries_to_extract_pd
    
    def process_percentiles(
        self,
        output_directory:os.PathLike,
        geom:geometry.multipolygon.MultiPolygon | geometry.polygon.Polygon = None,
        chunk_size:int = None,
        percentiles:list[float] = [], 
        no_data:list = None,
        region_id:str = "full_data",
        large_process:bool = True,
    ):
        """
        Process the high and low percentiles for the defined area
        """
        if large_process:
            raise ValueError("Use the process percentile chunk instead")
        
        self.percentiles_layer = {}
        concatenated_array = None

        # first extract the files data
        print(f"load the array for {self.dataset} over the study area: {region_id}...")
        t = 0
        for files in tqdm(self.dataframe_values["file_name"]):
            date = self.dataframe_values.loc[self.dataframe_values['file_name'] == files].index[0]
            t += 1
            scene = rio.open_rasterio(files)

            if geom is not None:
                clipped_scene = read_bbox(scene, geom)
                if concatenated_array is None:
                    concatenated_array = clipped_scene
                    concatenated_array["time"]  = date
                    
                else:
                    clipped_scene["time"] = date
                    concatenated_array_xr = np.concatenate([concatenated_array, clipped_scene.values], axis=0)

            elif chunk_size is not None:
                # split array in chunks and process over chunk then recombine
                raise ValueError("Not written yet")
    
        # no data check
        if no_data is None:
        # hardcoded values based on the various available datasets
            if self.dataset == "chirps":
                no_data = [-9999]
            if self.dataset == "NDVI":
                no_data = [254, 255]
            if self.dataset == "FAPAR":
                no_data = [254, 255]
            if self.dataset == "LST":
                no_data = [0]      


        for no_data_val in no_data:
            concatenated_array_xr = np.where(concatenated_array_xr != no_data_val, concatenated_array_xr, np.nan)

        for percentile_to_process in percentiles:
            percentiles_array = np.percentile(concatenated_array_xr, percentile_to_process, axis=0)

            # check if directory for output exist, else create
            if not os.path.isdir(output_directory):
                os.makedirs(output_directory)

            clipped_scene.values[0] = percentiles_array

            clipped_scene.rio.to_raster(str(output_directory) + "/" + f"{self.dataset}_p{percentile_to_process}_{region_id}.tif")

            self.percentiles_layer[percentile_to_process] = percentiles_array


    def process_tiff_files(self, tiff_files, chunk_size):
        # Create an empty list to store the DataArrays for each file
        data_arrays = []
        
        # Open all TIFF files and concatenate into a single DataArray along the time dimension
        print("Loading the information for quantile processing...")
        for tiff_file in tqdm(tiff_files):
            # Open the TIFF file as a DataArray
            da = xr.open_rasterio(tiff_file, chunks={'x': chunk_size, 'y': chunk_size})

            # Mask values equal to -9999
            masked_data = da.where(da != -9999, drop=False)
            data_arrays.append(masked_data)
            
        # Concatenate the DataArrays along the time dimension
        ds = xr.concat(data_arrays, dim='time')
 
        ds = ds.sel(band = 1)

        ds = ds.chunk({"time": -1, "x": chunk_size, "y": chunk_size})

        # Calculate mean along the time axis, ignoring -9999 values
        mean_data = ds.mean(dim='time')
        print(mean_data.values)
        # Calculate variance along the time axis, ignoring -9999 values
        variance_data = ds.var(dim='time')

        # Calculate standard deviation along the time axis, ignoring -9999 values
        std_data = ds.std(dim='time')

        # Calculate quantiles along the time axis, ignoring -9999 values
        quantiles_data = ds.quantile([0.9, 0.5, 0.1], dim='time')


        # actual computing
        print("Mean value processing")
        mean_data.compute()
        print("Variance value processing")
        variance_data.compute()
        print("Std_dev value processing")
        std_data.compute()
        print("quantiles value processing")
        quantiles_data.compute()

        return mean_data, variance_data, std_data, quantiles_data

    
    def process_statistics(
        self,
        start_time:str = 0,
        end_time:str = 0,
        chunk_size:int = 600,
        saving_directory:os.PathLike = None,
    ) -> List[xr.DataArray]:
        """
        Process quantiles from rioxarray
        """
        complete_tiff_list  = []
        processing_status = False
        if start_time == end_time == 0:
            # entire dataset considered
            filtered_dataframe = self.dataframe_values
        else:
            filtered_dataframe = filter_dataframe(self.dataframe_values, start_time, end_time)

        for files in filtered_dataframe.iterrows():
            if "aux" in str(files[1].file_name):
                pass
            else:
                complete_tiff_list.append(files[1].file_name)

        # here a check if relevant
        data_to_save = {"mean_data", "variance_data", "std_data", "quantiles_data"}
        for name in data_to_save:
            file_to_be_checked = f"{saving_directory}/{name}.tif"
            if os.path.isfile(file_to_be_checked):
                print(f"{file_to_be_checked} already processed")
                processing_status = True
                continue
            else:
                processing_status = False
                break

        if not processing_status:
            # processing step
            mean_data, variance_data, std_data, quantiles_data = self.process_tiff_files(complete_tiff_list, chunk_size)
            
            # saving time
            if not os.path.isdir(saving_directory):
                os.makedirs(saving_directory)

            # Create a dictionary of the DataArrays you want to save
            data_arrays = {
                "mean_data": mean_data,
                "variance_data": variance_data,
                "std_data": std_data,
                "quantiles_data": quantiles_data
            }

            # Loop through the data arrays and save them as GeoTIFF files
            for name, data_array in data_arrays.items():
                file_path = f"{saving_directory}/{name}.tif"
                data_array.rio.to_raster(file_path)
                print(f"Saved {name} as GeoTIFF: {file_path}")        

            

    def process_anomaly(
        self,
        start_time:str = 0,
        end_time:str = 0,
        chunk_size:int = 600,
    ) -> xr.Dataset:
        """
        Process dekadal anomalies based on S1
        """
        # Create an empty list to store the DataArrays for each file
        data_arrays = []
        complete_tiff_list  = []
        processing_status = False
        if start_time == end_time == 0:
            # entire dataset considered
            filtered_dataframe = self.dataframe_values
        else:
            filtered_dataframe = filter_dataframe(self.dataframe_values, start_time, end_time)

        for files in filtered_dataframe.iterrows():
            if "aux" in str(files[1].file_name):
                pass
            else:
                complete_tiff_list.append(files[1].file_name)
        # 
        for tiff_file in tqdm(complete_tiff_list):
            # Open the TIFF file as a DataArray
            da = xr.open_rasterio(tiff_file, chunks={'x': chunk_size, 'y': chunk_size})

            # Mask values equal to -9999
            masked_data = da.where(da != -9999, drop=False)
            data_arrays.append(masked_data)

        # Concatenate the DataArrays along the time dimension
        ds = xr.concat(data_arrays, dim='time')
        ds = ds.sel(band = 1)
        ds = ds.chunk({"time": -1, "x": chunk_size, "y": chunk_size})

        print(ds)
        # And now resample to dekadal values and get mean
        ds_grouped = ds.resample(time="10D").mean(dim="time")
        print(ds_grouped)
        ds_dekadal = ds_grouped.mean(dim="time")
        return ds_dekadal

    

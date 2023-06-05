import xarray as xr
from shapely import geometry
from equi7grid.equi7grid import Equi7Grid
import numpy as np
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
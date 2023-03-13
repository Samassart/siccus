from siccus_reader import chirps_reader
from pathlib import Path
from datetime import datetime
import rioxarray as rio
from shapely.geometry import mapping
import geopandas
import matplotlib.pyplot as plt
if __name__ == '__main__':
    chirps_directory = Path("/data/Drysat/Datasets/SPI")
    reader_for_chirps = chirps_reader(chirps_directory)

    # read a timeseries
    #timeseries_Muanza = reader_for_chirps.read_ts(
    #    lat=-19.080808,
    #    lon=35.146659,
    #    start_time=datetime(2015, 1, 4),
    #    end_time=datetime(2019,7, 10)
    #)

    # read a regional value from shapefile
    directory_of_selected_shapefiles = Path(
        "/data/Drysat/Datasets/shapefiles/Selected_districts.shp"
        )
    districts_shapefile = geopandas.read_file(
        directory_of_selected_shapefiles,
        crs="epsg:4326"
        )
    
    # case study of Mabote
    Mabote_geom = districts_shapefile[
        districts_shapefile.NAME_2 == "Massinga"
        ].geometry
    
    mabote_xarray_20170515 = reader_for_chirps.read_bbox(
        datetime_to_match=datetime(2017, 5, 15),
        geom=Mabote_geom
    )


    # test for NDVI
    ndvi_directory = Path("/data/Drysat/Datasets/NDVI")
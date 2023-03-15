# %%
from siccus_reader import chirps_reader, cgls_reader
from aesthetic_function import environmental_timeseries
from pathlib import Path
from datetime import datetime
import geopandas


# %% set paths
rainf_directory = Path("/data/Drysat/Datasets/SPI")
veget_directory = Path("/data/Drysat/Datasets/FAPAR")

# read a regional value from shapefile
directory_of_selected_shapefiles = Path(
    "/data/Drysat/Datasets/shapefiles/Selected_districts.shp"
    )

# %%
districts_shapefile = geopandas.read_file(
    directory_of_selected_shapefiles,
    crs="epsg:4326"
    )
# case study of Mabote
Mabote_geom = districts_shapefile[
    districts_shapefile.NAME_2 == "Massinga"
    ].geometry

# %% initialize the reader for chirps
reader_for_chirps = chirps_reader(rainf_directory)

# %% read a timeseries
timeseries_Muanza_chirps = reader_for_chirps.read_ts(
    lat=-19.10,
    lon=35.15,
    start_time=datetime(2001, 1, 4),
    end_time=datetime(2022,7, 10)
    )

# %% read a bounding box from Mabote
mabote_xarray_20170515 = reader_for_chirps.read_bbox(
    datetime_to_match=datetime(2017, 5, 15),
    geom=Mabote_geom
    )

# %% initialize the reader
reader_for_veget = cgls_reader(directory=veget_directory, dataset="FAPAR")

# %% read a timeseries
timeseries_Muanza_ndvi = reader_for_veget.read_ts(
    lat=-19.080808,
    lon=35.146659,
    start_time=datetime(2015, 1, 4),
    end_time=datetime(2019,7, 10)
    )

# %% same bbox approach for fapar/ndvi
mabote_xarray_20170515 = reader_for_veget.read_bbox(
    datetime_to_match=datetime(2017, 5, 15),
    geom=Mabote_geom
    )

# %% make aesthetically pleasing timeserie
g, axe = environmental_timeseries(
    timeseries_Muanza_chirps,
    "chirps"
)

g, axe = environmental_timeseries(
    timeseries_Muanza_ndvi,
    "FAPAR"
)
# %%

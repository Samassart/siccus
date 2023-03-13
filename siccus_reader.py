import os
from datetime import datetime
import pandas as pd
from temporal_function import (
    from_day_to_dekad,
    from_dekad_to_day,
    filter_dataframe,
    datetime_to_year_month_dekad
    )
import rioxarray as rio
import xarray as xr
from shapely import geometry
import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import mapping
import zipfile

class chirps_reader:
    """
    Class for reading and plotting chirps dataset
    """
    def __init__(
            self,
            directory:os.PathLike
    ) -> None:
        self.directory: directory
        # initialize datacube
        self.dataframe_values = {}
        for files in sorted(directory.glob("chirps*")):
            files_data = str(files.stem).split(".")
            datetime_string = "".join(files_data[2:4])
            dekad_string = files_data[-1]
            full_filename = files
            self.dataframe_values[int(datetime_string+dekad_string)] = (full_filename)
        self.dataframe_values = pd.DataFrame(
            {"file_name":self.dataframe_values.values()},
            index=self.dataframe_values.keys()
            )

    def read_ts(
            self,
            lat:int,
            lon:int,
            start_time:datetime,
            end_time:datetime
    ) -> pd.DataFrame:
        """
        Read and store a timeseries of the chirps dataset at a given lat/lon
        for given years/months/dekads period
        """
        # first convert the day into dekad*
        thresh_start = datetime_to_year_month_dekad(start_time)
        thresh_end = datetime_to_year_month_dekad(end_time)

        # filter the dataset 
        filtered_dataset = filter_dataframe(
                self.dataframe_values,
                thresh_start,
                thresh_end
        )

        timeseries_to_extract = {}
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
            file_to_extract_chirps = rio.open_rasterio(file[1].values[0])
            specific_pixel_chirps = file_to_extract_chirps.sel(
                x=lon,
                y=lat,
                method="nearest"
                ).values[0]

            # from dict to pandas    
            timeseries_to_extract[file_datetime_date] = specific_pixel_chirps
        timeseries_to_extract_pd = pd.DataFrame(
            {"Chirps":timeseries_to_extract.values()},
            index=timeseries_to_extract.keys()
            )

        return timeseries_to_extract_pd
    

    def read_bbox(
            self,
            datetime_to_match:datetime,
            geom:geometry.multipolygon.MultiPolygon | geometry.polygon.Polygon
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
        self.dataframe_values_dt = self.dataframe_values
        self.dataframe_values_dt.index = datetime_index_dt

        #
        idx = self.dataframe_values_dt.index[
            self.dataframe_values_dt.index.get_indexer(
                [datetime_to_match], method='nearest'
                    )
                ]

        # read the .tif file
        file_to_read = self.dataframe_values_dt.loc[idx]
        chirps_scene = rio.open_rasterio(file_to_read.file_name.values[0])
        chirps_scene.rio.set_spatial_dims(x_dim="x", y_dim="y", inplace=True)
        chirps_scene.rio.write_crs("epsg:4326", inplace=True)
        # crop with the (multi)polygon
        clipped_chirps_scene = chirps_scene.rio.clip(geom.values)
        
        return clipped_chirps_scene
    
class cgls_reader:
    """
    reader for various cgls based dataset
    Set the relevant dataset that you want
    Currently supported
    """
    def __init__(
            self,
            directory:os.PathLike,
            dataset:str
    ) -> None:
        self.directory: directory
        self.dataset = dataset
        self.dataframe_values = {}
        for file_dir in sorted(directory.glob(f"{self.dataset}*")):
            # cgls files are typically given in zip format when
            # downloaded. This loop checks and unzip the files if so.
            for files in file_dir.glob("*.zip"):
                output_file = str(files.parent / files.stem) 
                with zipfile.ZipFile(files,"r") as zip_ref:
                    zip_ref.extractall(output_file)
            for files in file_dir.glob("*"):    
                files_data = str(files.stem).split("_")
                # boilerplate test 
                if not files_data[2] == self.dataset:
                    raise NameError("Wrong dataset - check input and filename")
                
                datetime_string = files_data[3]
                version_dataset = files_data[-1]
                # build the string similar to chirps
                first_part_datetime_string = datetime_string[0:6]
                second_part_datetime_string =  str(int(datetime_string[6]) + 1)
                new_datetime_string = (
                    first_part_datetime_string
                    + second_part_datetime_string
                )
                self.dataframe_values[int(new_datetime_string)] = (files)

        self.dataframe_values = pd.DataFrame(
            {"file_name":self.dataframe_values.values()},
            index=self.dataframe_values.keys()
            )

    def read_ts(
            self,
            lat:int,
            lon:int,
            start_time:datetime,
            end_time:datetime
    ) -> pd.DataFrame:
        """
        Read and store a timeseries of a cgls dataset at a given lat/lon
        for given years/months/dekads period
        """
        # first convert the day into dekad*
        thresh_start = datetime_to_year_month_dekad(start_time)
        thresh_end = datetime_to_year_month_dekad(end_time)

        # filter the dataframe
        filtered_dataset = filter_dataframe(
                self.dataframe_values,
                thresh_start,
                thresh_end
            )

"""        # run the dataframe
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
            file_to_extract_chirps = rio.open_rasterio(file[1].values[0])
            specific_pixel_chirps = file_to_extract_chirps.sel(
                x=lon,
                y=lat,
                method="nearest"
                ).values[0]"""
from reader_functions import environmental_dataset
import pandas as pd
import os
import zipfile
from pathlib import Path
from temporal_function import (
    from_year_and_doy_to_datetime,
    from_day_to_dekad,
    )
from datetime import datetime


class modis_reader(environmental_dataset):
    """
    Class to read LST from the modis data
    """
    def __init__(
            self,
            directory:os.PathLike,
            dataset:str
    ) -> None:
        self.directory=directory
        self.dataset=dataset

        # initialize datacube
        self.dataframe_values = {}

        for files in sorted(directory.glob("MOD*")):
            files_data = str(files.stem).split("_")
            datetime_string = files_data[4]

            # .tif.aux boilerplate
            if str(datetime_string[-8:]) == ".tif.aux":
                continue
            
            # decrypt data
            year_value = int(datetime_string[3:7])
            doy_value = int(datetime_string[7:])
            datetime_file = from_year_and_doy_to_datetime(year_value, doy_value)
            dekad_of_file = str(from_day_to_dekad(datetime_file.day))

            # create the index
            datetime_index = datetime.strftime(datetime_file, "%Y%m") + dekad_of_file

            # append the index and filename
            self.dataframe_values[int(datetime_index)] = files

        # transform into dataframe
        self.dataframe_values = pd.DataFrame(
            {"file_name":self.dataframe_values.values()},
            index=self.dataframe_values.keys()
            )

class chirps_reader(environmental_dataset):
    """
    Class for reading and plotting chirps dataset
    """
    def __init__(
            self,
            directory:os.PathLike,
            dataset:str
    ) -> None:
        self.directory=directory
        self.dataset=dataset

        # initialize datacube
        self.dataframe_values = {}
        for files in sorted(directory.glob("chirps*")):
            files_data = str(files.stem).split(".")
            # aux filter 
            if files_data[-1] == "aux":
                continue
            datetime_string = "".join(files_data[2:4])
            dekad_string = files_data[-1]
            full_filename = files
            self.dataframe_values[int(datetime_string+dekad_string)] = (full_filename)
        self.dataframe_values = pd.DataFrame(
            {"file_name":self.dataframe_values.values()},
            index=self.dataframe_values.keys()
            )

class cgls_reader(environmental_dataset):
    """
    reader for various cgls based dataset
    Set the relevant dataset that you want
    Currently supported:
    NDVI
    FAPAR
    CHIRPS
    LST
    """
    def __init__(
            self,
            directory:os.PathLike,
            dataset:str 
    ) -> None:
        self.directory=directory
        self.dataset = dataset
        self.dataframe_values = {}
        for file_dir in sorted(directory.glob(f"{self.dataset}*")):
            # cgls files are typically given in zip format when
            # downloaded. This loop checks and unzip the files if so.
            for files in file_dir.glob("*"):
                output_file = str(files.parent / files.stem) 
                if ".zip" in str(files):
                    # for NDVI case 
                    with zipfile.ZipFile(files,"r") as zip_ref:
                        zip_ref.extractall(output_file)
                        # when extracted it can be deleted 
                        if list(Path(output_file).glob("*")) != 0:
                            os.remove(files) 
                else:
                    pass

            # locate the file in the paths
            for ndvi_folder in Path(file_dir).glob("*"):
                for date_folder in ndvi_folder.glob("*"):
                    for files in date_folder.glob(
                        f"*{self.dataset}-{self.dataset}*"
                        ):
                        files_data = str(files.stem).split("_")
                        # boilerplate tests
                        if len(files_data[2].split("-")) > 2:
                            continue
                        if not (files_data[2].split("-")[0] == self.dataset):
                            raise NameError("Wrong dataset - check input and filename")
      
                        datetime_string = files_data[3]
                        version_dataset = files_data[-1]
                        # build the string similar to chirps
                        if version_dataset.split(".")[-1] == "aux":
                            continue
                        first_part_datetime_string = datetime_string[0:6]

                        if self.dataset=="NDVI":
                            second_part_datetime_string =  str(int(datetime_string[6]) + 1)
                        else:
                            second_part_datetime_string =  str(int(datetime_string[6]))
                        new_datetime_string = (
                            first_part_datetime_string
                            + second_part_datetime_string
                        )
                        self.dataframe_values[int(new_datetime_string)] = (files)

                    # additional step for fapar
                    for version in ["RT0", "RT1", "RT2", "RT6"]:
                        for files in date_folder.glob(
                            f"*{version}-{self.dataset}*"
                            ):
                            files_data = str(files.stem).split("_")
                            # boilerplate test 
                            if not files_data[2].split("-")[0] == self.dataset:
                                raise NameError("Wrong dataset - check input and filename")
                            datetime_string = files_data[3]
                            version_dataset = files_data[-1]
                            # build the string similar to chirps
                            first_part_datetime_string = datetime_string[0:6]
                            if self.dataset=="NDVI":
                                second_part_datetime_string =  str(int(datetime_string[6]) + 1)
                            else:
                                second_part_datetime_string =  str(int(datetime_string[6]))
                            new_datetime_string = (
                                first_part_datetime_string
                                + second_part_datetime_string
                            )
                            self.dataframe_values[int(new_datetime_string)] = (files)

        self.dataframe_values = pd.DataFrame(
            {"file_name":self.dataframe_values.values()},
            index=self.dataframe_values.keys()
            ).sort_index()



import os
from datetime import datetime
import pandas as pd
from temporal_function import (
    from_day_to_dekad,
    from_dekad_to_day,
    filter_dataframe,
    datetime_to_year_month_dekad
    )
from spatial_function import read_bbox
from tqdm import tqdm
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

    def process_percentiles_chunk(
        self,
        file_type:str,
        output_directory:os.PathLike,
        geom:geometry.multipolygon.MultiPolygon | geometry.polygon.Polygon = None,
        chunk_size:int = None,
        percentiles:list[float] = [], 
        no_data:list = None,
        region_id:str = "full_data",
        ):
        if file_type == "tif":
            import rasterio
            percentiles = [10, 25, 50, 75, 90]

            # define the chunk size for reading the data from the raster layer
            chunk_size = 256

            # loop over each .tif layer in the list
            for layer_path in tqdm(self.dataframe_values["file_name"]):
                # open the .tif layer using rasterio
                with rasterio.open(layer_path) as src:
                    # calculate the total number of chunks for the raster layer
                    num_chunks = (src.height // chunk_size + 1) * (src.width // chunk_size + 1)
                    # initialize a list to store the percentile values for each chunk
                    percentile_list = []
                    # loop over each chunk of the raster layer
                    for i, window in tqdm(enumerate(src.block_windows(window_size=chunk_size, nodata=src.nodata)), total=num_chunks):
                        # get the window indices
                        row_start, col_start = window.ul[0], window.ul[1]
                        row_stop, col_stop = window.lr[0], window.lr[1]
                        # read the pixel values for the current chunk into a numpy array
                        pixel_array = src.read(window=window)
                        # reshape the numpy array to have two dimensions: pixel values and temporal values
                        pixel_array = pixel_array.reshape(pixel_array.shape[1], -1)
                        # calculate the percentiles along the second dimension
                        percentile_values = np.percentile(pixel_array, percentiles, axis=1)
                        # append the percentile values to the list for this layer
                        percentile_list.append(percentile_values)
                    # concatenate the percentile values for all chunks into a single numpy array
                    percentile_values = np.concatenate(percentile_list, axis=1)
                    # store the percentile values in a list or numpy array
                    # you can use the layer_path as the key if you want to keep track of which layer each percentile value corresponds to
                    # e.g., percentile_dict[layer_path] = percentile_values
                    # or, if you just want a list of percentile values:
                    print(percentile_values)
    def compare_scene_to_percentiles_map(
            self,
            map_to_categorize:xr.DataArray,
            percentile_directories:os.PathLike,
            percentiles_to_use:list[float | int], # put in order
            region_id:str = "full_data",
            geom:geometry.multipolygon.MultiPolygon | geometry.polygon.Polygon = None,
            cat_val:int = 1
        ) -> xr.DataArray:
        """
        Take a scene, and create a categorical drought map
        """
        # sort the percentile, key to the process 
        percentiles_to_use.sort(reverse = True)

        # initialise multicategorical drought map variable and dict of category
        categorical_drought_map = None
        category_for_plot = {}
        # first check if the percentile map exist, if not, call the percentile function first
        for idx, percentile_analysed in enumerate(percentiles_to_use):
            # take care of the category dictionary
            category_for_plot[idx+1] = percentile_analysed

            if os.path.isdir(percentile_directories):
                file_for_perc = str(percentile_directories) + "/" + f"{self.dataset}_p{percentile_analysed}_{region_id}.tif"
                if os.path.isfile(file_for_perc):
                    percentile_scene = rio.open_rasterio(file_for_perc)
                else:
                    print(f"Percentiles {percentile_analysed} not found, start the processing")
                    if geom is None:
                        raise ValueError("Percentiles not processed, please add the geometry for the processing.")
                    else:
                        self.process_percentiles(
                            output_directory=percentile_directories,
                            geom=geom,
                            percentiles=[percentile_analysed],
                            no_data=[255, 241],
                            region_id = region_id
                            )

                    percentile_scene = rio.open_rasterio(file_for_perc)

            # Compare the scenes
            if categorical_drought_map is None:
                categorical_drought_map = xr.where(map_to_categorize < percentile_scene, 1, 0)
            else:
                categorical_drought_map_to_add = xr.where(map_to_categorize < percentile_scene, 1, 0)
                categorical_drought_map += categorical_drought_map_to_add

        # and finally 0 to nan
        categorical_drought_map = categorical_drought_map.where(categorical_drought_map > 0.1)
        return categorical_drought_map, category_for_plot
    

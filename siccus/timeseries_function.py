import seaborn as sns
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.dates as md
import pandas as pd
from cmcrameri import cm
from pytesmo.time_series import anomaly


def anomaly_plotting(
        variable_dataset:pd.DataFrame, 
        window_size:int = 30) -> pd.DataFrame:
    """
    Plot soil moisture anomalies from a pandas dataframe
    """
    anomaly_dataframe = anomaly.calc_anomaly(
        variable_dataset[variable_dataset.columns[0]],
        window_size=window_size
        ).to_frame(name=variable_dataset.columns[0])
    return anomaly_dataframe

def annual_plot(
    variable_dataset:pd.DataFrame,
    ):
    # Split the dataframe into different years
    year_groups = variable_dataset.groupby(pd.Grouper(freq='Y'))

    # Create the plot
    fig, ax = plt.subplots(figsize=(8, 5))
    sns.set_style("whitegrid")

    # Plot each year's data in grey, except for the last year in blue
    for year, data in year_groups:
        if year == list(year_groups.groups.keys())[-2]:
            sns.lineplot(x=data.index.dayofyear, y=data[variable_dataset.columns[0]], color='blue', label=str(year), ax=ax)
        else:
            sns.lineplot(x=data.index.dayofyear, y=data[variable_dataset.columns[0]], color='grey', alpha=0.5, ax=ax)

    # Add labels and legend
    ax.set_xlabel('Day of the Year')
    ax.set_ylabel('Value')
    ax.set_title(variable_dataset.columns[0])
    ax.legend(title='Year')
    plt.show()


def environmental_timeseries(
        dataframe_dict:list,
        percentile_plotting:list = None,
        z_score_plotting:bool = False,
        savepath:os.PathLike = None,
        mono_map:bool = True
):
    """
    Get a bunch of timeseries dataframe and 
    plot them in a clean way
    """
    nb_lines_to_plot = len(dataframe_dict)
    variable_names = []
    for individual_dataframe in dataframe_dict:
        variable_names.append(individual_dataframe.columns[0])

    # set plot characteristics
    list_of_colours = [
        "#B6443F", # redish
        "#367588", # azure
        "#00A86B", # greenish
        "#FA7A35", # orangish
        "#663399", # purple
    ]
    palette_default = sns.color_palette(list_of_colours, 5)
    sns.set_theme(
        style="whitegrid",
        rc = {'lines.linewidth': 0.59},),
    sns.set_palette(palette_default)
    g, axe = plt.subplots(figsize = (20, 5))
    min_extent_x = None
    max_extent_x = None
    min_extent_y = None
    max_extent_y = None
    threshold_val = None
    # iterate through the data to plot
    dataset_name_list = []
    for idx, dataframe in enumerate(dataframe_dict):
        # get var name
        dataframe_name = dataframe.columns[0]

        # get Z-score
        if z_score_plotting:
            Mean_for_z=dataframe[f"{dataframe_name}"].mean()
            Std_for_z=dataframe[f"{dataframe_name}"].std()

            Z_score = (dataframe[f"{dataframe_name}"] - Mean_for_z) / Std_for_z
            dataframe["Z_score"] = Z_score

            # prepare the mask
            mask = dataframe["Z_score"] < -1
            dataframe["masked_values"] = dataframe[dataframe_name].mask(mask)
            dataframe["z_score_hue"] = mask
            # 


        dataset_name_list.append(dataframe_name)
        dataframe["hue"] = idx
        if not z_score_plotting:
            g = sns.lineplot(
                data=dataframe,
                x=dataframe.index,
                y=dataframe_name,
                linewidth=0.95,
                color=palette_default[idx],
                )
        if z_score_plotting:
            print("To implement")
            break

        # find extent of timeseries
        min_extent_x_to_test = dataframe.index[0]
        max_extent_x_to_test = dataframe.index[-1]
        min_extent_y_to_test = np.nanmin(dataframe[f"{dataframe_name}"])
        max_extent_y_to_test = np.nanmax(dataframe[f"{dataframe_name}"])

        # x-axis extent
        if min_extent_x is None and max_extent_x is None:
            min_extent_x = min_extent_x_to_test
            max_extent_x = max_extent_x_to_test

        else:
            # if longer extent, replace
            if min_extent_x > min_extent_x_to_test:
                min_extent_x = min_extent_x_to_test
            if max_extent_x < max_extent_x_to_test:
                max_extent_x = max_extent_x_to_test

        # y-axis extent
        if min_extent_y is None and max_extent_y is None:
            min_extent_y = min_extent_y_to_test
            max_extent_y = max_extent_y_to_test
        else:
            # if longer extent, replace
            if min_extent_y < min_extent_y_to_test:
                min_extent_y = min_extent_y_to_test
            if max_extent_y > max_extent_y_to_test:
                max_extent_y = max_extent_y_to_test

        # clean percentile limit thresholding
        high_percentile = np.nanpercentile(dataframe[f"{dataframe_name}"], 95)
        low_percentile = np.nanpercentile(dataframe[f"{dataframe_name}"], 5)

        threshold_percentile_high = max_extent_y - high_percentile
        threshold_percentile_low = low_percentile - min_extent_y

        max_y_value = max_extent_y + threshold_percentile_high
        min_y_value = min_extent_y - threshold_percentile_low


        if threshold_val is None:
            definite_max_y_value = max_y_value
            definite_min_y_value = min_y_value
            threshold_val = "already one dataset"
        else: 
            if max_y_value > definite_max_y_value:
                definite_max_y_value = max_y_value
            if min_y_value < definite_min_y_value:
                definite_min_y_value = min_y_value



    # check impossible combination of param
    if (nb_lines_to_plot>1) and (not percentile_plotting is None):
        raise ValueError("Cannot plot multiple timeseries with percentiles")
    else:
        # calculate the dataframe percentiles
        percentile_to_plot = {}
        for percentiles in percentile_plotting:
            percentile_to_plot[percentiles] = np.nanpercentile(dataframe[f"{dataframe_name}"], percentiles)

    try:
        axe.set_ylim(   
            definite_min_y_value, 
            definite_max_y_value,
            ) 
        axe.set_xlim(
            min_extent_x,
            max_extent_x
            )
    except:
        pass
    
    plt.xlabel('Date', fontsize=20)
    plt.ylabel(f'{dataframe_name}', fontsize = 20)

    colourmap = cm.lajolla_r
    len_of_pc_list = len(percentile_to_plot)
    incremental_increase_cmap = np.floor((len(colourmap.colors) / len_of_pc_list))
    for index, percentile_line in enumerate(percentile_to_plot):

        index_of_color = int(((index+1) * incremental_increase_cmap) - 1)
        colour_line = colourmap.colors[index_of_color]
        axe.axhline(y =percentile_to_plot[percentile_line], linestyle = "-", color=colour_line, linewidth= 1.8)


    # legend data management
    leg = plt.legend(
        dataset_name_list,
        loc='upper right',
        bbox_to_anchor=(1, 1.13),
        fontsize=14
        )
    for idx, dataframe in enumerate(dataframe_dict):
        leg.legendHandles[idx].set_color(palette_default[idx])

    # set ticks param
    # Set major and minor date tick locators
    maj_loc = md.MonthLocator(bymonth=np.arange(1,12,6))
    axe.xaxis.set_major_locator(maj_loc)
    min_loc = md.MonthLocator()
    axe.xaxis.set_minor_locator(min_loc)
    # Set major date tick formatter
    zfmts = ['', '%b\n%Y', '%b', '%b-%d', '%H:%M', '%H:%M']
    maj_fmt = md.ConciseDateFormatter(maj_loc, zero_formats=zfmts, show_offset=False)
    axe.xaxis.set_major_formatter(maj_fmt)
    axe.figure.autofmt_xdate(rotation=45, ha='center')

    props = dict(boxstyle='round', facecolor='white', alpha=0.45,edgecolor="C7")
       
    if mono_map:
        return g, axe
    else:
        return axe

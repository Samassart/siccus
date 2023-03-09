import pandas as pd

def from_day_to_dekad(day:int) -> int:
    # dekad analysis
    dekad_to_day_dict = {
        1:[*range(1, 11, 1)],
        2:[*range(11, 21, 1)],
        3:[*range(21, 32, 1)]
        }
    
    # test each three potential dekad
    the_right_dekad = False
    dekad_to_test = 1
    while not the_right_dekad:
        if not day in dekad_to_day_dict[dekad_to_test]:
            dekad_to_test += 1
        elif day in dekad_to_day_dict[dekad_to_test]:
            return dekad_to_test
        
        if dekad_to_test > 10:
            raise ValueError("Bug")


def filter_dataframe(
        dataframe:pd.DataFrame,
        start_threshold:str,
        end_threshold:str
    ) -> pd.DataFrame:
    # filter the dekad
    half_filtered_dataframe = dataframe[
        dataframe.index > int(start_threshold)
        ]
    filtered_dataframe = half_filtered_dataframe[
        half_filtered_dataframe.index < int(end_threshold)
        ]

    return filtered_dataframe

def from_dekad_to_day(dekad:int) ->int:
    if dekad==1:
        day = 5
    elif dekad==2:
        day = 15
    elif dekad==3:
        day = 25
    else:
        raise ValueError("Dekad is likely wrong")
    return day

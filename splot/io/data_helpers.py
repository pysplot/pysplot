import numpy as np
import pandas as pd


def _transform_to_skycoord(dataframe, frame, units, coord_kwargs={}):
    """
    Converts data to astropy SkyCoord object of sunpy frame with astropy.units.
    Determines xyz or lat/lon by number of columns in shape.
    """
    data = dataframe.values
    obstimes = dataframe.index
    coord = None
    if not isinstance(units, list):
        units = [units] * data.shape[1]
    units = [v if v is not None else 1 for v in units]

    if frame is not None:
        if data.shape[1]==3:
            coord = frame(data[:,0]*units[0], data[:,1]*units[1], data[:,2]*units[2], obstime=obstimes, **coord_kwargs)
        elif data.shape[1]==2:
            coord = frame(data[:,0]*units[0], data[:,1]*units[1], obstime=obstimes, **coord_kwargs)
    return coord

def _transform_into_quantity(dataframe, units):
    data = dataframe.values
    quantity = None
    if units is not None:
        quantity = data * units
    return quantity


def _convert_coordinates(location, desired_coord=None, coord_kwargs={}, keep_column_names=False):
    if (location.skycoord is not None) and (desired_coord is not None) and (location.coord != desired_coord):
        transformed_skycoord = location.skycoord.transform_to(desired_coord(obstime=location.data.index, **coord_kwargs))
        location.converted_skycoord = transformed_skycoord
        location.coord = desired_coord
        skycoord_cols = transformed_skycoord.get_representation_component_names().keys()#transformed_skycoord._data.__dict__.keys()
        for col,scol in zip(location.columns, skycoord_cols):
            #try:
            vals = transformed_skycoord.__getattr__(scol).value#_data.__dict__[scol].value
            #except Exception:
            #    vals = transformed_skycoord.__getattr__(scol)
            location.data[col] = vals
        if not keep_column_names:
            location.data = location.data.rename(columns={c: sc for c,sc in zip(location.columns, skycoord_cols)})
            location.columns = skycoord_cols        
    return location

def _convert_units(science, desired_units=None):
    if (science.quantity is not None) and (desired_units is not None) and (science.units != desired_units):
        transformed_quantity = science.quantity.to(desired_units).value
        science.converted_quantity = transformed_quantity
        science.units = desired_units
        for i,col in enumerate(science.columns):
            science.data[col] = transformed_quantity[:,i]
    return science


def _match_data_cadence(location, science, data_cadence=None, interpolation_method='time', combine_method='mean'):
    """
    match cadence of location and science objects
    """
    ## set combine method function
    combine_method_dict = {
        'mean': np.nanmean, 
        'median': np.nanmedian,
    }

    if combine_method in combine_method_dict:
        combine_method_function = combine_method_dict[combine_method]
    #else:
        #TODO: return and throw error

    location_data = location.data


    if science is not None:
        ## TODO: change units
        #if desired_units is not None:
        #    if science_data.units is not None and science_data.units!=desired_units:
        #        science_data.convert to desired_units
        science_data = science.data
    
    merged_data = location_data.join(science_data, how='outer')

    if not data_cadence:
        for col in location.columns:
            merged_data[col] = merged_data[col].interpolate(interpolation_method)
        merged_data = merged_data.reindex(science_data.index)
    else:
        merged_data = merged_data.resample(data_cadence).apply(combine_method_function)
    return merged_data

def _transform_into_dataframe(list_of_data, column_names=[], combine_axis='columns'):
    """
    change into pandas dataframe, set time as index
    """
    combined_dataframe = pd.DataFrame()
    combined_columns = []
    for data, columns in zip(list_of_data, column_names):
        dataframe = pd.DataFrame()
        dataframe.index = pd.to_datetime(data['x'])
        if len(data['y'].shape)==1:
            ncol=1
        else:
            ncol = data['y'].shape[1]
        ##TODO: check column and data length match.
        if ncol==1:
            dataframe[columns[0]] = data['y']
        else:
            for c,colname in enumerate(columns):
                dataframe[colname] = data['y'][:,c]
        combined_dataframe = pd.concat([combined_dataframe, dataframe], axis=combine_axis)
        for col in columns:
            combined_columns.append(col)
    return combined_dataframe, combined_columns
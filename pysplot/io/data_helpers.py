import inspect
import numpy as np
import pandas as pd
import astropy.coordinates
import sunpy.coordinates
import astropy.units

def _transform_into_skycoord(dataframe, frame, units, coord_kwargs={}):
    """
    Convert the input data from a Pandas DataFrame to an Astropy SkyCoord object
    of a specified coordinate frame with appropriate units.

    Parameters
    ----------
    dataframe : pandas.DataFrame
        A DataFrame containing coordinate data, where each row represents a point and
        each column corresponds to a coordinate (e.g., x, y, z or latitude, longitude).
        The DataFrame's index should represent observation times.
    
    frame : sunpy.coordinates.frame or astropy.coordinates.frame
        The coordinate frame class to use for creating the SkyCoord object. Typically a SunPy or Astropy frame.
    
    units : astropy.units
        The units for the coordinates in the DataFrame. If a single unit object is provided,
        it will be applied to all coordinates. If a list is provided, it should have
        the same length as the number of columns in `dataframe`.
    
    coord_kwargs : dict, optional
        Additional keyword arguments to pass to the frame constructor, such as metadata
        or other frame-specific parameters (default is an empty dictionary).

    Returns
    -------
    coord : astropy.coordinates.SkyCoord
        An Astropy SkyCoord object created using the input data, units, and frame.
        The object represents the coordinates in the specified frame, with the 
        appropriate unit conversion applied.
    
    Notes
    -----
    The function assumes that the number of columns in `dataframe` corresponds to
    the dimensions of the coordinate system:
        - 3 columns -> xyz coordinates
        - 2 columns -> latitude/longitude or similar.
    If the unit for a coordinate is not provided (None), it will default to unit 1.
    
    Example
    -------
    >>> from astropy.coordinates import GCRS
    >>> from astropy.units import AU
    >>> import pandas as pd
    >>> df = pd.DataFrame([[1, 2, 3], [4, 5, 6]], columns=['x', 'y', 'z'], index=['t1', 't2'])
    >>> coord = _transform_to_skycoord(df, GCRS, u.AU)
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
    """
    Convert the values in a Pandas DataFrame to an Astropy Quantity object by applying the specified units.

    Parameters
    ----------
    dataframe : pandas.DataFrame
        A DataFrame containing numerical data to be converted into an Astropy Quantity. Each entry
        in the DataFrame will be multiplied by the specified units.
    
    units : str or `astropy.units.Unit` or `None`
        The units to apply to the data. If not `None`, each value in the DataFrame will be multiplied
        by this unit. If `None`, the data will not be converted and the function returns `None`.
    
    Returns
    -------
    quantity : `astropy.units.Quantity` or `None`
        An Astropy `Quantity` object with the same shape as the input DataFrame. Each entry in the 
        DataFrame is multiplied by the specified units. If `units` is `None`, the function returns `None`.
    
    Notes
    -----
    This function assumes that the provided `units` are compatible with the data in the DataFrame.
    
    Example
    -------
    >>> import pandas as pd
    >>> from astropy import units as u
    >>> df = pd.DataFrame([[1, 2], [3, 4]], columns=['a', 'b'])
    >>> result = _transform_into_quantity(df, u.meter)
    >>> print(result)
    """
    data = dataframe.values
    quantity = None
    if units is not None:
        quantity = data * units
    return quantity


def _convert_coordinates(spatial, desired_coord=None, coord_kwargs={}):
    """
    Convert the coordinates of a given location to a desired coordinate frame, if applicable.

    Parameters
    ----------
    spatial : pysplot.SpatialData object
        An object that contains `data` (a DataFrame) and a `skycoord` (an Astropy `SkyCoord` object).
        The `data` represents the coordinates, and the `skycoord` contains the original coordinates 
        in a particular frame.

    desired_coord : astropy.coordinates.frame, optional
        The desired coordinate frame to transform the `SkyCoord` object to. This should be a valid frame 
        class, such as `GCRS`, `ICRS`, `HeliographicStonyhurst`, etc. If `None`, no transformation is applied.

    coord_kwargs : dict, optional
        Additional keyword arguments to pass to the `desired_coord` frame constructor (default is an empty dictionary).

    Returns
    -------
    transformed_df : pandas.DataFrame
        A DataFrame containing the transformed coordinates in the desired coordinate frame. The column names 
        will be updated according to the transformed frame's components.

    transformed_skycoord : astropy.coordinates.SkyCoord or None
        The transformed `SkyCoord` object corresponding to the new coordinate frame. If no transformation 
        is needed, this will be the original `SkyCoord` object from `location`.

    Notes
    -----
    - The transformation occurs only if the `location.skycoord` is not `None`, `desired_coord` is provided,
      and the current coordinate frame (`location.coord`) differs from the target (`desired_coord`).
    - The function assumes that the `location` object contains valid data with corresponding `SkyCoord` attributes.
    - If no transformation is required, the original data and `SkyCoord` are returned as-is.

    Example
    -------
    >>> from sunpy.coordinates import HeliocentricInertial, GeocentricSolarMagnetospheric
    >>> location = SpatialData(location, coord=HeliocentricInertial)
    >>> transformed_df, transformed_skycoord = _convert_coordinates(location, GeocentricSolarMagnetospheric)
    """
    transformed_df = pd.DataFrame()
    transformed_df.index = spatial.data.index
    if (spatial.skycoord is not None) and (desired_coord is not None) and (spatial.coord != desired_coord):
        transformed_skycoord = spatial.skycoord.transform_to(desired_coord(obstime=spatial.data.index, **coord_kwargs))
        ## TODO: have choice of what type of representation for output (cartesian, cyclindrical, etc)
        skycoord_cols = transformed_skycoord.get_representation_component_names().keys()
        for col,scol in zip(spatial.columns, skycoord_cols):
            vals = transformed_skycoord.__getattr__(scol).value
            transformed_df[col] = vals
        transformed_df.rename(columns={c: sc for c,sc in zip(spatial.columns, skycoord_cols)}, inplace=True)
    else:
        transformed_df = spatial.data
        transformed_skycoord = spatial.skycoord
    return transformed_df, transformed_skycoord


def _convert_units(science, desired_units=None):
    """
    Converts the units of a `science.quantity` object to the desired units.

    This function takes a `science` object (which is expected to be a 
    `science.quantity` type) and transforms its values into the specified 
    `desired_units`. If no desired units are provided, it defaults to the 
    current units of the `science` object.

    Parameters:
    -----------
    science : object
        An object containing a `quantity` attribute (e.g., astropy.quantity) 
        with data and units. The object should also have a `columns` attribute 
        and a `data` attribute, which holds the raw data values in a 
        `pandas.DataFrame`.
        
    desired_units : str, optional
        The units to which the `science.quantity` should be converted. If 
        `None`, the current units of the `science` object will be used.

    Returns:
    --------
    transformed_df : pandas.DataFrame
        A DataFrame with the same index as the original `science.data`, 
        containing the transformed values in the desired units. If no units 
        transformation is needed, the original data is returned unchanged.

    Notes:
    ------
    - If the `science` object does not have units (`science.units is None`), 
      the original data is returned without any conversion.
    - If the units of the `science.quantity` are already equal to the desired 
      units, no transformation will occur.
    """
    transformed_df = pd.DataFrame()
    transformed_df.index = science.data.index
    if science.units is not None:
        if (science.units != desired_units):
            output_units = desired_units
        else:
            output_units = science.units
        transformed_quantity = science.quantity.to(output_units).value
        for i,col in enumerate(science.columns):
            transformed_df[col] = transformed_quantity[:,i]
    else:
        transformed_df = science.data
    return transformed_df


def _match_data_cadence(location_data, science_data, data_cadence=None, interpolation_method='time', combine_method='mean'):
    """
    Match the cadence (time intervals) of `location_data` and `science_data` by resampling or interpolating 
    the data as necessary matching their time intervals. 
    If the cadences are mismatched, it either resamples or interpolates the `location_data` to match the 
    cadence of the `science_data`, or applies a specified combine method to downsample the data. 

    Parameters
    ----------
    location_data : pandas.DataFrame
        A DataFrame containing the location data, with its index representing the time or cadence of the data.
    
    science_data : pandas.DataFrame, optional
        A DataFrame containing the science data to be aligned with `location_data`. Its index represents 
        the time or cadence for the science data. If `None`, the function will only adjust `location_data`.

    data_cadence : str, optional
        The desired cadence (time frequency) to which the data should be resampled. This is a valid Pandas 
        offset string (e.g., 'H' for hourly, 'D' for daily). If `None`, the data will be interpolated 
        based on the time intervals of `science_data`.

    interpolation_method : str, optional
        The interpolation method to use when `data_cadence` is `None`. It can be:
        - 'time' : Interpolate based on time indices.
        - Other methods supported by Pandas' `.interpolate()` can also be used (e.g., 'linear', 'polynomial', etc.).
        Default is 'time'.

    combine_method : str, optional
        The method to use when resampling the data. It determines how to combine data when multiple entries 
        exist within the same resampling interval. Options are:
        - 'mean' : Take the mean of the values in the resampling interval.
        - 'median' : Take the median of the values in the resampling interval.
        Default is 'mean'.

    Returns
    -------
    merged_data : pandas.DataFrame
        A DataFrame with the location and science data aligned according to the specified cadence, 
        with appropriate interpolation or resampling applied.
    
    Notes
    -----
    - If `science_data` is `None`, the function will only adjust the cadence of `location_data`.
    - The function uses `np.nanmean` or `np.nanmedian` for resampling when the cadence is specified.
    - Interpolation is performed on `location_data` based on time when no specific cadence is provided.
    - The function assumes that the index of both `location_data` and `science_data` represents time.

    Example
    -------
    >>> location_data = pd.DataFrame({'x': [1, 2, 3]}, index=pd.to_datetime(['2023-01-01', '2023-01-02', '2023-01-03']))
    >>> science_data = pd.DataFrame({'y': [4, 5, 6]}, index=pd.to_datetime(['2023-01-01', '2023-01-03', '2023-01-05']))
    >>> merged_data = _match_data_cadence(location_data, science_data, data_cadence='1T', combine_method='mean')
    """
    ## set combine method function
    combine_method_dict = {
        'mean': np.nanmean, 
        'median': np.nanmedian,
    }
    combine_method_function = combine_method_dict[combine_method]
    if science_data is not None:
        merged_data = location_data.join(science_data, how='outer')
    else:
        merged_data = location_data
    if not data_cadence:
        if science_data is not None:
            for col in location_data.columns:
                merged_data[col] = merged_data[col].interpolate(interpolation_method)
            merged_data = merged_data.reindex(science_data.index)
    else:
        merged_data = merged_data.resample(data_cadence).apply(combine_method_function)
    return merged_data


def _transform_into_dataframe(list_of_data, column_names=[], combine_axis='columns'):
    """
    Convert a list of data objects into a single Pandas DataFrame, with time as the index.

    This function processes a list of data dictionaries and concatenates them into a single DataFrame. 
    Each data dictionary should contain 'x' (time) and 'y' (values) keys, where 'x' represents the 
    time values (which will be set as the DataFrame index) and 'y' contains the data values. The function 
    supports data with both single-column and multi-column 'y' values, and allows combining data 
    along either rows or columns.

    Parameters
    ----------
    list_of_data : list of dict
        A list of dictionaries (output of pySPEDAS), where each contains two keys:
        - 'x': A sequence of time values (e.g., a list, numpy array, or pandas Series).
        - 'y': A sequence of data values corresponding to the 'x' times. It can either be a 1D array 
          or a 2D array (with multiple columns of data).
    
    column_names : list of list of str, optional
        A list of lists, where each sublist contains the column names for the corresponding data in 
        `list_of_data`. If not provided, default column names will be used.
    
    combine_axis : {'columns', 'rows'}, optional
        Specifies whether the data should be concatenated along columns ('columns') or rows ('rows').
        Default is 'columns'.
    
    Returns
    -------
    combined_dataframe : pandas.DataFrame
        A DataFrame with the concatenated data, with time as the index. The values from the input data
        are arranged according to the specified `combine_axis`.
    
    combined_columns : list of str
        A list of column names for the resulting DataFrame, in the order they appear after concatenation.
    
    Notes
    -----
    - The time values ('x' in each data dictionary) are converted into pandas `DatetimeIndex`.
    - If `y` is a 1D array, it is treated as a single column of data. If `y` is a 2D array, each column
      in `y` is treated as a separate column in the resulting DataFrame.
    - The `column_names` list should have the same length as `list_of_data`. Each sublist in `column_names`
      corresponds to the column names for a particular data dictionary in the list.
    - The `combine_axis` parameter allows for flexibility in how data is combined:
        - `'columns'` will concatenate the DataFrames along the columns.
        - `'rows'` will concatenate the DataFrames along the rows.

    Example
    -------
    >>> data1 = {'x': ['2023-01-01', '2023-01-02'], 'y': [[1, 2], [3, 4]]}
    >>> data2 = {'x': ['2023-01-01', '2023-01-02'], 'y': [[5, 6], [7, 8]]}
    >>> list_of_data = [data1, data2]
    >>> column_names = [['A', 'B'], ['C', 'D']]
    >>> df, columns = _transform_into_dataframe(list_of_data, column_names, combine_axis='columns')
    >>> print(df)
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


def _validate_data(data):
    """Validates the format of input data."""
    if not isinstance(data, dict):
        raise TypeError("Input data must be a dictionary.")
    if 'x' not in data:
        raise KeyError("Input data must contain the 'x' key.")
    if 'y' not in data:
        raise KeyError("Each dictionary in the location list must contain the 'y' key.")


def _validate_coord(coord):
    """Validates the coordinate system (either from sunpy or astropy)."""
    if coord is not None:
        if coord in [v[1] for v in inspect.getmembers(astropy.coordinates, inspect.isclass)]:
            return
        if coord in [v[1] for v in inspect.getmembers(sunpy.coordinates, inspect.isclass)]:
            return
        raise TypeError("Desired coordinate must be a valid coordinate system from sunpy.coordinates or astropy.coordinates.")


def _validate_units(units):
    """Validates the units as from astropy.units."""
    if units is not None:
        if type(units) in [v[1] for v in inspect.getmembers(astropy.units, inspect.isclass)]:
            return
        return TypeError("Desired units must be a class within astropy.units.")
    

def _validate_interpolation_method(interpolation_method):
    """Validates the interpolation method."""
    ##TODO: change to not hardcoding valid methods
    valid_methods = ['linear', 'time', 'index', 'polynomial', 'nearest', 'spline', 'barycentric', 'pchip']
    if interpolation_method not in valid_methods:
        raise ValueError(f"Invalid interpolation method. Expected one of {valid_methods}, got '{interpolation_method}'.")
import inspect
import numpy as np
import pandas as pd
import astropy.coordinates
import sunpy.coordinates
import astropy.units

def _frame_has_obstime(frame):
    try:
        if 'obstime' in frame.frame_attributes:
            return True
        else:
            return False
    except Exception:
        if 'obstime' in frame.__dict__['frame_attributes']:
            return True
        else:
            return False


def _transform_into_skycoord(dataframe, frame, units, coord_kwargs={}):
    """
    Convert a Pandas DataFrame into an Astropy SkyCoord object for a specified coordinate frame,
    applying appropriate units to the coordinate data.

    Parameters
    ----------
    dataframe : pandas.DataFrame
        A DataFrame containing coordinate data, where each row represents a point (e.g., an observation),
        and each column corresponds to a coordinate (e.g., x, y, z, or latitude, longitude).
        The DataFrame's index should represent the observation times (timestamps).

    frame : sunpy.coordinates.frame or astropy.coordinates.frame
        The coordinate frame class to use for creating the SkyCoord object. 
        Typically a SunPy or Astropy frame, such as ICRS, Heliocentric, etc.

    units : astropy.units or list of astropy.units
        The units for the coordinates in the DataFrame. If a single unit is provided, it is applied to all coordinates.
        If a list of units is provided, it should have the same length as the number of columns in `dataframe`.

    coord_kwargs : dict, optional
        Additional keyword arguments to pass to the frame constructor, such as metadata or other frame-specific parameters
        (default is an empty dictionary).

    Returns
    -------
    coord : astropy.coordinates.SkyCoord
        An Astropy SkyCoord object created using the input data, units, and frame. 
        The object represents the coordinates in the specified frame with the appropriate unit conversion applied.

    Notes
    -----
    - The function assumes that the number of columns in `dataframe` corresponds to the dimensionality of the coordinate system:
      - 3 columns → xyz coordinates
      - 2 columns → latitude/longitude or similar
    
    Example
    -------
    >>> from astropy.coordinates import GCRS
    >>> from astropy.units import AU
    >>> import pandas as pd
    >>> df = pd.DataFrame([[1, 2, 3], [4, 5, 6]], columns=['x', 'y', 'z'], index=['t1', 't2'])
    >>> coord = _transform_into_skycoord(df, GCRS, AU)
    """
    data = dataframe.values
    obstimes = dataframe.index
    coord = None
    if not isinstance(units, list):
        units = [units] * data.shape[1]
    units = [v if v is not None else 1 for v in units]
    if frame is not None:
        if data.shape[1]==3:
            if _frame_has_obstime(frame)==True:
                coord = frame(data[:,0]*units[0], data[:,1]*units[1], data[:,2]*units[2], obstime=obstimes, **coord_kwargs)
            else:
                coord = frame(data[:,0]*units[0], data[:,1]*units[1], data[:,2]*units[2], **coord_kwargs)
        elif data.shape[1]==2:
            if _frame_has_obstime(frame)==True:
                coord = frame(data[:,0]*units[0], data[:,1]*units[1], obstime=obstimes, **coord_kwargs)
            else:
                coord = frame(data[:,0]*units[0], data[:,1]*units[1], **coord_kwargs)
    return coord


def _transform_into_quantity(dataframe, units):
    """
    Convert the numerical values in a Pandas DataFrame into an Astropy Quantity object by applying the specified units.

    Parameters
    ----------
    dataframe : pandas.DataFrame
        A DataFrame containing numerical data to be converted into an Astropy `Quantity`.
        Each entry in the DataFrame will be multiplied by the specified units.

    units : `astropy.units.Unit`, or `None`
        The units to apply to the data. If not `None`, each value in the DataFrame will be multiplied
        by this unit. If `None`, the function will return `None` without applying any units to the data.

    Returns
    -------
    quantity : `astropy.units.Quantity` or `None`
        An Astropy `Quantity` object with the same shape as the input DataFrame. Each entry in the 
        DataFrame is multiplied by the specified units. If `units` is `None`, the function returns `None`.

    Notes
    -----
    - This function assumes that the `units` provided are compatible with the data in the DataFrame.
    - The returned `Quantity` will have the same shape as the original DataFrame.

    Example
    -------
    >>> import pandas as pd
    >>> from astropy import units as u
    >>> df = pd.DataFrame([[1, 2], [3, 4]], columns=['a', 'b'])
    >>> result = _transform_into_quantity(df, u.meter)
    >>> print(result)
    <Quantity [[1. 2.] [3. 4.]] m>
    """
    data = dataframe.values
    quantity = None
    if units is not None:
        quantity = data * units
    return quantity

def _adjust_skycoord(skycoord, desired_coord, desired_coord_kwargs, obstime=None):
    """
    Transform an Astropy SkyCoord object to a specified coordinate frame at a given observation time.

    Parameters
    ----------
    skycoord : astropy.coordinates.SkyCoord
        An Astropy `SkyCoord` object representing the initial coordinates. This will be transformed 
        into the desired coordinate frame.

    obstime : `astropy.time.Time`
        The observation time for the transformation. This will be used to create the desired coordinate frame 
        and to ensure the transformation is applied at the correct time.

    desired_coord : astropy.coordinates.frame
        The desired coordinate frame class (e.g., GCRS, ICRS, etc.) to which the SkyCoord object should be transformed.

    desired_coord_kwargs : dict
        Additional keyword arguments to be passed to the constructor of the `desired_coord` frame. These may include 
        metadata or other frame-specific parameters.

    Returns
    -------
    transformed_skycoord : astropy.coordinates.SkyCoord
        A transformed `SkyCoord` object in the desired coordinate frame at the specified observation time.

    desired_frame : astropy.coordinates.frame
        The desired coordinate frame object created using the given `obstime` and `desired_coord_kwargs`.

    Notes
    -----
    - This function performs a coordinate transformation using the `transform_to()` method of Astropy's `SkyCoord`.
    - The function assumes that the desired coordinate frame requires an observation time (`obstime`).
    - The returned `transformed_skycoord` is a `SkyCoord` object that represents the input `skycoord` in the 
      new coordinate frame at the specified `obstime`.

    Example
    -------
    >>> from astropy.coordinates import GCRS, ICRS
    >>> from astropy import units as u
    >>> from astropy.time import Time
    >>> import astropy.coordinates as coord
    >>> skycoord = coord.SkyCoord(1, 2, 3, unit=u.kpc)
    >>> obstime = Time('2025-03-25')
    >>> transformed, frame = _adjust_skycoord(skycoord, obstime, GCRS, {'other_param': 'value'})
    >>> print(transformed)
    <SkyCoord (GCRS): (x, y, z) in kpc [1.0, 2.0, 3.0] ...>
    """
    if _frame_has_obstime(desired_coord):
        desired_frame = desired_coord(obstime=obstime, **desired_coord_kwargs)
    else:
        desired_frame = desired_coord(**desired_coord_kwargs)
    transformed_skycoord = skycoord.transform_to(desired_frame)
    return transformed_skycoord, desired_frame

def _rename_dataframe_with_skycoord_columns(df, skycoord, frame, spatial_columns, desired_units):
    """
    Rename and update the columns of a Pandas DataFrame with coordinate data from an Astropy SkyCoord object, 
    transforming the data to the specified units and frame.

    Parameters
    ----------
    df : pandas.DataFrame
        A DataFrame containing spatial data (e.g., coordinates) that needs to be updated with values from the SkyCoord object.
        The DataFrame should have columns corresponding to spatial coordinates (e.g., x, y, z).

    skycoord : astropy.coordinates.SkyCoord
        An Astropy `SkyCoord` object containing the coordinate data that will be used to update the DataFrame.
        The SkyCoord object should represent coordinates in a given frame and representation.

    frame : astropy.coordinates.frame
        The coordinate frame class (e.g., ICRS, GCRS, etc.) to which the `skycoord` object will be converted.

    spatial_columns : list of str
        A list of column names in the DataFrame that correspond to the spatial coordinates (e.g., x, y, z, etc.) 
        that will be replaced with values from the SkyCoord object.

    desired_units : str, `astropy.units.Unit`, or list of `astropy.units.Unit`
        The units to apply to the coordinate data. If a single unit is provided, it will be applied to all columns. 
        If a list of units is provided, it should have the same length as the number of spatial columns in the DataFrame.
        If `None` is provided for a unit, no transformation will be applied to that column.

    Returns
    -------
    df : pandas.DataFrame
        The updated DataFrame with the renamed columns and transformed coordinate values from the SkyCoord object.
    
    skycoord_cols : list of str
        The names of the SkyCoord representation components corresponding to the updated DataFrame columns.

    Notes
    -----
    - This function assumes that the number of `spatial_columns` matches the number of components in the SkyCoord representation.
    - If the length of `desired_units` does not match the length of the SkyCoord representation components, the function will print a warning and ignore the unit transformation.
    - The transformation to the specified `frame` is handled by converting the SkyCoord object to the desired representation before extracting values.
    - The units of the SkyCoord representation will be adjusted based on the provided `desired_units` for each component.

    Example
    -------
    >>> import pandas as pd
    >>> from astropy.coordinates import SkyCoord, ICRS
    >>> from astropy import units as u
    >>> # Example DataFrame with spatial columns
    >>> df = pd.DataFrame([[1, 2, 3], [4, 5, 6]], columns=['x', 'y', 'z'])
    >>> skycoord = SkyCoord(1, 2, 3, unit=u.kpc, frame=ICRS)
    >>> frame = ICRS
    >>> desired_units = u.parsec
    >>> df_updated, skycoord_columns = _rename_dataframe_with_skycoord_columns(df, skycoord, frame, ['x', 'y', 'z'], desired_units)
    >>> print(df_updated)
    """
    skycoord_representation = skycoord.represent_as(frame.get_representation_cls())
    skycoord_cols = skycoord.representation_info[frame.get_representation_cls()]['names']
    if not isinstance(desired_units, list):
        desired_units = [desired_units] * len(skycoord_cols)
    else:
        if len(desired_units)!=len(skycoord_cols):
            print("Length of desired units does not match length of skycoord data columns. Ignoring unit transformation.")
            desired_units = [None] * len(skycoord_cols)
    for col,scol,dunit, srattv in zip(spatial_columns, skycoord_cols, desired_units, skycoord_representation.__dict__.values()):
        try:
            scol_col = skycoord_representation.__getattribute__(scol)
        except AttributeError:
            ## get skycoord column brute force through attributes
            scol_col = srattv
        if dunit is not None:
            if dunit != scol_col.unit:
                scol_col = scol_col.to(dunit)
        vals = scol_col.value
        df[col] = vals
    df.rename(columns={c: sc for c,sc in zip(spatial_columns, skycoord_cols)}, inplace=True)
    return df, skycoord_cols


def _convert_coordinates(spatial, desired_coord=None, coord_kwargs={}, desired_units=None):
    """
    Convert the coordinates of a given location to a desired coordinate frame, if applicable, 
    and update the corresponding DataFrame and SkyCoord object.

    Parameters
    ----------
    spatial : pysplot.SpatialData
        An object that contains `data` (a Pandas DataFrame) and a `skycoord` (an Astropy `SkyCoord` object).
        The `data` represents the spatial coordinates, while the `skycoord` contains the original coordinates 
        in a particular frame.

    desired_coord : astropy.coordinates.frame, optional
        The desired coordinate frame to which the `SkyCoord` object should be transformed. 
        This should be a valid frame class such as `GCRS`, `ICRS`, `HeliographicStonyhurst`, etc. 
        If `None`, no transformation is applied, and the original coordinates are retained.

    coord_kwargs : dict, optional
        Additional keyword arguments to pass to the `desired_coord` frame constructor. This allows for frame-specific 
        parameters to be provided (default is an empty dictionary).

    desired_units : `astropy.units.Unit`, optional
        The units to apply to the transformed coordinates. If `None`, no unit transformation is applied 
        (default is `None`).

    Returns
    -------
    transformed_df : pandas.DataFrame
        A DataFrame containing the transformed coordinates in the desired coordinate frame. 
        The column names of the DataFrame will be updated to match the components of the transformed frame.

    transformed_skycoord : astropy.coordinates.SkyCoord or None
        The transformed `SkyCoord` object corresponding to the new coordinate frame. If no transformation 
        is applied (i.e., when the current coordinate frame matches the desired frame), this will be the 
        original `SkyCoord` object from `spatial`.

    Notes
    -----
    - Transformation occurs only if the `spatial.skycoord` is not `None`, a `desired_coord` is provided, 
      and the current coordinate frame (`spatial.coord`) differs from the target (`desired_coord`).
    - If no transformation is required, the function will return the original DataFrame and SkyCoord object.
    - The function assumes that the `spatial` object contains valid data with corresponding `SkyCoord` attributes.
    - If the coordinate transformation is applied, the resulting DataFrame will have updated column names 
      according to the components of the transformed coordinate frame.

    Example
    -------
    >>> from sunpy.coordinates import HeliocentricInertial, GeocentricSolarMagnetospheric
    >>> location = SpatialData(location, coord=HeliocentricInertial)
    >>> transformed_df, transformed_skycoord = _convert_coordinates(location, GeocentricSolarMagnetospheric)
    """
    transformed_df = pd.DataFrame()
    transformed_df.index = spatial.data.index
    if (spatial.skycoord is not None) and (desired_coord is not None) and (spatial.coord != desired_coord):
        transformed_skycoord, desired_frame = _adjust_skycoord(spatial.skycoord, desired_coord, coord_kwargs, obstime=spatial.data.index)
        transformed_df, _ = _rename_dataframe_with_skycoord_columns(transformed_df, transformed_skycoord, desired_frame, spatial.columns, desired_units)
    else:
        transformed_df = spatial.data
        transformed_skycoord = spatial.skycoord
    return transformed_df



def _convert_units(science, desired_units=None):
    """
    Converts the units of a `science.quantity` object to the desired units and returns the transformed data.

    This function takes a `science` object, which is expected to contain a `quantity` attribute (e.g., 
    an `astropy.Quantity` object), and transforms its values into the specified `desired_units`. If no 
    `desired_units` are provided, the current units of the `science` object are used.

    Parameters
    ----------
    science : object
        An object containing a `quantity` attribute (e.g., an `astropy.Quantity` object) with data and units.
        It should also have a `columns` attribute (representing the column names) and a `data` attribute, 
        which holds the raw data values in a `pandas.DataFrame`.

    desired_units : `astropy.units.Unit`, optional
        The units to which the `science.quantity` should be converted. If `None`, the current units of the 
        `science` object will be retained. If a unit is specified, the values in the `science.quantity` object 
        will be converted to this unit.

    Returns
    -------
    transformed_df : pandas.DataFrame
        A DataFrame with the same index as the original `science.data`, containing the transformed values in 
        the desired units. If no units transformation is needed (i.e., the current and desired units are the same),
        the original data is returned unchanged.

    Notes
    -----
    - If the `science` object does not have units (i.e., `science.units is None`), the function returns 
      the original `science.data` without any unit conversion.
    - If the current units of the `science.quantity` object are the same as the `desired_units`, no transformation 
      occurs, and the original values are returned unchanged.
    - This function assumes that the `science.quantity` object has the appropriate unit conversion functionality 
      (e.g., through Astropy's `Quantity.to()` method).

    Example
    -------
    >>> from astropy import units as u
    >>> import pandas as pd
    >>> # Assume `science` is an object with a `quantity` and `data` attributes
    >>> science = SomeScienceObject(quantity=some_quantity, data=pd.DataFrame([[1, 2], [3, 4]]), columns=['a', 'b'])
    >>> transformed_df = _convert_units(science, desired_units=u.meter)
    >>> print(transformed_df)
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


def match_data_cadence(interp_data, standard_data, data_cadence=None, interpolate_columns='all', interpolation_method='time', combine_method='mean'):
    """
    Match the cadence (time intervals) of `interp_data` (location data) and `standard_data` (science data) 
    by resampling or interpolating the data as necessary to align their time intervals. 
    If the cadences are mismatched, this function either resamples or interpolates the `interp_data` to match 
    the cadence of the `standard_data`, or applies a specified combine method to downsample the data.

    Parameters
    ----------
    interp_data : pandas.DataFrame
        A DataFrame containing the location data, with its index representing the time or cadence of the data.
    
    standard_data : pandas.DataFrame, optional
        A DataFrame containing the science data to be aligned with `interp_data`. Its index represents 
        the time or cadence for the science data. If `None`, the function will only adjust `interp_data` to the 
        specified cadence.

    data_cadence : str, optional
        The desired cadence (time frequency) to which the data should be resampled. This is a valid Pandas 
        offset string (e.g., 'H' for hourly, 'D' for daily). If `None`, the data will be interpolated 
        based on the time intervals of `standard_data`.

    interpolate_columns : str or list of str, optional
        The columns of `interp_data` to be interpolated. By default, it will interpolate all columns in `interp_data`. 
        If a specific list of columns is provided, only those columns will be interpolated.

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
    - If `standard_data` is `None`, the function will only adjust the cadence of `interp_data`.
    - The function uses `np.nanmean` or `np.nanmedian` for resampling when `data_cadence` is specified.
    - If `data_cadence` is `None`, interpolation is performed on `interp_data` based on the time intervals of `standard_data`.
    - The function assumes that the index of both `interp_data` and `standard_data` represents time.
    - The interpolation and resampling methods allow for flexible handling of time series data with different cadences.
    - If `interpolate_columns` is specified as 'all', all columns in `interp_data` will be interpolated. 
      Otherwise, only the specified columns will be interpolated.

    Example
    -------
    >>> location_data = pd.DataFrame({'x': [1, 2, 3]}, index=pd.to_datetime(['2023-01-01', '2023-01-02', '2023-01-03']))
    >>> science_data = pd.DataFrame({'y': [4, 5, 6]}, index=pd.to_datetime(['2023-01-01', '2023-01-03', '2023-01-05']))
    >>> merged_data = match_data_cadence(location_data, science_data, data_cadence='1T', combine_method='mean')
    """
    ## set combine method function
    combine_method_dict = {
        'mean': np.nanmean, 
        'median': np.nanmedian,
    }
    combine_method_function = combine_method_dict[combine_method]
    if standard_data is not None:
        merged_data = interp_data.join(standard_data, how='outer')
    else:
        merged_data = interp_data
    if not data_cadence:
        if standard_data is not None:
            if interpolate_columns=='all':
                interpolate_columns = interp_data.columns
            #else:
                # check that interpolate_columns is list   
            for col in interpolate_columns:
                merged_data[col] = merged_data[col].interpolate(interpolation_method)
            merged_data = merged_data.reindex(standard_data.index)
    else:
        merged_data = merged_data.resample(data_cadence).apply(combine_method_function)
    return merged_data


def _transform_into_dataframe(list_of_data, column_names=[], combine_axis='columns'):
    """
    Convert a list of data dictionaries into a single Pandas DataFrame, with time as the index.

    This function processes a list of data dictionaries, each containing time ('x') and corresponding 
    data values ('y'). The data is then concatenated into a single DataFrame, with the time values as 
    the index. The function supports both single-column and multi-column data, and allows the data to 
    be combined along rows or columns.

    Parameters
    ----------
    list_of_data : list of dict
        A list of dictionaries, where each dictionary represents a dataset with two keys:
        - 'x': A sequence of time values (e.g., a list, numpy array, or pandas Series).
        - 'y': A sequence of data values corresponding to 'x'. 'y' can be either a 1D array 
          (single-column data) or a 2D array (multi-column data).
    
    column_names : list of list of str, optional
        A list of lists, where each sublist contains the column names for the corresponding data in 
        `list_of_data`. Each sublist should have the same number of elements as the number of columns in 
        the corresponding 'y' data. If not provided, default column names will be used.
    
    combine_axis : {'columns', 'rows'}, optional
        Specifies whether the data should be concatenated along columns ('columns') or rows ('rows').
        Default is 'columns'. 
        - 'columns' means the data from each dictionary will be concatenated as new columns.
        - 'rows' means the data from each dictionary will be stacked as new rows.

    Returns
    -------
    combined_dataframe : pandas.DataFrame
        A DataFrame with the concatenated data, with time as the index. The data from the input dictionaries 
        is arranged according to the specified `combine_axis`. The DataFrame's index is set to the time values 
        ('x') from the dictionaries.

    combined_columns : list of str
        A list of column names for the resulting DataFrame, in the order they appear after concatenation.

    Notes
    -----
    - The time values ('x' in each data dictionary) are converted into a pandas `DatetimeIndex`.
    - If 'y' is a 1D array, it is treated as a single column of data. If 'y' is a 2D array, each column in 'y' 
      is treated as a separate column in the resulting DataFrame.
    - The `column_names` list should match the length of `list_of_data`. Each sublist in `column_names` 
      corresponds to the column names for a particular data dictionary in the list.
    - The `combine_axis` parameter determines how the data is combined:
        - `'columns'`: Concatenate data along columns (add new columns for each dataset).
        - `'rows'`: Concatenate data along rows (stack datasets as additional rows).
        
    Example
    -------
    >>> data1 = {'x': ['2023-01-01', '2023-01-02'], 'y': [[1, 2], [3, 4]]}
    >>> data2 = {'x': ['2023-01-01', '2023-01-02'], 'y': [[5, 6], [7, 8]]}
    >>> list_of_data = [data1, data2]
    >>> column_names = [['A', 'B'], ['C', 'D']]
    >>> df, columns = _transform_into_dataframe(list_of_data, column_names, combine_axis='columns')
    >>> print(df)
                    A  B  C  D
    2023-01-01  1  2  5  6
    2023-01-02  3  4  7  8
    """
    combined_dataframe = pd.DataFrame()
    combined_columns = []
    for data, columns in zip(list_of_data, column_names):
        dataframe = pd.DataFrame()
        dataframe.index = pd.to_datetime(data['x'])
        if len(np.shape(data['y']))==1:
            ncol=1
        else:
            ncol = np.shape(data['y'])[1]
        if ncol==1:
            dataframe[columns[0]] = data['y']
        else:
            data['y'] = np.array(data['y'])
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

def _check_column_length(data, column_names, prefix=''):
    """Checks that columns provided are the correct length given the data. Otherwise create default column names."""
    column_names_checked = []
    if len(data) > len(column_names):
        n_column_names_missing = len(data) - len(column_names)
        for _ in range(n_column_names_missing):
            column_names.append([])
    for d, dat in enumerate(data):
        dcols = column_names[d]
        if len(np.shape(dat['y']))==1:
            ncols=1
        else:
            ncols=np.shape(dat['y'])[1]
        dcols_length = len(dcols)
        if dcols_length!=ncols:
            dcols = [f'{prefix}_data{d}_{i}' for i in range(1,ncols+1)]
        column_names_checked.append(dcols)
    return column_names_checked

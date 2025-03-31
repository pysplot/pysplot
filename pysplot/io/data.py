from pysplot.io.data_helpers import _transform_into_dataframe, _transform_into_skycoord, _transform_into_quantity, _convert_coordinates, _convert_units, match_data_cadence, _validate_coord, _validate_data, _validate_units, _validate_interpolation_method, _adjust_skycoord, _rename_dataframe_with_skycoord_columns
import pandas as pd

class SpatialData:
    """
    A class to represent spatial data with associated time, location, and coordinate system information.

    This class is designed to store and manipulate location data, convert it into a Pandas DataFrame, and 
    represent the data in a specific coordinate system using Astropy's `SkyCoord`. It supports resampling, 
    interpolation, and transforming the data to different coordinate frames, allowing for efficient spatial data 
    analysis and manipulation.

    Attributes
    ----------
    data : pandas.DataFrame
        A DataFrame containing the location data, with time as the index. Each row represents a data point 
        with its associated time and spatial coordinates.
    
    columns : list of str
        A list of column names corresponding to the data in the `data` DataFrame. These columns represent 
        the different spatial coordinate components (e.g., 'x', 'y', 'z').

    coord : astropy.coordinates.frame, optional
        The coordinate frame used for transforming the data (e.g., GCRS, ICRS, Heliocentric). If `None`, 
        the data remains in its original frame. This frame will be used for transforming the location 
        data into a different coordinate system.
    
    units : astropy.units or list of astropy.units, optional
        The units of the coordinates in the `data` DataFrame, applied during transformation. If provided, 
        the data will be converted into the specified units when transforming into the desired coordinate system.
    
    interpolation_method : str, optional
        The method used for interpolating the data if necessary (e.g., 'linear'). This interpolation method 
        will be applied when resampling the data.

    skycoord : astropy.coordinates.SkyCoord
        An Astropy `SkyCoord` object representing the location data in the specified coordinate frame. 
        This object provides methods to manipulate the spatial coordinates and convert them between frames.
    
    Methods
    -------
    __init__(spatial, coord=None, units=None, combine_axis='columns', interpolation_method='linear', spatial_columns_names=[], coord_kwargs={}, drop_duplicates=True)
        Initializes the `SpatialData` object by transforming the provided location data into a DataFrame 
        and converting it into the specified coordinate system. This method validates input data, transforms 
        coordinates, and initializes the attributes of the object.

    Notes
    -----
    - The `spatial` parameter is expected to be a list of dictionaries containing 'x' (time) and 'y' (coordinate values). 
      If a single object is provided, it will be automatically converted into a list.
    - The `coord` parameter should specify a valid Astropy frame, such as `ICRS`, `GCRS`, `Heliocentric`, etc.
    - The `spatial_columns_names` parameter can be used to specify custom column names for the DataFrame. If not provided, 
      default column names will be assigned.
    - If no coordinate system (`coord`) is provided, the data will remain in its original frame and will not be transformed.
    - The `combine_axis` parameter controls how the data is combined when creating the DataFrame. It can concatenate data 
      along rows ('rows') or columns ('columns').
    - The `drop_duplicates` flag determines whether duplicate data points should be removed from the DataFrame.
    - The `units` parameter allows the transformation of coordinates into the specified physical units (e.g., AU, km, etc.).

    Example
    -------
    >>> spatial_data = [{'x': [0, 1, 2], 'y': [[1, 2], [3, 4], [5, 6]]}]
    >>> spatial_data_obj = SpatialData(spatial_data, coord=sunpy.coordinates.Heliocentric, units=u.AU, combine_axis='columns')
    >>> print(spatial_data_obj.data)
    >>> print(spatial_data_obj.skycoord)

    """
    def __init__(self, spatial, coord=None, units=None, combine_axis='columns', interpolation_method='linear', spatial_columns_names=[], coord_kwargs={}, drop_duplicates=True):        
        """
        Initializes the `SpatialData` object by transforming the provided location data into a DataFrame 
        and converting it into the specified coordinate system.

        Parameters
        ----------
        spatial : list of dicts
            A list of dictionaries, where each dictionary contains 'x' (time) and 'y' (coordinate values). 
            If a single dictionary is provided, it is automatically converted into a list.
        
        coord : astropy.coordinates.frame, optional
            The coordinate system to transform the spatial data into. Valid options include Astropy coordinate frames 
            such as `ICRS`, `GCRS`, `Heliocentric`, etc. If `None`, no transformation is performed.

        units : astropy.units or list of astropy.units, optional
            The units to apply to the coordinates during transformation. If `None`, the data will remain in the original units.

        combine_axis : str, optional
            Specifies whether to combine the data along rows ('rows') or columns ('columns'). Default is `'columns'`.
        
        interpolation_method : str, optional
            The interpolation method to use when resampling the data. Default is `'linear'`.
        
        spatial_columns_names : list of str, optional
            A list of column names for the DataFrame. If not provided, default names will be used.
        
        coord_kwargs : dict, optional
            Additional keyword arguments to pass when transforming coordinates (e.g., frame-specific parameters).
        
        drop_duplicates : bool, optional
            If `True`, duplicate data points are removed from the DataFrame. Default is `True`.
        """
        if not isinstance(spatial, list):
            spatial = [spatial]
            spatial_columns_names = [spatial_columns_names]
        
        # do data validation checks
        for s in spatial:
            _validate_data(s)
        if coord is not None:
            _validate_coord(coord)
        if units is not None:
            _validate_units(units)
        _validate_interpolation_method(interpolation_method)
        
        dataframe, columns = _transform_into_dataframe(spatial, column_names=spatial_columns_names, combine_axis=combine_axis)
        if drop_duplicates:
            dataframe = dataframe.drop_duplicates(keep='first')
        skycoord = _transform_into_skycoord(dataframe, coord, units, coord_kwargs=coord_kwargs)
        self.data = dataframe
        self.columns = columns
        self.coord = coord
        self.units = units
        self.interpolation_method = interpolation_method
        self.skycoord = skycoord


class ScienceData:
    """
    A class to represent scientific data with associated values and units.

    This class is designed to hold scientific data, convert it into a Pandas DataFrame, and provide 
    a representation of the data as an Astropy `Quantity` object, with optional unit conversion. 
    It allows for unit-aware operations on the data and facilitates manipulation and analysis.

    Attributes
    ----------
    data : pandas.DataFrame
        A DataFrame containing the scientific data, with time or index as the index and the corresponding 
        values in the columns. This DataFrame forms the core structure of the `ScienceData` object.

    columns : list of str
        A list of column names corresponding to the data in the DataFrame. These names describe the nature of 
        the data and are either provided explicitly or generated automatically.

    units : str or list of str, optional
        The units to apply to the scientific data. This can be a single unit applied to all columns, 
        or a list of units matching the number of columns. If no units are provided, the data is treated as unitless.

    quantity : astropy.units.Quantity
        An Astropy `Quantity` object representing the scientific data in the specified units. This allows for 
        unit-aware operations on the data, such as unit conversions or arithmetic operations while respecting the units.

    Methods
    -------
    __init__(science, units=None, science_columns_names=[], drop_duplicates=True)
        Initializes the `ScienceData` object by transforming the provided scientific data into a DataFrame 
        and converting it into the specified units as an Astropy `Quantity`.

    Notes
    -----
    - The class assumes that the input `science` data is a list of dictionaries containing 'x' (time) 
      and 'y' (values) keys. If a single dataset is provided, it will be automatically converted into a list.
    - The `units` parameter specifies the units for the data, which will be applied to each column in the 
      DataFrame. If no units are provided, the data is treated as unitless.
    - The `science_columns_names` parameter can be used to specify the column names for the DataFrame. 
      If not provided, default column names will be generated.
    - The `quantity` attribute represents the data in unit-aware format using Astropy's `Quantity` class. 
      It is designed to facilitate unit conversion and unit-aware operations on the data.

    Example
    -------
    >>> science_data = [{'x': [0, 1, 2], 'y': [[10, 20], [30, 40], [50, 60]]}]
    >>> science = ScienceData(science_data, units='km/s', science_columns_names=[['Velocity X', 'Velocity Y']])
    >>> print(science.data)
       Velocity X  Velocity Y
    0          10          20
    1          30          40
    2          50          60
    >>> print(science.quantity)
    <Quantity [10., 30., 50.] km / s>
    
    Notes on Example:
    - The `science_data` input is a list of dictionaries where 'x' represents time or index values, 
      and 'y' represents the associated data values.
    - The `units` are specified as `'km/s'`, which are applied to all columns in the DataFrame.
    - The `science_columns_names` parameter is used to assign custom column names (`'Velocity X'`, `'Velocity Y'`).
    - The resulting `data` attribute contains the data as a DataFrame, while `quantity` contains the data in an 
      Astropy `Quantity` object with the specified units.
    """
    
    def __init__(self, science, units=None, science_columns_names=[], drop_duplicates=True):
        """
        Initializes the `ScienceData` object by converting the provided scientific data into a DataFrame 
        and converting it into the specified units as an Astropy `Quantity`.

        Parameters
        ----------
        science : list of dicts
            A list of dictionaries, where each dictionary contains 'x' (time) and 'y' (values) keys. 
            If a single dictionary is provided, it is automatically converted into a list.
        
        units : str or list of str, optional
            The units to apply to the scientific data. If not provided, the data is treated as unitless.
        
        science_columns_names : list of str, optional
            A list of column names for the DataFrame. If not provided, default names will be generated.
        
        drop_duplicates : bool, optional
            If `True`, duplicate rows in the data will be dropped. Default is `True`.
        """
        ## TODO: add coordinates to science values
        if not isinstance(science, list):
            science = [science]
            science_columns_names = [science_columns_names]

        # do data validation checks
        for s in science:
            _validate_data(s)
        if units is not None:
            _validate_units(units)

        dataframe, columns = _transform_into_dataframe(science, column_names=science_columns_names)
        if drop_duplicates:
            dataframe = dataframe.drop_duplicates(keep='first')
        self.data = dataframe
        self.columns = columns
        self.units = units
        self.quantity = _transform_into_quantity(dataframe, units)



class SpatialTimeData:
    """
    A class to represent spatiotemporal data by combining location and science data for analysis and plotting.

    This class integrates spatial data (e.g., position) with science data (e.g., measurements) over time, 
    allowing for advanced analysis and visualization of spatiotemporal data. It supports transforming both location 
    and science data, resampling data to match cadences, and performing unit conversions. The combined dataset is 
    ready for spatiotemporal plotting and further analysis.

    Attributes
    ----------
    spatial : SpatialData
        The spatial data, represented as a `SpatialData` object, which typically contains spatial coordinates or positions over time.
    
    science : ScienceData, optional
        The science data, represented as a `ScienceData` object, which typically contains measurements or scientific quantities 
        associated with the location data. Defaults to `None` if no science data is provided.
    
    interpolation_method : str, optional
        The method to use for interpolating the location data to match the science data's time cadence. 
        Default is `'linear'`.
    
    data_cadence : str, optional
        The time cadence to which the data should be resampled. If `None`, the data cadence is not changed.
    
    combine_method : str, optional
        The method used to combine location and science data when merging. Can be either `'mean'` or `'median'`.
        Default is `'mean'`.
    
    desired_science_units : astropy.units, optional
        The units to which the science data should be converted. If not provided, the original units are retained.
    
    desired_coord : astropy.coordinates or sunpy.coordinates, optional
        The coordinate system to which the location data should be transformed. If not provided, the original coordinate system is retained.
    
    desired_coord_kwargs : dict, optional
        Additional keyword arguments to pass when transforming the location data to the desired coordinate system.

    data : pandas.DataFrame
        A DataFrame containing the merged and processed spatiotemporal data, combining both location and science data.
    
    coord : astropy.coordinates.frame
        The desired coordinate system of the location data after transformation.
    
    skycoord : astropy.coordinates.SkyCoord
        An Astropy `SkyCoord` object representing the transformed location data in the desired coordinate system.
    
    spatial_columns : list of str
        The column names of the spatial data after transformation.
    
    science_columns : list of str, optional
        The column names of the science data after transformation. If no science data is provided, this will be `None`.

    Methods
    -------
    __init__(spatial, science=None, interpolation_method='linear', data_cadence=None, 
             combine_method='mean', desired_science_units=None, desired_coord=None, desired_coord_kwargs={})
        Initializes the `SpatialTimeData` object by processing location and science data, resampling them to 
        match cadences, and performing unit and coordinate transformations as needed.
    
    _combine_data(interpolation_method, data_cadence, combine_method, desired_coord, 
                  desired_science_units, desired_coord_kwargs)
        Prepares the location and science data for spatiotemporal plotting by combining and transforming them.
    
    transform_coord(desired_coord, desired_coord_kwargs={}, desired_units=None)
        Transforms the coordinate system of the spatial data to the specified `desired_coord` and updates the object with the new transformation.

    Notes
    -----
    - The `spatial` and `science` data should be iterable (e.g., a list or pandas DataFrame).
    - The `interpolation_method` determines how to interpolate the location data to match the science data time.
    - The `combine_method` dictates how merged datasets are combined: `'mean'` computes the mean, while 
      `'median'` computes the median.
    - If no `science` data is provided, only the spatial data will be processed.
    - The `desired_science_units` and `desired_coord` allow for transformations and unit conversions of the science 
      and location data, respectively.
    - The resulting merged dataset is stored in the `data` attribute, ready for analysis or plotting.
    - The `transform_coord` method can be used to update the coordinate system of the location data.

    Example
    -------
    >>> spatial_data_object = SpatialData(...)
    >>> science_data_object = ScienceData(...)
    >>> splot_data_obj = SpatialTimeData(spatial_data_object, science_data_object=science_data, 
    >>>                                      interpolation_method='linear', 
    >>>                                      data_cadence='1T', combine_method='mean')
    >>> print(splot_data_obj.data)  # Merged and processed spatiotemporal data
    >>> print(splot_data_obj.skycoord)  # Transformed location data in the desired coordinate system

    """
    def __init__(self, spatial, science=None, interpolation_method='linear', data_cadence=None, combine_method='mean', desired_science_units=None, desired_coord=None, desired_coord_kwargs={}, desired_spatial_units=None):
        
        # do data validation checks
        if not isinstance(spatial, SpatialData):
            raise TypeError("Input should be a SpatialData object.")
        if science is not None:
            if not isinstance(science, ScienceData):
                raise TypeError("Input should be a ScienceData object.")
        if combine_method.lower() not in ['mean', 'median']:
            raise TypeError("Combine method {combine_method} is not an acceptable method. Must be 'mean' or 'median'.")
        _validate_units(desired_spatial_units)
        _validate_units(desired_science_units)
        _validate_coord(desired_coord)

        self.spatial = spatial
        self.science = science
        self.interpolation_method = interpolation_method
        self._combine_data(interpolation_method, data_cadence, combine_method, desired_coord, desired_spatial_units, desired_coord_kwargs, desired_science_units)
        

    def _combine_data(self, interpolation_method, data_cadence, combine_method, desired_coord, desired_spatial_units, desired_coord_kwargs, desired_science_units):
        """
        Prepare data for spatiotemporal plotting by transforming and combining spatial and science data.

        This method applies transformations like coordinate conversions, unit conversions, and matching the 
        time cadence between location and science datasets. The final merged and processed data is ready 
        for analysis and visualization.

        Parameters
        ----------
        interpolation_method : str
            The interpolation method used to resample the location data to match the science data cadence.
        
        data_cadence : str or None
            The target time cadence for resampling the merged data. If `None`, no resampling is applied.
        
        combine_method : str
            The method used to combine the location and science data when merging them. 
            Can be 'mean' or 'median'.
        
        desired_coord : astropy.coordinates or sunpy.coordinates, optional
            The coordinate system to which the location data will be transformed.
        
        desired_science_units : astropy.units, optional
            The units for the science data. If `None`, the original units are kept.
        
        desired_coord_kwargs : dict, optional
            Additional keyword arguments passed to the transformation function for the coordinate system.
        
        Returns
        -------
        self : SpatialTimeData
            The `SpatialTimeData` object, with the processed and merged data stored in the `data` attribute.
        """
        if desired_coord is None:
            desired_coord = self.spatial.coord
        transformed_location_df = _convert_coordinates(self.spatial, desired_coord=desired_coord, coord_kwargs=desired_coord_kwargs, desired_units=desired_spatial_units)
        transformed_location_columns = transformed_location_df.columns
        if desired_spatial_units is None:
            desired_spatial_units = self.spatial.units
        if self.science is not None:
            if desired_science_units is None:
                desired_science_units = self.science.units
            transformed_science_df = _convert_units(self.science, desired_units=desired_science_units)
            transformed_science_columns = transformed_science_df.columns
        else:
            transformed_science_df = None
            transformed_science_columns = None
        merged_dataframe = match_data_cadence(transformed_location_df, transformed_science_df, interpolation_method=interpolation_method, data_cadence=data_cadence, combine_method=combine_method)
        merged_skycoord = _transform_into_skycoord(merged_dataframe[transformed_location_columns], desired_coord, desired_spatial_units, coord_kwargs=desired_coord_kwargs)
        self.data = merged_dataframe
        self.spatial_units = desired_spatial_units
        self.science_units = desired_science_units
        self.coord = desired_coord
        self.skycoord = merged_skycoord
        self.spatial_columns = transformed_location_columns
        self.science_columns = transformed_science_columns
        return self
    
    def transform_coord(self, desired_coord, desired_coord_kwargs={}, desired_units=None):
        """
        Transform the coordinate system of the spatial data to the specified `desired_coord`.

        This method allows users to update the coordinate system of the location data and transforms the 
        associated spatial data accordingly.

        Parameters
        ----------
        desired_coord : astropy.coordinates.frame
            The desired coordinate frame for the spatial data.
        
        desired_coord_kwargs : dict, optional
            Additional keyword arguments passed to the transformation function for the coordinate system.
        
        desired_units : astropy.units, optional
            The units for the spatial data. If `None`, the original units are retained.
        
        Returns
        -------
        self : SpatialTimeData
            The `SpatialTimeData` object with the updated coordinate system and transformed data.
        """
        current_skycoord = self.skycoord
        obstime = self.data.index
        transformed_skycoord, desired_frame = _adjust_skycoord(current_skycoord, desired_coord, desired_coord_kwargs, obstime=obstime)
        transformed_df, spatial_columns = _rename_dataframe_with_skycoord_columns(self.data, transformed_skycoord, desired_frame, self.spatial_columns, desired_units)
        self.data = transformed_df
        self.coord = desired_coord
        self.skycoord = transformed_skycoord
        self.spatial_columns = spatial_columns
        return self
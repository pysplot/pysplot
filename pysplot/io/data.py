from pysplot.io.data_helpers import _transform_into_dataframe, _transform_into_skycoord, _transform_into_quantity, _convert_coordinates, _convert_units, _match_data_cadence, _validate_coord, _validate_data, _validate_units, _validate_interpolation_method


class SpatialData:
    """
    A class to represent spatial data with associated time, location, and coordinate system information.

    This class is designed to hold location data, convert it into a Pandas DataFrame, and provide a 
    representation of the data in a specific coordinate system using Astropy's `SkyCoord`. It allows 
    for resampling, interpolation, and transforming the data to different coordinate frames.

    Attributes
    ----------
    data : pandas.DataFrame
        A DataFrame containing the location data, with time as the index.
    
    columns : list of str
        A list of column names corresponding to the data in the DataFrame.
    
    coord : sunpy.coordinate or astropy.coordinate, optional
        The coordinate frame used for transforming the data (e.g., GCRS, ICRS).
    
    units : astropy.units or list of astropy.units, optional
        The units of the coordinates in the DataFrame, applied to the data during transformation.
    
    interpolation_method : str, optional
        The method used for interpolating the data if necessary (e.g., linear).
    
    skycoord : astropy.coordinates.SkyCoord
        An Astropy `SkyCoord` object representing the location data in the specified coordinate frame.
    
    Methods
    -------
    __init__(spatial, coord=None, units=None, combine_axis='columns', interpolation_method='linear', spatial_columns_names=[], coord_kwargs={})
        Initializes the SpatialData object by transforming the provided location data into a DataFrame 
        and converting it into the specified coordinate system.
    
    Notes
    -----
    - The class assumes that the input `location` data is a list of dictionaries containing 'x' (time) 
      and 'y' (coordinate values).
    - If `location` is provided as a single object, it is automatically converted into a list.
    - The `coord` parameter should specify a valid Astropy frame, such as `ICRS`, `GCRS`, etc.
    - The `spatial_columns_names` parameter can be used to specify column names for the DataFrame.
    - If no coordinate system (`coord`) is provided, the data will be left in its original frame.
    - The `combine_axis` parameter allows data to be combined either by rows or columns when creating 
      the DataFrame.

    Example
    -------
    >>> spatial_data = [{'x': [0, 1, 2], 'y': [[1, 2], [3, 4], [5, 6]]}]
    >>> spatial_data_obj = SpatialData(spatial_data, coord=sunpy.coordinates.Heliocentric, units=u.AU, combine_axis='columns')
    >>> print(spatial_data_obj.data)
    >>> print(spatial_data_obj.skycoord)
    """
    def __init__(self, spatial, coord=None, units=None, combine_axis='columns', interpolation_method='linear', spatial_columns_names=[], coord_kwargs={}):        
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

    Attributes
    ----------
    data : pandas.DataFrame
        A DataFrame containing the scientific data, with time as the index and the values in the columns.
    
    columns : list of str
        A list of column names corresponding to the data in the DataFrame.
    
    units : str or list of str, optional
        The units to apply to the scientific data. This can be a single unit applied to all columns, 
        or a list of units matching the number of columns.
    
    quantity : astropy.units.Quantity
        An Astropy `Quantity` object representing the scientific data in the specified units. This allows 
        for unit-aware operations on the data.

    Methods
    -------
    __init__(science, units=None, science_columns_names=[])
        Initializes the ScienceData object by transforming the provided scientific data into a DataFrame 
        and converting it into the specified units as an Astropy `Quantity`.
    
    Notes
    -----
    - The class assumes that the input `science` data is a list of dictionaries containing 'x' (time) 
      and 'y' (values) keys. If a single dataset is provided, it will be automatically converted into a list.
    - The `units` parameter specifies the units for the data, which will be applied to each column in the 
      DataFrame. If no units are provided, the data is treated as unitless.
    - The `science_columns_names` parameter can be used to specify the column names for the DataFrame.
    - The `quantity` attribute represents the data in unit-aware format using Astropy's `Quantity` class.

    Example
    -------
    >>> science_data = [{'x': [0, 1, 2], 'y': [[10, 20], [30, 40], [50, 60]]}]
    >>> science = ScienceData(science_data, units='km/s', science_columns_names=[['Velocity X', 'Velocity Y']])
    >>> print(science.data)
    >>> print(science.quantity)
    """
    def __init__(self, science, units=None, science_columns_names=[]):
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
        self.data = dataframe
        self.columns = columns
        self.units = units
        self.quantity = _transform_into_quantity(dataframe, units)


class SpatialTimeData:
    """
    A class to represent spatiotemporal data by combining location and science data for analysis and plotting.

    This class integrates spatial data (e.g., position) with science data (e.g., measurements) over time.
    It allows for transformation of both location and science data, resampling the data to match cadences,
    and performing unit conversion. The combined dataset is ready for spatiotemporal plotting and further analysis.

    Attributes
    ----------
    spatial : object
        Spatial data, in the form of a SpatialData object, which typically contains spatial coordinates or positions over time.
    
    science : object, optional
        Science data, in the form of the ScienceData object, which typically contains measurements or scientific quantities associated with the 
        location data. Defaults to `None` if no science data is provided.
    
    interpolation_method : str, optional
        The method to use for interpolating location data to match the science data's time cadence. 
        Default is `'linear'`.
    
    data_cadence : str, optional
        The time cadence to which the data should be resampled. If `None`, the data cadence is not changed.
    
    combine_method : str, optional
        The method used to combine location and science data when merging. Can be either `'mean'` or `'median'`.
        Default is `'mean'`.
    
    desired_units : astropy.units, optional
        The units to which the science data should be converted.
    
    desired_coord : astropy.coordinates or sunpy.coordinates, optional
        The coordinate system to which the location data should be transformed.
    
    desired_coord_kwargs : dict, optional
        Additional keyword arguments to pass when transforming the location data to the desired coordinate system.

    Methods
    -------
    __init__(spatial, science=None, interpolation_method='linear', data_cadence=None, 
             combine_method='mean', desired_units=None, desired_coord=None, desired_coord_kwargs={})
        Initializes the SpatialTimeData object by processing location and science data, resampling them to 
        match cadences, and performing unit and coordinate transformations as needed.

    _combine_data(interpolation_method, data_cadence, combine_method, desired_coord, 
                  desired_units, desired_coord_kwargs)
        Prepares the location and science data for spatiotemporal plotting by combining and transforming them.
    
    Notes
    -----
    - The `spatial` and `science` data should be iterable (e.g., a list or pandas DataFrame).
    - The `interpolation_method` determines how to interpolate the location data to match the science data time.
    - The `combine_method` dictates how merged datasets are combined: `'mean'` computes the mean, while 
      `'median'` computes the median.
    - If no `science` data is provided, only the spatial data will be processed.
    - The `desired_units` and `desired_coord` allow for transformations and unit conversions of the science 
      and location data, respectively.
    - The resulting merged dataset is stored in the `data` attribute, ready for analysis or plotting.

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
    def __init__(self, spatial, science=None, interpolation_method='linear', data_cadence=None, combine_method='mean', desired_units=None, desired_coord=None, desired_coord_kwargs={}):
        
        # do data validation checks
        if not isinstance(spatial, SpatialData):
            raise TypeError("Input should be a SpatialData object.")
        if science is not None:
            if not isinstance(science, ScienceData):
                raise TypeError("Input should be a ScienceData object.")
        if combine_method.lower() not in ['mean', 'median']:
            raise TypeError("Combine method {combine_method} is not an acceptable method. Must be 'mean' or 'median'.")
        _validate_units(desired_units)
        _validate_coord(desired_coord)

        self.spatial = spatial
        self.science = science
        self.interpolation_method = interpolation_method
        self._combine_data(interpolation_method, data_cadence, combine_method, desired_coord, desired_units, desired_coord_kwargs)
        

    def _combine_data(self, interpolation_method, data_cadence, combine_method, desired_coord, desired_units, desired_coord_kwargs):
        """
        Prepare data read from pyspedas for spatiotemporal plotting.

        This method processes the location and science data by applying transformations such as coordinate 
        conversions, unit conversions, and matching the time cadence between location and science datasets.
        The final processed data is ready for analysis and plotting.

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
        
        desired_units : astropy.units, optional
            The units for the science data. If `None`, the original units are kept.
        
        desired_coord_kwargs : dict, optional
            Additional keyword arguments passed to the transformation function for the coordinate system.
        
        Returns
        -------
        self : SpatialTimeData
            The SpatialTimeData object, with the processed and merged data stored in the `data` attribute.
        """
        if desired_coord is None:
            desired_coord = self.spatial.coord
        transformed_location_df, transformed_location_skycoord = _convert_coordinates(self.spatial, desired_coord=desired_coord, coord_kwargs=desired_coord_kwargs)
        transformed_location_columns = transformed_location_df.columns
        if self.science is not None:
            if desired_units is None:
                desired_units = self.science.units
            transformed_science_df = _convert_units(self.science, desired_units=desired_units)
            transformed_science_columns = transformed_science_df.columns
        else:
            transformed_science_df = None
            transformed_science_columns = None
        merged_dataframe = _match_data_cadence(transformed_location_df, transformed_science_df, interpolation_method=interpolation_method, data_cadence=data_cadence, combine_method=combine_method)
        self.data = merged_dataframe
        self.units = desired_units
        self.coord = desired_coord
        self.skycoord = transformed_location_skycoord
        self.spatial_columns = transformed_location_columns
        self.science_columns = transformed_science_columns
        return self
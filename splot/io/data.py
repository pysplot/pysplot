#from pysplot.io.data_helpers import _transform_into_dataframe, _transform_to_skycoord, _transform_into_quantity, _convert_coordinates, _convert_units, _match_data_cadence
from data_helpers import _transform_into_dataframe, _transform_to_skycoord, _transform_into_quantity, _convert_coordinates, _convert_units, _match_data_cadence


class SpatialData:
    def __init__(self, location, coord=None, units=None, combine_axis='columns', interpolation_method='linear', location_columns_names=[], coord_kwargs={}):        
        if not isinstance(location, list):
            location = [location]
            location_columns_names = [location_columns_names]
        dataframe, columns = _transform_into_dataframe(location, column_names=location_columns_names, combine_axis=combine_axis)
        skycoord = _transform_to_skycoord(dataframe, coord, units, coord_kwargs=coord_kwargs)
        self.data = dataframe
        self.columns = columns
        self.coord = coord
        self.units = units
        self.interpolation_method = interpolation_method
        self.skycoord = skycoord


class ScienceData:
    def __init__(self, science, units=None, science_columns_names=[]):
        if not isinstance(science, list):
            science = [science]
            science_columns_names = [science_columns_names]
        dataframe, columns = _transform_into_dataframe(science, column_names=science_columns_names)
        self.data = dataframe
        self.columns = columns
        self.units = units
        self.quantity = _transform_into_quantity(dataframe, units)


class SpatialTimeData:
    def __init__(self, location, science=None, interpolation_method='linear', data_cadence=None, combine_method='mean', desired_units=None, desired_coord=None, desired_coord_kwargs={}):
        self.location = location
        self.science = science
        self.interpolation_method = interpolation_method
        self._combine_data(interpolation_method, data_cadence, combine_method, desired_coord, desired_units, desired_coord_kwargs)
        

    def _combine_data(self, interpolation_method, data_cadence, combine_method, desired_coord, desired_units, desired_coord_kwargs):
        """
        Prepare data read from pyspedas for spatiotemporal plotting.

        Inputs
        location_data: iterable of spatial data
        science_data: iterable or list of iterables of spatial data
        location_interpolation_method (str): interpolation method at which to sample location data to match science data, if applicable.
        data_cadence (str or None): time cadence at which to up or downsample final dataset.
        science_columns_prefix (str): prefix to apply to science_data columns in final dataset. Will be superceded if science_columns_names exist.
        science_columns_names (list): list of names of science_data columns to name in the final dataset. Must match shape of science_data
        location_columns_prefix (str): prefix to apply to location_data columns in final dataset. Will be superceded if location_columns_names exist.
        location_columns_names (list): list of names of location_data columns to name in the final dataset. Must match shape of location_data.
        combine_method (str): name of method of how to join merged dataset. Must be 'mean' or 'median'.

        Output
        Merged dataset as a `pd.DataFrame` of location_data and science_data.
        """
        # check arguments
        #err, msg = _check_data_args(self.location, self.science)
        #if (err is not None) and (msg is not None):
        #    return _print_error_msg(err, msg)
        self.location = _convert_coordinates(self.location, desired_coord=desired_coord, coord_kwargs=desired_coord_kwargs)
        if self.science is not None:
            self.science = _convert_units(self.science, desired_units=desired_units)
        merged_dataframe = _match_data_cadence(self.location, self.science, interpolation_method=interpolation_method, data_cadence=data_cadence, combine_method=combine_method)
        self.data = merged_dataframe
        self.units = desired_units
        self.coord = desired_coord
        return self
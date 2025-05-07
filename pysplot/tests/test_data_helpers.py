import pytest
import pandas as pd
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord, ICRS
from astropy.time import Time
from sunpy.coordinates import HeliocentricInertial, HeliocentricEarthEcliptic
from pysplot.io.data_helpers import (
    _transform_into_skycoord,
    _transform_into_quantity,
    _adjust_skycoord,
    _rename_dataframe_with_skycoord_columns,
    _convert_units,
    _convert_coordinates,
    _validate_data,
    _validate_coord,
    _validate_units,
    _validate_interpolation_method,
    _transform_into_dataframe,
    _check_column_length,
    match_data_cadence
)
from pysplot.io.data import SpatialData

@pytest.fixture
def sample_spatial_data():
    return [{'x': ['2020-01-01','2020-01-02', '2020-01-03'], 'y': [[1, 2, 3], [4, 5, 6], [7, 8, 9]]}]

@pytest.fixture
def sample_science_data():
    return [{'x': ['2020-01-01 01:00','2020-01-02 03:00', '2020-01-03 02:00'], 'y': [[10, 20], [30, 40], [50, 60]]}]


def test_transform_into_skycoord_3d():
    df = pd.DataFrame([[1, 2, 3], [4, 5, 6]], columns=["x", "y", "z"], index=pd.date_range("2022-01-01", periods=2))
    coord = _transform_into_skycoord(df, ICRS, u.deg)
    assert isinstance(coord, ICRS)
    assert coord.shape[0] == 2

def test_transform_into_skycoord_2d():
    df = pd.DataFrame([[1, 2], [3, 4]], columns=["lon", "lat"], index=pd.date_range("2020", periods=2))
    coord = _transform_into_skycoord(df, ICRS, [u.deg, u.deg])
    assert isinstance(coord, ICRS)

def test_transform_into_quantity():
    df = pd.DataFrame([[1, 2], [3, 4]], columns=["a", "b"])
    q = _transform_into_quantity(df, u.km)
    assert q.shape == df.shape
    assert isinstance(q, u.Quantity)

def test_adjust_skycoord():
    obstime = Time('2025-03-25')
    sc = SkyCoord(1, 2, 3, unit=u.AU, frame=HeliocentricInertial(representation_type='cartesian',  obstime=obstime) )
    transformed, frame = _adjust_skycoord(sc, HeliocentricEarthEcliptic, {}, obstime)
    assert isinstance(transformed, SkyCoord)
    assert isinstance(frame, HeliocentricEarthEcliptic)

def test_rename_dataframe_with_skycoord_columns():
    df = pd.DataFrame([[1, 2, 3], [4, 5, 6]], columns=['x', 'y', 'z'], index=pd.date_range("2022-01-01", periods=2))
    frame = ICRS(representation_type='cartesian')
    sc = SkyCoord(1, 2, 3, unit=u.kpc, frame=frame)
    df, cols = _rename_dataframe_with_skycoord_columns(df.copy(), sc, frame, ['x', 'y', 'z'], u.AU)
    assert all(col in df.columns for col in cols)

def test_convert_coordinates(sample_spatial_data):
    spatial_obj = SpatialData(sample_spatial_data, coord=ICRS, units=u.deg)
    df = _convert_coordinates(spatial_obj, desired_coord=ICRS, desired_units=u.AU)
    assert isinstance(df, pd.DataFrame)
    assert df.shape[0] == spatial_obj.data.shape[0]


def test_convert_units():
    class MockScience:
        quantity = np.array([[1, 2], [3, 4]]) * u.kg
        units = u.kg
        data = pd.DataFrame([[1, 2], [3, 4]], columns=["a", "b"])
        columns = ["a", "b"]

    science = MockScience()
    result = _convert_units(science, desired_units=u.g)
    assert isinstance(result, pd.DataFrame)
    assert np.allclose(result.values, science.data.values * 1000)

def test_transform_into_dataframe(sample_spatial_data):
    df, cols = _transform_into_dataframe(sample_spatial_data, column_names=[['x','y','z']])
    assert isinstance(df, pd.DataFrame)
    assert list(df.columns) == cols
    assert df.shape == (3, 3)

def test_transform_into_dataframe_columns():
    list_of_data = [
        {'x': ['2023-01-01', '2023-01-02'], 'y': [[1, 2], [3, 4]]},
        {'x': ['2023-01-01', '2023-01-02'], 'y': [[5, 6], [7, 8]]}
    ]
    column_names = [['A', 'B'], ['C', 'D']]
    df, cols = _transform_into_dataframe(list_of_data, column_names, combine_axis='columns')
    assert list(df.columns) == cols
    assert df.shape == (2, 4)

def test_validate_data_valid():
    valid_data = {'x': [1, 2], 'y': [3, 4]}
    _validate_data(valid_data)  # should not raise

def test_validate_data_invalid():
    with pytest.raises(TypeError):
        _validate_data([1, 2])
    with pytest.raises(KeyError):
        _validate_data({'y': [1, 2]})
    with pytest.raises(KeyError):
        _validate_data({'x': [1, 2]})

def test_validate_interpolation_method_valid():
    _validate_interpolation_method('linear')  # Should pass

def test_validate_interpolation_method_invalid():
    with pytest.raises(ValueError):
        _validate_interpolation_method('invalid')

def test_validate_units():
    assert _validate_units(u.km) is None
    assert _validate_units(None) is None

def test_validate_coord():
    coord = ICRS
    assert _validate_coord(coord) is None

def test_check_column_length():
    data = [{'x': [1,1], 'y': [[1, 2],[3,4]]}]
    column_names = ['foo']
    output = _check_column_length(data, column_names, prefix='test')
    assert output[0][0].startswith('test_')

    column_names = [['A', 'B']]
    result = _check_column_length(data, column_names, prefix='test')
    assert result == column_names

    # no column names provided
    result = _check_column_length(data, [[]], prefix='test')
    assert isinstance(result, list)
    assert len(result[0]) == 2

def test_match_data_cadence(sample_spatial_data, sample_science_data):
    df1, _ = _transform_into_dataframe(sample_spatial_data)
    df2, _ = _transform_into_dataframe(sample_science_data)
    df1.index = pd.date_range("2022-01-01", periods=3, freq='1min')
    df2.index = pd.date_range("2022-01-01", periods=3, freq='2min')
    matched = match_data_cadence(df1, df2, interpolation_method='linear', data_cadence='1min')
    assert isinstance(matched, pd.DataFrame)
    assert len(matched) >= 3


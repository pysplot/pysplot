import pytest
import pandas as pd
import astropy.units as u
from sunpy.coordinates import HeliocentricInertial, Heliocentric, GeocentricSolarEcliptic, GeocentricEarthEquatorial
from pysplot.io.data import SpatialData, ScienceData, SpatialTimeData


@pytest.fixture
def sample_spatial_data():
    return [{'x': [0, 1, 2], 'y': [[1, 2, 3], [4, 5, 6], [7, 8, 9]]}]


@pytest.fixture
def sample_science_data():
    return [{'x': [0, 1, 2], 'y': [[10, 20], [30, 40], [50, 60]]}]


def test_spatial_data_init(sample_spatial_data):
    obj = SpatialData(
        spatial=sample_spatial_data,
        coord=HeliocentricInertial,
        units=[u.deg, u.deg, u.AU],
        spatial_columns_names=[['lon', 'lat', 'distance']]
    )

    assert isinstance(obj.data, pd.DataFrame)
    assert obj.coord is HeliocentricInertial
    assert obj.units == [u.deg, u.deg, u.AU]
    assert 'lon' in obj.data.columns 
    assert hasattr(obj, 'skycoord')


def test_science_data_init(sample_science_data):
    obj = ScienceData(
        science=sample_science_data,
        units=u.km/u.s,
        science_columns_names=[['Vx', 'Vy']]
    )

    assert isinstance(obj.data, pd.DataFrame)
    assert obj.units == u.km/u.s
    assert hasattr(obj, 'quantity')
    assert 'Vx' in obj.data.columns


def test_spatial_time_data_init(sample_spatial_data, sample_science_data):
    spatial_obj = SpatialData(
        spatial=sample_spatial_data,
        coord=HeliocentricInertial,
        units=u.AU,
        spatial_columns_names=[['X', 'Y', 'Z']],
        coord_kwargs = {'representation_type':'cartesian'}
    )

    science_obj = ScienceData(
        science=sample_science_data,
        units=u.km/u.s,
        science_columns_names=[['Vx', 'Vy']]
    )

    spt_data = SpatialTimeData(
        spatial=spatial_obj,
        science=science_obj,
        interpolation_method='linear',
        data_cadence=None,
        combine_method='mean',
        desired_coord=Heliocentric,
        desired_science_units=u.m/u.s,
        desired_spatial_units=[u.deg, u.deg, u.km],
        desired_coord_kwargs = {'observer':'Earth', 'representation_type':'spherical'}
    )

    assert isinstance(spt_data.data, pd.DataFrame)
    assert spt_data.coord == Heliocentric
    assert hasattr(spt_data, 'skycoord')
    assert 'Vx' in spt_data.data.columns


def test_transform_coord_method(sample_spatial_data, sample_science_data):
    spatial_obj = SpatialData(
        spatial=sample_spatial_data,
        coord=GeocentricSolarEcliptic,
        units=u.AU,
        spatial_columns_names=[['X', 'Y', 'Z']],
        coord_kwargs={'representation_type':'cartesian'}
    )
    science_obj = ScienceData(
        science=sample_science_data,
        units=u.km/u.s,
        science_columns_names=[['Vx', 'Vy']]
    )

    spt_data = SpatialTimeData(
        spatial=spatial_obj,
        science=science_obj,
        desired_coord=GeocentricSolarEcliptic,
        desired_spatial_units=u.km,
        desired_science_units=u.km/u.s,
        desired_coord_kwargs={'representation_type':'cartesian'}
    )

    new_data = spt_data.transform_coord(GeocentricEarthEquatorial)
    assert isinstance(new_data.data, pd.DataFrame)
    assert hasattr(new_data, 'skycoord')
    assert new_data.coord == GeocentricEarthEquatorial
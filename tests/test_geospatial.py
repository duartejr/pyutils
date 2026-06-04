"""Unit tests for pyutils.geospatial module."""

import numpy as np
import pytest
from pyutils.geospatial import (
    InverseDistanceWeighting,
    ThiessenPolygon,
    XarrayNetCDFHandler,
    ShapefileHandler,
    MapRenderer,
)
from pyutils.core import ValidationError


# ------------------------------------------------------------------ #
# Fixtures                                                             #
# ------------------------------------------------------------------ #

@pytest.fixture
def three_stations():
    return {
        "x": np.array([0.0, 10.0, 20.0]),
        "y": np.array([0.0, 10.0, 20.0]),
        "z": np.array([100.0, 150.0, 120.0]),
    }


@pytest.fixture
def target_grid():
    return {
        "xi": np.linspace(0, 20, 5),
        "yi": np.linspace(0, 20, 5),
    }


# ------------------------------------------------------------------ #
# InverseDistanceWeighting                                             #
# ------------------------------------------------------------------ #

class TestIDW:
    def test_output_shape(self, three_stations, target_grid):
        idw = InverseDistanceWeighting()
        result = idw.interpolate(**three_stations, **target_grid)
        assert result.shape == (5, 5)

    def test_exact_hit_returns_station_value(self):
        idw = InverseDistanceWeighting()
        x = np.array([0.0, 10.0])
        y = np.array([0.0,  0.0])
        z = np.array([100.0, 200.0])
        result = idw.interpolate(x, y, z, xi=np.array([0.0]), yi=np.array([0.0]))
        assert abs(result[0, 0] - 100.0) < 1e-6

    def test_interpolated_values_within_range(self, three_stations, target_grid):
        idw = InverseDistanceWeighting()
        result = idw.interpolate(**three_stations, **target_grid)
        assert np.all(result >= three_stations["z"].min())
        assert np.all(result <= three_stations["z"].max())

    def test_higher_power_increases_locality(self, three_stations, target_grid):
        result_p1 = InverseDistanceWeighting(power=1).interpolate(**three_stations, **target_grid)
        result_p4 = InverseDistanceWeighting(power=4).interpolate(**three_stations, **target_grid)
        # Higher power → sharper transitions → greater variance
        assert np.std(result_p4) >= np.std(result_p1)

    def test_mismatched_xyz_raises(self, target_grid):
        idw = InverseDistanceWeighting()
        with pytest.raises(ValidationError):
            idw.interpolate(
                x=np.array([0.0, 1.0]),
                y=np.array([0.0]),          # length mismatch
                z=np.array([1.0, 2.0]),
                **target_grid,
            )

    def test_invalid_power_raises(self):
        with pytest.raises(ValidationError):
            InverseDistanceWeighting(power=0)

    def test_validate_returns_true(self):
        assert InverseDistanceWeighting().validate() is True


# ------------------------------------------------------------------ #
# ThiessenPolygon                                                      #
# ------------------------------------------------------------------ #

class TestThiessen:
    def test_output_shape(self, three_stations, target_grid):
        th = ThiessenPolygon()
        result = th.interpolate(**three_stations, **target_grid)
        assert result.shape == (5, 5)

    def test_only_station_values_returned(self, three_stations, target_grid):
        th = ThiessenPolygon()
        result = th.interpolate(**three_stations, **target_grid)
        unique_vals = set(np.unique(result))
        station_vals = set(three_stations["z"])
        assert unique_vals.issubset(station_vals)

    def test_nearest_station_assigned(self):
        th = ThiessenPolygon()
        x = np.array([0.0, 100.0])
        y = np.array([0.0,   0.0])
        z = np.array([10.0,  20.0])
        result = th.interpolate(x, y, z, xi=np.array([5.0]), yi=np.array([0.0]))
        assert result[0, 0] == 10.0   # closer to first station

    def test_validate_returns_true(self):
        assert ThiessenPolygon().validate() is True


# ------------------------------------------------------------------ #
# XarrayNetCDFHandler                                                  #
# ------------------------------------------------------------------ #

class TestXarrayNetCDFHandler:
    def test_instantiation(self):
        handler = XarrayNetCDFHandler()
        assert handler.name == "xarray_handler"

    def test_instantiation_with_chunks(self):
        handler = XarrayNetCDFHandler(chunks={"time": 50})
        assert handler.chunks == {"time": 50}

    def test_validate_without_dataset_is_false(self):
        handler = XarrayNetCDFHandler()
        assert handler.validate() is False

    def test_read_missing_file_raises(self, tmp_path):
        handler = XarrayNetCDFHandler()
        with pytest.raises(ValidationError):
            handler.read(str(tmp_path / "nonexistent.nc"))

    def test_write_without_dataset_raises(self, tmp_path):
        handler = XarrayNetCDFHandler()
        with pytest.raises(ValidationError):
            handler.write(str(tmp_path / "out.nc"))

    def test_read_write_roundtrip(self, tmp_path):
        """Write a minimal NetCDF and read it back."""
        import xarray as xr
        ds = xr.Dataset({"temperature": (["time"], [15.0, 16.0, 17.0])})
        path = str(tmp_path / "test.nc")
        ds.to_netcdf(path)

        handler = XarrayNetCDFHandler()
        loaded = handler.read(path)
        assert "temperature" in loaded.data_vars
        assert handler.validate() is True

    def test_repr(self):
        assert "XarrayNetCDFHandler" in repr(XarrayNetCDFHandler())


# ------------------------------------------------------------------ #
# ShapefileHandler                                                     #
# ------------------------------------------------------------------ #

class TestShapefileHandler:
    def test_instantiation(self):
        shp = ShapefileHandler()
        assert shp.crs == "EPSG:4674"

    def test_custom_crs(self):
        shp = ShapefileHandler(crs="EPSG:4326")
        assert shp.crs == "EPSG:4326"

    def test_read_missing_file_raises(self):
        shp = ShapefileHandler()
        with pytest.raises(ValidationError):
            shp.read_shapefile("/nonexistent/path.shp")

    def test_write_without_geodataframe_raises(self, tmp_path):
        shp = ShapefileHandler()
        with pytest.raises(ValidationError):
            shp.write_shapefile(str(tmp_path / "out.shp"))

    def test_compute_areas_with_geodataframe(self):
        import geopandas as gpd
        from shapely.geometry import Polygon
        gdf = gpd.GeoDataFrame(
            geometry=[Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])],
            crs="EPSG:4326",
        )
        shp = ShapefileHandler()
        areas = shp.compute_areas(gdf)  # uses default EPSG:6933
        assert len(areas) == 1
        assert areas[0] > 0


# ------------------------------------------------------------------ #
# MapRenderer                                                          #
# ------------------------------------------------------------------ #

class TestMapRenderer:
    def test_instantiation(self):
        renderer = MapRenderer()
        assert renderer.figsize == (12, 8)

    def test_custom_figsize(self):
        renderer = MapRenderer(figsize=(10, 6))
        assert renderer.figsize == (10, 6)

    def test_create_figure_returns_fig_ax(self):
        import matplotlib
        matplotlib.use("Agg")
        renderer = MapRenderer()
        fig, ax = renderer.create_figure()
        assert fig is not None and ax is not None
        renderer.close()

    def test_set_title(self):
        import matplotlib
        matplotlib.use("Agg")
        renderer = MapRenderer()
        renderer.create_figure()
        renderer.set_title("Test Map")
        assert renderer.ax.get_title() == "Test Map"
        renderer.close()

    def test_close_clears_fig(self):
        import matplotlib
        matplotlib.use("Agg")
        renderer = MapRenderer()
        renderer.create_figure()
        renderer.close()
        assert renderer.fig is None
        assert renderer.ax is None

    def test_save_creates_file(self, tmp_path):
        import matplotlib
        matplotlib.use("Agg")
        renderer = MapRenderer()
        renderer.create_figure()
        path = str(tmp_path / "map.png")
        renderer.save(path)
        assert (tmp_path / "map.png").exists()
        renderer.close()

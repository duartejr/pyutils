"""Shapefile operations and geometry utilities.

Consolidates shp_area.py, buffer.py, and related shapefile operations
using geopandas and shapely for robust, modern geometry handling.
"""

from typing import Optional, Tuple, Union
from pathlib import Path
import numpy as np
import geopandas as gpd
import shapely.ops as ops
from shapely.geometry import Point, Polygon, MultiPolygon, shape
from shapely.ops import unary_union
import pyproj
from functools import partial

from pyutils.core import ValidationError


class ShapefileHandler:
    """Handle shapefile I/O and spatial operations.

    Consolidates shp_area.py and buffer.py with modern geopandas interface.
    Enables robust geometry operations, coordinate transformation, and spatial analysis.
    """

    def __init__(self, crs: str = "EPSG:4674"):
        """Initialize shapefile handler.

        Parameters
        ----------
        crs : str, optional
            Default coordinate reference system. Default is EPSG:4674 (SIRGAS 2000 - Brazil).
        """
        self.crs = crs
        self.gdf: Optional[gpd.GeoDataFrame] = None

    def read_shapefile(self, path: str) -> gpd.GeoDataFrame:
        """Read shapefile using geopandas.

        Parameters
        ----------
        path : str
            Path to shapefile (.shp file).

        Returns
        -------
        gpd.GeoDataFrame
            Loaded geodataframe with all attributes and geometries.

        Raises
        ------
        ValidationError
            If file not found or cannot be read.
        """
        path = Path(path)
        if not path.exists():
            raise ValidationError(f"Shapefile not found: {path}")

        try:
            self.gdf = gpd.read_file(path)
            return self.gdf
        except Exception as e:
            raise ValidationError(f"Failed to read shapefile: {e}")

    def write_shapefile(
        self,
        path: str,
        gdf: Optional[gpd.GeoDataFrame] = None,
    ) -> None:
        """Write geodataframe to shapefile.

        Parameters
        ----------
        path : str
            Output file path (.shp).
        gdf : gpd.GeoDataFrame, optional
            GeoDataFrame to write. If None, uses current gdf.

        Raises
        ------
        ValidationError
            If no geodataframe available or write fails.
        """
        if gdf is None:
            if self.gdf is None:
                raise ValidationError("No geodataframe loaded")
            gdf = self.gdf

        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)

        try:
            gdf.to_file(path, driver="ESRI Shapefile")
        except Exception as e:
            raise ValidationError(f"Failed to write shapefile: {e}")

    def compute_area(
        self,
        geometry: Optional[Polygon] = None,
        to_crs: str = "EPSG:6933",
    ) -> float:
        """Compute polygon area in projected coordinate system.

        Projects geometry to equal-area projection (AEA - Albers Equal Area)
        before computing area to get accurate results.

        Parameters
        ----------
        geometry : Polygon, optional
            Geometry to compute area. If None, uses first geometry in current gdf.
        to_crs : str, optional
            Target CRS for projection. Default is EPSG:6933 (WGS 84 / NSIDC EASE-Grid 2.0 Global — equal-area).

        Returns
        -------
        float
            Area in square meters.

        Raises
        ------
        ValidationError
            If no geometry available or computation fails.
        """
        if geometry is None:
            if self.gdf is None or len(self.gdf) == 0:
                raise ValidationError("No geometry available")
            geometry = self.gdf.geometry.iloc[0]

        try:
            # Create temporary GeoDataFrame for projection
            gdf_temp = gpd.GeoDataFrame(
                geometry=[geometry], crs=self.crs
            )

            # Project to equal-area CRS
            gdf_proj = gdf_temp.to_crs(to_crs)

            return gdf_proj.geometry.iloc[0].area
        except Exception as e:
            raise ValidationError(f"Failed to compute area: {e}")

    def compute_areas(
        self,
        gdf: Optional[gpd.GeoDataFrame] = None,
        to_crs: str = "EPSG:6933",
    ) -> np.ndarray:
        """Compute areas for all geometries in geodataframe.

        Parameters
        ----------
        gdf : gpd.GeoDataFrame, optional
            GeoDataFrame to compute areas. If None, uses current gdf.
        to_crs : str, optional
            Target CRS for projection.

        Returns
        -------
        np.ndarray
            Areas in square meters for each geometry.

        Raises
        ------
        ValidationError
            If no geodataframe available.
        """
        if gdf is None:
            if self.gdf is None:
                raise ValidationError("No geodataframe loaded")
            gdf = self.gdf

        try:
            gdf_proj = gdf.to_crs(to_crs)
            return gdf_proj.geometry.area.values
        except Exception as e:
            raise ValidationError(f"Failed to compute areas: {e}")

    def buffer_geometries(
        self,
        distance: float,
        gdf: Optional[gpd.GeoDataFrame] = None,
        output_path: Optional[str] = None,
    ) -> gpd.GeoDataFrame:
        """Create buffers around geometries.

        Parameters
        ----------
        distance : float
            Buffer distance in current CRS units.
        gdf : gpd.GeoDataFrame, optional
            GeoDataFrame to buffer. If None, uses current gdf.
        output_path : str, optional
            Path to save buffered shapefile.

        Returns
        -------
        gpd.GeoDataFrame
            New GeoDataFrame with buffered geometries.

        Raises
        ------
        ValidationError
            If no geodataframe available or operation fails.
        """
        if gdf is None:
            if self.gdf is None:
                raise ValidationError("No geodataframe loaded")
            gdf = self.gdf.copy()
        else:
            gdf = gdf.copy()

        try:
            gdf["geometry"] = gdf.geometry.buffer(distance)

            if output_path:
                self.write_shapefile(output_path, gdf)

            return gdf
        except Exception as e:
            raise ValidationError(f"Failed to buffer geometries: {e}")

    def spatial_join(
        self,
        gdf1: gpd.GeoDataFrame,
        gdf2: gpd.GeoDataFrame,
        how: str = "left",
        predicate: str = "intersects",
    ) -> gpd.GeoDataFrame:
        """Perform spatial join between two geodataframes.

        Parameters
        ----------
        gdf1 : gpd.GeoDataFrame
            Left geodataframe.
        gdf2 : gpd.GeoDataFrame
            Right geodataframe.
        how : str, optional
            Join type ('left', 'right', 'inner', 'outer').
        predicate : str, optional
            Spatial relationship ('intersects', 'contains', 'within', etc.).

        Returns
        -------
        gpd.GeoDataFrame
            Joined geodataframe.

        Raises
        ------
        ValidationError
            If join operation fails.
        """
        try:
            return gpd.sjoin(
                gdf1, gdf2, how=how, predicate=predicate
            )
        except Exception as e:
            raise ValidationError(f"Spatial join failed: {e}")

    def dissolve(
        self,
        by: Optional[str] = None,
        gdf: Optional[gpd.GeoDataFrame] = None,
        aggfunc: str = "first",
    ) -> gpd.GeoDataFrame:
        """Dissolve geometries based on attribute.

        Parameters
        ----------
        by : str, optional
            Column name to dissolve by. If None, dissolves all geometries.
        gdf : gpd.GeoDataFrame, optional
            GeoDataFrame to dissolve. If None, uses current gdf.
        aggfunc : str, optional
            Aggregation function for other columns.

        Returns
        -------
        gpd.GeoDataFrame
            Dissolved geodataframe.

        Raises
        ------
        ValidationError
            If dissolve fails.
        """
        if gdf is None:
            if self.gdf is None:
                raise ValidationError("No geodataframe loaded")
            gdf = self.gdf

        try:
            return gdf.dissolve(by=by, aggfunc=aggfunc)
        except Exception as e:
            raise ValidationError(f"Dissolve failed: {e}")

    def clip(
        self,
        gdf_to_clip: gpd.GeoDataFrame,
        clip_geom: Union[Polygon, gpd.GeoDataFrame],
    ) -> gpd.GeoDataFrame:
        """Clip geodataframe to geometry or another geodataframe.

        Parameters
        ----------
        gdf_to_clip : gpd.GeoDataFrame
            GeoDataFrame to clip.
        clip_geom : Polygon or gpd.GeoDataFrame
            Clipping geometry or geodataframe.

        Returns
        -------
        gpd.GeoDataFrame
            Clipped geodataframe.

        Raises
        ------
        ValidationError
            If clipping fails.
        """
        try:
            if isinstance(clip_geom, gpd.GeoDataFrame):
                clip_geom = unary_union(clip_geom.geometry)

            return gpd.clip(gdf_to_clip, clip_geom)
        except Exception as e:
            raise ValidationError(f"Clipping failed: {e}")

    def reproject(
        self,
        to_crs: str,
        gdf: Optional[gpd.GeoDataFrame] = None,
    ) -> gpd.GeoDataFrame:
        """Reproject geodataframe to different CRS.

        Parameters
        ----------
        to_crs : str
            Target CRS (e.g., 'EPSG:4326').
        gdf : gpd.GeoDataFrame, optional
            GeoDataFrame to reproject. If None, uses current gdf.

        Returns
        -------
        gpd.GeoDataFrame
            Reprojected geodataframe.

        Raises
        ------
        ValidationError
            If reprojection fails.
        """
        if gdf is None:
            if self.gdf is None:
                raise ValidationError("No geodataframe loaded")
            gdf = self.gdf

        try:
            return gdf.to_crs(to_crs)
        except Exception as e:
            raise ValidationError(f"Reprojection failed: {e}")

    def get_vertices(
        self,
        attribute: Optional[str] = None,
        value: Optional[str] = None,
        gdf: Optional[gpd.GeoDataFrame] = None,
    ) -> np.ndarray:
        """Extract polygon vertices as a numpy array.

        Parameters
        ----------
        attribute : str, optional
            Column name to filter by (e.g. 'basin'). If None, uses first feature.
        value : str, optional
            Value to match in the attribute column.
        gdf : gpd.GeoDataFrame, optional
            GeoDataFrame to extract from. If None, uses current gdf.

        Returns
        -------
        np.ndarray
            Array of shape (N, 2) with (x, y) coordinates of polygon vertices.

        Raises
        ------
        ValidationError
            If no geometry is available or the attribute/value is not found.

        Examples
        --------
        >>> shp = ShapefileHandler()
        >>> shp.read_shapefile("basin.shp")
        >>> vertices = shp.get_vertices()
        >>> vertices.shape
        (243, 2)
        """
        if gdf is None:
            if self.gdf is None:
                raise ValidationError("No geodataframe loaded")
            gdf = self.gdf

        try:
            if attribute is not None and value is not None:
                subset = gdf[gdf[attribute] == value]
                if len(subset) == 0:
                    raise ValidationError(
                        f"No feature with {attribute}='{value}'"
                    )
                geom = subset.geometry.iloc[0]
            else:
                geom = gdf.geometry.iloc[0]

            if hasattr(geom, "geoms"):
                # MultiPolygon — concatenate all rings
                vertices = np.concatenate(
                    [np.array(part.exterior.coords) for part in geom.geoms]
                )
            else:
                vertices = np.array(geom.exterior.coords)

            return vertices
        except ValidationError:
            raise
        except Exception as e:
            raise ValidationError(f"Failed to extract vertices: {e}")

    def points_in_polygon(
        self,
        lat: np.ndarray,
        lon: np.ndarray,
        buffer: float = 0.0,
        gdf: Optional[gpd.GeoDataFrame] = None,
    ) -> np.ndarray:
        """Return grid coordinates that fall inside the loaded shapefile.

        Parameters
        ----------
        lat : np.ndarray
            1-D array of latitude values.
        lon : np.ndarray
            1-D array of longitude values.
        buffer : float, optional
            Buffer distance in CRS units applied to geometries before
            the point-in-polygon test. Default is 0 (no buffer).
        gdf : gpd.GeoDataFrame, optional
            GeoDataFrame to test against. If None, uses current gdf.

        Returns
        -------
        np.ndarray
            Array of shape (M, 2) with (lon, lat) of points inside the polygon,
            where M ≤ len(lat) × len(lon).

        Raises
        ------
        ValidationError
            If no geodataframe available or operation fails.

        Examples
        --------
        >>> shp = ShapefileHandler()
        >>> shp.read_shapefile("basin.shp")
        >>> inside = shp.points_in_polygon(lat=lats, lon=lons, buffer=0.01)
        """
        if gdf is None:
            if self.gdf is None:
                raise ValidationError("No geodataframe loaded")
            gdf = self.gdf

        try:
            geom = gdf.geometry
            if buffer > 0:
                geom = geom.buffer(buffer)
            merged = unary_union(geom)

            xi, yi = np.meshgrid(lon, lat)
            flat_lon = xi.ravel()
            flat_lat = yi.ravel()

            points_gdf = gpd.GeoDataFrame(
                geometry=gpd.points_from_xy(flat_lon, flat_lat),
                crs=self.crs,
            )

            inside_mask = points_gdf.geometry.within(merged)
            inside_coords = np.column_stack(
                [flat_lon[inside_mask], flat_lat[inside_mask]]
            )

            return inside_coords
        except Exception as e:
            raise ValidationError(f"points_in_polygon failed: {e}")

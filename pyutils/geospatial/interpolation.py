"""Spatial interpolation methods.

Consolidates idw.py and thiessen.py implementations with improved
performance and consistency.
"""

from typing import Tuple, Optional
import numpy as np
from scipy.spatial import Voronoi
from pyutils.core import SpatialInterpolator, ValidationError


class InverseDistanceWeighting(SpatialInterpolator):
    """Inverse Distance Weighting (IDW) spatial interpolator.

    Uses weighted average of nearby observations, with weights inversely
    proportional to distance. Suitable for meteorological station data.
    Consolidates and improves previous idw.py implementation.

    References
    ----------
    Shepard, D. 1968. A two-dimensional interpolation function for
    irregularly-spaced data. Proceedings of the 23rd National Conference
    of the Association for Computing Machinery, pp. 517-524.
    """

    def __init__(self, power: float = 2.0, min_neighbors: int = 1):
        """Initialize IDW interpolator.

        Parameters
        ----------
        power : float, optional
            Power parameter (typically 1-3). Default is 2 (inverse square distance).
        min_neighbors : int, optional
            Minimum neighbors to use. Default is 1.

        Raises
        ------
        ValidationError
            If parameters are invalid.
        """
        super().__init__(
            name="Inverse Distance Weighting",
            description=f"IDW interpolation with power={power}",
        )
        if power <= 0:
            raise ValidationError("Power must be positive")
        if min_neighbors < 1:
            raise ValidationError("min_neighbors must be >= 1")

        self.power = power
        self.min_neighbors = min_neighbors

    def interpolate(
        self,
        x: np.ndarray,
        y: np.ndarray,
        z: np.ndarray,
        xi: np.ndarray,
        yi: np.ndarray,
    ) -> np.ndarray:
        """Interpolate values using IDW.

        Parameters
        ----------
        x, y : np.ndarray
            Coordinates of source points (UTM or lat/lon).
        z : np.ndarray
            Values at source points.
        xi, yi : np.ndarray
            Target grid coordinates (1D arrays for mesh, or 2D for irregular).

        Returns
        -------
        np.ndarray
            Interpolated values at target locations (same shape as xi).

        Raises
        ------
        ValidationError
            If input data is invalid.
        """
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        z = np.asarray(z, dtype=float)
        xi = np.asarray(xi, dtype=float)
        yi = np.asarray(yi, dtype=float)

        if not (len(x) == len(y) == len(z)):
            raise ValidationError("x, y, z must have same length")

        if len(z) < self.min_neighbors:
            raise ValidationError(
                f"Need at least {self.min_neighbors} points"
            )

        # Handle 1D grids (create mesh)
        if xi.ndim == 1 and yi.ndim == 1:
            xi_mesh, yi_mesh = np.meshgrid(xi, yi)
            original_shape = xi_mesh.shape
        else:
            xi_mesh = xi.ravel()
            yi_mesh = yi.ravel()
            original_shape = xi.shape

        zi = np.zeros(xi_mesh.ravel().shape)

        # Vectorized distance calculation
        for i, (xi_val, yi_val) in enumerate(
            zip(xi_mesh.ravel(), yi_mesh.ravel())
        ):
            distances = np.sqrt((x - xi_val) ** 2 + (y - yi_val) ** 2)

            # Check for exact hits
            exact_hit = distances < 1e-10
            if np.any(exact_hit):
                zi[i] = np.mean(z[exact_hit])
            else:
                # Use inverse distance weighting
                weights = 1.0 / (distances ** self.power)
                weights /= np.sum(weights)
                zi[i] = np.sum(z * weights)

        return zi.reshape(original_shape)

    def validate(self) -> bool:
        """Validate interpolator parameters."""
        return self.power > 0 and self.min_neighbors >= 1


class ThiessenPolygon(SpatialInterpolator):
    """Thiessen polygon (Voronoi) spatial interpolator.

    Assigns each grid point to nearest station. Useful for precipitation
    estimation and areal rainfall computation. Consolidates previous thiessen.py.

    References
    ----------
    Thiessen, A.H. 1911. Precipitation averages for large areas.
    Monthly Weather Review, 39(7), 1082-1089.
    """

    def __init__(self):
        """Initialize Thiessen polygon interpolator."""
        super().__init__(
            name="Thiessen Polygon",
            description="Voronoi-based nearest neighbor interpolation",
        )

    def interpolate(
        self,
        x: np.ndarray,
        y: np.ndarray,
        z: np.ndarray,
        xi: np.ndarray,
        yi: np.ndarray,
    ) -> np.ndarray:
        """Interpolate using Thiessen polygons.

        Parameters
        ----------
        x, y : np.ndarray
            Station coordinates.
        z : np.ndarray
            Values at stations.
        xi, yi : np.ndarray
            Target grid coordinates.

        Returns
        -------
        np.ndarray
            Values assigned from nearest station.

        Raises
        ------
        ValidationError
            If input data is invalid.
        """
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        z = np.asarray(z, dtype=float)
        xi = np.asarray(xi, dtype=float)
        yi = np.asarray(yi, dtype=float)

        if not (len(x) == len(y) == len(z)):
            raise ValidationError("x, y, z must have same length")

        # Handle 1D grids
        if xi.ndim == 1 and yi.ndim == 1:
            xi_mesh, yi_mesh = np.meshgrid(xi, yi)
            original_shape = xi_mesh.shape
        else:
            xi_mesh = xi.ravel()
            yi_mesh = yi.ravel()
            original_shape = xi.shape

        # Stack coordinates for distance calculation
        points = np.column_stack([x, y])
        grid = np.column_stack([xi_mesh.ravel(), yi_mesh.ravel()])

        # Find nearest neighbor
        distances = np.sqrt(((grid[:, np.newaxis] - points) ** 2).sum(axis=2))
        nearest_idx = np.argmin(distances, axis=1)
        zi = z[nearest_idx]

        return zi.reshape(original_shape)

    def compute_polygon_areas(
        self,
        x: np.ndarray,
        y: np.ndarray,
        bounds: Optional[Tuple[float, float, float, float]] = None,
    ) -> np.ndarray:
        """Compute areas of Thiessen polygons.

        Parameters
        ----------
        x, y : np.ndarray
            Station coordinates.
        bounds : tuple, optional
            (xmin, xmax, ymin, ymax) for bounded area. If None, uses convex hull.

        Returns
        -------
        np.ndarray
            Area of each polygon.

        Raises
        ------
        ValidationError
            If computation fails.
        """
        try:
            points = np.column_stack([x, y])
            vor = Voronoi(points)

            areas = np.zeros(len(points))

            for point_idx, region_idx in enumerate(vor.point_region):
                region = vor.regions[region_idx]

                if -1 in region:
                    # Infinite region - skip or estimate
                    areas[point_idx] = np.inf
                else:
                    vertices = vor.vertices[region]
                    if len(vertices) > 2:
                        # Shoelace formula for polygon area
                        areas[point_idx] = 0.5 * np.abs(
                            np.dot(
                                vertices[:, 0],
                                np.roll(vertices[:, 1], 1)
                            )
                            - np.dot(
                                vertices[:, 1],
                                np.roll(vertices[:, 0], 1)
                            )
                        )

            return areas
        except Exception as e:
            raise ValidationError(f"Failed to compute polygon areas: {e}")

    def validate(self) -> bool:
        """Validate interpolator state."""
        return True

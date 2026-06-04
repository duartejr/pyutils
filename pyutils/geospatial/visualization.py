"""Cartographic visualization utilities.

Consolidates plot_maps.py and custom_colorbar.py functionality for
creating thematic maps and customizing colorbars.
"""

from typing import Optional, List, Tuple, Union
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import ListedColormap, BoundaryNorm
import geopandas as gpd


class MapRenderer:
    """Create thematic maps with custom styling and colorbars.

    Consolidates plot_maps.py and custom_colorbar.py for integrated
    cartographic visualization with matplotlib and geopandas.
    """

    def __init__(self, figsize: Tuple[float, float] = (12, 8)):
        """Initialize map renderer.

        Parameters
        ----------
        figsize : tuple, optional
            Figure size (width, height) in inches. Default is (12, 8).
        """
        self.figsize = figsize
        self.fig = None
        self.ax = None

    def create_figure(self, figsize: Optional[Tuple[float, float]] = None) -> Tuple:
        """Create matplotlib figure and axis.

        Parameters
        ----------
        figsize : tuple, optional
            Figure size. If None, uses instance figsize.

        Returns
        -------
        tuple
            (figure, axis) objects.
        """
        figsize = figsize or self.figsize
        self.fig, self.ax = plt.subplots(figsize=figsize)
        return self.fig, self.ax

    def plot_geodataframe(
        self,
        gdf: gpd.GeoDataFrame,
        column: Optional[str] = None,
        cmap: str = "viridis",
        edgecolor: str = "k",
        alpha: float = 0.8,
        legend: bool = True,
        ax: Optional[plt.Axes] = None,
    ) -> plt.Axes:
        """Plot geodataframe with optional column-based coloring.

        Parameters
        ----------
        gdf : gpd.GeoDataFrame
            GeoDataFrame to plot.
        column : str, optional
            Column name for color mapping.
        cmap : str, optional
            Colormap name. Default is 'viridis'.
        edgecolor : str, optional
            Edge color for geometries.
        alpha : float, optional
            Transparency (0-1).
        legend : bool, optional
            Show legend.
        ax : plt.Axes, optional
            Axes to plot on. If None, creates new figure.

        Returns
        -------
        plt.Axes
            Axes with plotted data.
        """
        if ax is None:
            self.create_figure()
            ax = self.ax

        gdf.plot(
            ax=ax,
            column=column,
            cmap=cmap,
            edgecolor=edgecolor,
            alpha=alpha,
            legend=legend,
        )

        ax.set_aspect("equal")
        ax.grid(alpha=0.3)

        return ax

    def add_colorbar(
        self,
        cmap: Union[str, ListedColormap],
        bounds: Optional[List[float]] = None,
        label: str = "",
        orientation: str = "vertical",
        shrink: float = 0.8,
        extend: str = "neither",
    ) -> None:
        """Add custom colorbar to figure.

        Parameters
        ----------
        cmap : str or ListedColormap
            Colormap to use.
        bounds : list, optional
            Color boundaries for discrete colorbar.
        label : str, optional
            Colorbar label.
        orientation : str, optional
            'vertical' or 'horizontal'.
        shrink : float, optional
            Shrink factor for colorbar (0-1).
        extend : str, optional
            Extend style ('neither', 'both', 'min', 'max').
        """
        if self.ax is None:
            raise ValueError("No axes available. Create figure first.")

        if isinstance(cmap, str):
            cmap = plt.get_cmap(cmap)

        if bounds is not None:
            norm = BoundaryNorm(bounds, cmap.N)
        else:
            norm = None

        cbar = plt.colorbar(
            mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
            ax=self.ax,
            orientation=orientation,
            shrink=shrink,
            extend=extend,
        )

        if label:
            cbar.set_label(label)

    def customize_colorbar(
        self,
        cmap: Union[str, ListedColormap],
        bounds: List[float],
        label: str = "",
        ticks_label: Optional[List[str]] = None,
        orientation: str = "vertical",
        position: Optional[List[float]] = None,
    ) -> None:
        """Create heavily customized colorbar (legacy custom_colorbar.py equivalent).

        Parameters
        ----------
        cmap : str or ListedColormap
            Colormap.
        bounds : list
            Color boundaries.
        label : str, optional
            Colorbar label.
        ticks_label : list, optional
            Custom tick labels.
        orientation : str, optional
            'vertical' or 'horizontal'.
        position : list, optional
            Colorbar position [left, bottom, width, height].
        """
        if self.fig is None:
            self.create_figure()

        if isinstance(cmap, str):
            cmap = plt.get_cmap(cmap)

        # Create custom axes for colorbar
        if position:
            cax = self.fig.add_axes(position)
        else:
            cax = self.fig.add_axes([0.92, 0.15, 0.02, 0.7])

        # Compute boundaries
        if len(bounds) > 1:
            step = abs(bounds[1] - bounds[0])
        else:
            step = 1

        boundaries = np.array([bounds[0] - step] + list(bounds) + [bounds[-1] + step])
        ticks = np.array(boundaries)

        norm = BoundaryNorm(bounds, cmap.N)
        cbar = mpl.colorbar.ColorbarBase(
            cax,
            cmap=cmap,
            norm=norm,
            boundaries=boundaries,
            extend="both",
            extendfrac="auto",
            ticks=ticks,
            spacing="uniform",
            orientation=orientation,
        )

        if ticks_label and len(ticks_label) > 0:
            if orientation == "vertical":
                cax.set_yticklabels(ticks_label)
            else:
                cax.set_xticklabels(ticks_label)

        if label:
            cbar.set_label(label)

    def set_title(self, title: str, fontsize: int = 14) -> None:
        """Set figure title.

        Parameters
        ----------
        title : str
            Title text.
        fontsize : int, optional
            Font size.
        """
        if self.ax is None:
            raise ValueError("No axes available")
        self.ax.set_title(title, fontsize=fontsize)

    def set_labels(
        self,
        xlabel: str = "",
        ylabel: str = "",
        fontsize: int = 12,
    ) -> None:
        """Set axis labels.

        Parameters
        ----------
        xlabel : str, optional
            X-axis label.
        ylabel : str, optional
            Y-axis label.
        fontsize : int, optional
            Font size.
        """
        if self.ax is None:
            raise ValueError("No axes available")
        if xlabel:
            self.ax.set_xlabel(xlabel, fontsize=fontsize)
        if ylabel:
            self.ax.set_ylabel(ylabel, fontsize=fontsize)

    def remove_axis(self) -> None:
        """Remove axis spines and ticks for clean map display."""
        if self.ax is None:
            raise ValueError("No axes available")
        self.ax.set_xticks([])
        self.ax.set_yticks([])
        for spine in self.ax.spines.values():
            spine.set_visible(False)

    def save(self, path: str, dpi: int = 300, bbox_inches: str = "tight") -> None:
        """Save figure to file.

        Parameters
        ----------
        path : str
            Output file path.
        dpi : int, optional
            Resolution in dots per inch.
        bbox_inches : str, optional
            Bounding box setting ('tight' for minimal padding).
        """
        if self.fig is None:
            raise ValueError("No figure to save")
        self.fig.savefig(path, dpi=dpi, bbox_inches=bbox_inches)

    def show(self) -> None:
        """Display figure."""
        if self.fig is None:
            raise ValueError("No figure to show")
        plt.show()

    def close(self) -> None:
        """Close figure."""
        if self.fig is not None:
            plt.close(self.fig)
            self.fig = None
            self.ax = None

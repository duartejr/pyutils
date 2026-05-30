"""
Setup configuration for pyutils.
"""
from setuptools import setup

# For flat-layout packages, list packages and modules separately
setup(
    packages=["pyeto", "rv"],
    py_modules=[
        "pyutils",  # Main package module
        "anom_pr", "bh_thorthwaite", "bh_thorthwaite2", "buffer", "ccsthorthwaite",
        "corr_med", "custom_colorbar", "eto", "evapotranspiration", "getvertshp",
        "idw", "isinpoly", "isinshp", "kriging", "plot_maps", "plot_maps_old",
        "read_nc", "rvGamma", "shp_area", "spi", "thiessen", "vazoes_minimas"
    ]
)

"""
Setup configuration for pyutils.
This file is for backward compatibility with older pip versions.
Primary configuration is in pyproject.toml.
"""
from setuptools import setup

setup(
    packages=["pyeto", "rv"],
    py_modules=[
        "bh_thorthwaite", "bh_thorthwaite2", "buffer", "ccsthorthwaite",
        "corr_med", "custom_colorbar", "eto", "evapotranspiration", "getvertshp",
        "idw", "isinpoly", "isinshp", "plot_maps",
        "read_nc", "rvGamma", "shp_area", "spi", "thiessen", "vazoes_minimas"
    ]
)

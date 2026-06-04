"""Physical and meteorological constants used across pyutils modules."""

import math

# Astronomical constants
SOLAR_CONSTANT = 0.0820  # MJ m^-2 min^-1 (solar irradiance at top of atmosphere)
STEFAN_BOLTZMANN = 2.042e-10  # MJ K^-4 m^-2 day^-1

# Atmospheric constants
ATMOSPHERIC_PRESSURE_SEA_LEVEL = 101.3  # kPa at sea level
GAS_CONSTANT_DRY_AIR = 0.287  # kJ kg^-1 K^-1
SPECIFIC_HEAT_AIR = 1.013  # kJ kg^-1 K^-1

# Thermodynamic constants
LATENT_HEAT_VAPORIZATION = 2.45  # MJ kg^-1
PSYCHROMETRIC_CONSTANT = 0.066  # kPa K^-1
WATER_DENSITY = 1000  # kg m^-3
SPECIFIC_HEAT_WATER = 4.18  # kJ kg^-1 K^-1

# Geometric constants
EARTH_RADIUS = 6371000  # meters (mean radius)
DEGREES_TO_RADIANS = math.pi / 180.0
RADIANS_TO_DEGREES = 180.0 / math.pi

# Time constants
DAYS_PER_YEAR = 365.25
HOURS_PER_DAY = 24
MINUTES_PER_HOUR = 60
SECONDS_PER_MINUTE = 60
SECONDS_PER_DAY = HOURS_PER_DAY * MINUTES_PER_HOUR * SECONDS_PER_MINUTE

# Precipitation analysis
EPSILON = 1e-10  # Small value for numerical stability in precipitation calculations
MIN_VALID_PRECIPITATION = 0.0  # Minimum valid precipitation value (mm)
MAX_VALID_PRECIPITATION = 1000.0  # Maximum reasonable daily precipitation (mm)

# Temperature thresholds
ABS_ZERO_CELSIUS = -273.15  # Absolute zero in Celsius
MIN_VALID_TEMPERATURE = -100.0  # Minimum reasonable temperature (°C)
MAX_VALID_TEMPERATURE = 60.0  # Maximum reasonable temperature (°C)
FREEZING_POINT = 0.0  # Celsius
BOILING_POINT = 100.0  # Celsius

# Gamma distribution parameters (for SPI)
GAMMA_SHAPE_MIN = 0.001  # Minimum shape parameter
GAMMA_SHAPE_MAX = 100.0  # Maximum shape parameter
GAMMA_SCALE_MIN = 0.001  # Minimum scale parameter
GAMMA_SCALE_MAX = 1000.0  # Maximum scale parameter

# Missing data indicator
MISSING_VALUE = -9999  # Standard missing data value
MISSING_VALUE_TOLERANCE = 1e-6  # Tolerance for float comparison with missing value

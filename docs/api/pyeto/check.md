# pyeto._check — Input Validation

Internal validation functions that enforce physical constraints on meteorological inputs. These functions raise `ValueError` when a value falls outside its valid range.

!!! note "Internal API"
    These functions are used internally by `pyeto.fao` functions. You only need to call them directly if you are building functions that accept the same inputs and want to enforce the same constraints.

---

::: pyeto._check.check_doy

::: pyeto._check.check_latitude_rad

::: pyeto._check.check_sol_dec_rad

::: pyeto._check.check_sunset_hour_angle_rad

::: pyeto._check.check_day_hours

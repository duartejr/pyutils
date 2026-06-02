# spi — Standardized Precipitation Index

Computes the SPI by fitting a Gamma distribution per season to accumulated precipitation, then transforming to standard normal variates.

Negative SPI indicates drought; positive indicates surplus.

## Reference

McKee, T. B., Doesken, N. J., & Kleist, J. (1993). *The relationship of drought frequency and duration to time scales.* Proceedings of the 8th Conference on Applied Climatology, 179-184.

---

## Main function

::: spi.spi
    options:
      show_signature_annotations: true

---

## Internal helpers

These helpers are documented for transparency. They are not part of the public API.

::: spi._accumulate_precipitation

::: spi._remove_zeros_and_nans

::: spi._compute_zero_probability

::: spi._gamma_to_standard_normal

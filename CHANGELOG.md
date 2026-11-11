# Changelog

## 0.1.0 (unreleased)

First version.

- `solve()`: TDOA on fine timestamps, RSSI multilateration with an unknown transmit power, and a centroid fallback, each with an accuracy radius from the covariance.
- **TDOA robustness:**
  - With 4 or more gateways, one gateway whose clock or reflection doesn't fit is left out.
  - When 3 gateways leave two equally good answers, the result is their midpoint, with a radius that reaches both.
- **Geometry checks:** gateways in a line give a position on the line, with a radius that reaches both sides. Gateways within 10 m of each other count as one place.
- **Input validation:** positions that are not finite, out of range, or on Null Island are ignored, as are RSSI values that are not numbers.
- **ChirpStack v4:**
  - `solve_events()`, `receptions()` and `gateway_positions()` read `rxInfo[]`.
  - `gps_time_ns()` parses `fineTimeSinceGpsEpoch` without losing nanoseconds.
- **Calibration:** keyword arguments for the path-loss exponent, shadowing, reference RSSI and timestamp noise.
- **`lorawan-locate` CLI:** ChirpStack events in, one JSON line per device out.
- **`[h3]` extra:** `lorawan_locate.coverage`, with H3 coverage cells, exact merging into parent cells, and GeoJSON output.
- Accuracy is measured on synthetic data only (`tests/test_accuracy.py`).

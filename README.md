# lorawan-locate

Estimates where a LoRaWAN device is from the gateways that heard it. It reads ChirpStack v4 uplink events directly, or plain gateway positions with RSSI and fine timestamps. It returns a position, an accuracy radius, and the method it used.

The core package uses the standard library only and needs no service or account. H3 coverage cells are an optional extra.

> **Accuracy has only been measured on synthetic data so far.** The [numbers below](#accuracy) come from simulated gateways and devices. They have not been checked against real receptions from devices at surveyed positions. If you have such data and can share results, please open an issue.

## Install

```bash
pip install git+https://github.com/pulzar-net/lorawan-locate          # solver and CLI
pip install 'lorawan-locate[h3] @ git+https://github.com/pulzar-net/lorawan-locate'   # plus coverage cells
```

It needs Python 3.9 or newer. The `[h3]` extra needs 3.10 or newer for current h3 wheels. It is not on PyPI yet.

## Quick start

Pass a list of ChirpStack `up` events from one device, as JSON from an MQTT or HTTP integration:

```python
import json
import lorawan_locate as ll

events = [json.loads(line) for line in open('uplinks.jsonl')]
print(ll.solve_events(events))
# {'latitude': 45.0110016, 'longitude': 10.0169989, 'accuracyM': 25.0, 'method': 'tdoa',
#  'gatewayIds': ['aa555a0000000001', 'aa555a0000000002', 'aa555a0000000003', 'aa555a0000000004'],
#  'uplinks': 1}
```

Or from the command line, with one JSON object, an array, or JSON lines, in a file or on stdin. It prints one line per `devEui`, with `location: null` when the device can't be placed:

```bash
lorawan-locate uplinks.jsonl
lorawan-locate --gateways gateways.json < uplinks.json
```

## Input

### ChirpStack v4 events

`solve_events(events, positions=None, **calibration)` uses each event's `time` and `rxInfo[]`. It reads these `rxInfo` fields:

| Field | Used for |
|---|---|
| `gatewayId` | Which gateway heard the uplink. Case is ignored. A gateway listed once per antenna counts as one reception: the loudest antenna's RSSI and SNR, and any antenna's fine timestamp. |
| `rssi` | RSSI multilateration and the centroid. A reception without it is skipped. |
| `fineTimeSinceGpsEpoch` | TDOA, as `"1400000000.000006043s"` or `{"seconds", "nanos"}`. It is parsed as text to the nanosecond, because a float would lose the digits that matter. |
| `location.latitude`, `location.longitude` | The gateway's position, used when `positions` is not given. Null Island (0, 0) is ignored. |

The location in `rxInfo` is whatever the gateway was configured with, and whoever runs that gateway chose it. If a gateway could lie to you, for example on a shared network, pass `positions` yourself. Then only the gateways in it are used, at the positions you give.

### Plain receptions

`solve(observations, positions, **calibration)` is the same solver without ChirpStack:

```python
ll.solve(
    [  # uplinks of one device, newest first
        [{'gateway_id': 'gw1', 'rssi': -97, 'fine_ns': None},
         {'gateway_id': 'gw2', 'rssi': -112, 'fine_ns': None},
         {'gateway_id': 'gw3', 'rssi': -104, 'fine_ns': None}],
    ],
    {'gw1': (45.000, 10.000), 'gw2': (45.000, 10.040), 'gw3': (45.025, 10.000)},
)
# {'latitude': 45.0090531, 'longitude': 10.0016279, 'accuracyM': 1567.9, 'method': 'rssi',
#  'gatewayIds': ['gw1', 'gw2', 'gw3'], 'uplinks': 1}
```

- **`positions`:** `{gateway_id: (latitude, longitude)}` in degrees. Gateways missing from it are ignored, as are entries that are not finite numbers, are out of range, or sit on Null Island.
- **`rssi`:** in dBm. A reception whose RSSI is missing or not a finite number is skipped.
- **`fine_ns`:** the arrival time in nanoseconds on a clock shared by every gateway: GPS time for ChirpStack. Only differences between gateways matter. `None` means no timestamp. `ll.gps_time_ns()` converts ChirpStack's text form.
- **`snr`:** optional and not used by the solver.

## Methods, and when each applies

`solve()` tries these in order and uses the first that works:

1. **TDOA (time difference of arrival).** Applies when 3 or more gateways at distinct places put a fine timestamp on the same uplink. That needs GPS-disciplined concentrators that produce fine timestamps, such as SX1303-based gateways. Each pair of arrival times puts the device on a hyperbola, and the solver finds where the hyperbolas cross. It uses Levenberg-Marquardt least squares from several starting points, inside and outside the gateways. TDOA uses the newest uplink that qualifies.
   - **Timestamp sanity:** a spread of arrival times wider than light could cover across the gateways means a clock is wrong. With 4 or more gateways, the solver retries leaving one out. It uses the smaller fit when exactly one of those fits works, or, with 5 or more gateways, the one that fits best. With 4, any 3 gateways fit exactly, so a merely poor fit can't say which gateway is wrong. It stands, with a radius widened by its residuals. With 3, a broken clock means no TDOA, and the solver falls back to RSSI.
   - **Two answers:** 3 gateways can leave two places that fit equally well. With the device outside the gateways, the hyperbolas often cross twice. When the solver finds both, it returns the midpoint, with a radius that reaches both.
2. **RSSI multilateration.** Applies with 3 or more gateways at distinct places that have RSSI. It fits a log-distance path-loss model, RSSI = A − 10·n·log10(d), solving for the position and A together. A is the RSSI at 1 m, which depends on the device's transmit power and antenna, so nobody knows it in advance. A weak prior on A makes 3 gateways enough. Each gateway's RSSI is averaged over every uplink given. That averages out fast fading, but not shadowing, which belongs to the place.
3. **Centroid.** Applies with 2 gateways at distinct places, or when the RSSI solve fails. It returns the midpoint weighted by modelled range, closer to the louder gateway. Its radius is at least half the distance between the gateways: it means "somewhere around these two".

Fewer than 2 distinct places returns `None`. Several gateways on one mast, within 10 m of each other, count as one place.

**Gateways in a line,** such as along a road, can't tell a device from its mirror image across the line. A TDOA or RSSI solve from them is moved onto the line, with a radius that reaches both sides. "In a line" means the gateways stray from it by less than 5% of its length.

### Output

| Key | |
|---|---|
| `latitude`, `longitude` | Degrees, 7 decimal places. |
| `accuracyM` | About a 1-sigma radial error in metres, never below 25. It is the square root of the trace of the position covariance, widened when residuals or the problems above call for it. The device is inside it most of the time, not always: see the next section. |
| `method` | `tdoa`, `rssi` or `centroid`. |
| `gatewayIds` | The gateways the solve used. A gateway left out for a bad clock is not listed. |
| `uplinks` | How many uplinks it used: 1 for TDOA, all of them for RSSI. |

## Accuracy

**Every number here comes from synthetic data.** To reproduce the table, run `python tests/test_accuracy.py`. The test suite also holds each row to bounds slightly looser than these, so a change that makes the solver worse fails CI.

Each row is 300 trials. The gateways sit on a 3 km square, and the device is placed at random inside the square unless the row says otherwise. The solver runs with its default calibration and is not told the real noise.

| Scenario | Median error | 90th percentile | Within accuracyM | Within 2x | Method used |
|---|---|---|---|---|---|
| TDOA, 4 gw, 10 ns | 3 m | 5 m | 100% | 100% | tdoa 300 |
| TDOA, 4 gw, 30 ns | 8 m | 15 m | 100% | 100% | tdoa 300 |
| TDOA, 4 gw, 100 ns | 26 m | 50 m | 58% | 94% | tdoa 300 |
| TDOA, 3 gw (triangle), 30 ns | 9 m | 18 m | 99% | 100% | tdoa 300 |
| TDOA, 4 gw, 30 ns, one link up to 500 ns late | 37 m | 76 m | 48% | 89% | tdoa 300 |
| TDOA, 5 gw, 30 ns, one link up to 500 ns late | 11 m | 42 m | 86% | 93% | tdoa 300 |
| TDOA, 5 gw, 30 ns, one clock 2 us off | 9 m | 17 m | 98% | 100% | tdoa 300 |
| TDOA, 3 gw along a road, 30 ns | 802 m | 1370 m | 92% | 100% | rssi 14, tdoa 286 |
| TDOA, 3 gw (triangle), 30 ns, outside | 96 m | 2172 m | 59% | 88% | rssi 13, tdoa 287 |
| TDOA, 4 gw, 30 ns, outside | 46 m | 127 m | 72% | 94% | tdoa 300 |
| RSSI, 5 gw, 4 dB | 521 m | 1147 m | 93% | 100% | rssi 300 |
| RSSI, 5 gw, 6 dB | 764 m | 1796 m | 85% | 99% | rssi 300 |
| RSSI, 5 gw, 8 dB | 1018 m | 2085 m | 83% | 99% | rssi 300 |
| RSSI, 3 gw, 6 dB | 974 m | 1941 m | 92% | 99% | rssi 300 |
| RSSI, 5 gw, 6 dB, n=2.5 actual | 706 m | 1318 m | 90% | 99% | rssi 300 |
| RSSI, 5 gw, 6 dB, n=3.5 actual | 1372 m | 3861 m | 82% | 100% | rssi 300 |
| RSSI, 5 gw, 6 dB, outside | 2332 m | 3872 m | 47% | 61% | rssi 300 |
| RSSI, 3 gw along a road, 6 dB | 870 m | 1733 m | 95% | 100% | rssi 300 |
| Centroid, 2 gw, 6 dB | 609 m | 999 m | 100% | 100% | centroid 300 |

### How the data was simulated

**TDOA rows:**
- Each timestamp gets Gaussian noise; that is the "ns" in each row.
- "One link up to 500 ns late": one random gateway per uplink also gets a reflection that arrives up to 500 ns late. This is a crude stand-in for multipath in a built-up area.
- "One clock 2 us off": one gateway's clock is 2 µs out.
- "Outside": the device is 1–3 km beyond the gateways.
- "Along a road": 3 gateways in a 3 km line with a 60 m bend, and the device up to 1.5 km to either side.

**RSSI rows:**
- Shadowing (the "dB") is fixed for each device–gateway link, plus 3 dB of fast fading on each of 10 uplinks.
- The transmit power varies by ±3 dB around −12 dBm at 1 m.
- In "n=… actual", the real path-loss exponent differs from the 3.0 the solver assumes.

### What the table shows

- **TDOA is good to tens of metres with clean timestamps and the device inside the gateways.** Timestamp quality sets the error.
- **With 4 gateways, one late reflection roughly quadruples the error,** and the radius is then too small about half the time. With 5 or more, the solver can drop the bad gateway, and the error mostly recovers.
- **TDOA with 3 gateways and the device outside is often ambiguous.** The 90th percentile is about 2 km. The radius covers that ambiguity 88% of the time at 2×, not always.
- **RSSI positions are good to about a kilometre** with gateways 3 km apart and realistic shadowing. That is enough to find a device's area, not its exact spot.
- **Outside the gateways, the RSSI radius is overconfident.** The device is within 2× the radius only 61% of the time.
- **Not simulated:** errors in gateway positions, antenna patterns, fixed delays that differ between gateway models, altitude, and correlated multipath. Real data will be worse than this table in ways it can't show.

## Calibration

`solve()` and `solve_events()` take a site's calibration as keyword arguments:

| Argument | Default | |
|---|---|---|
| `path_loss_exponent` | 3.0 | n in the path-loss model; about 2 in open country, 3–4 in towns. |
| `shadowing_db` | 6.0 | Standard deviation of shadowing. |
| `reference_rssi` | −15.0 | Expected RSSI at 1 m, in dBm. It is solved for anyway; this is the prior. |
| `timestamp_sigma_ns` | 30.0 | What a fine timestamp is good to. |

Fit them from devices whose position you know.

A gateway's cable and RF-chain delay (about 5 ns per metre of coax) shifts every one of its fine timestamps by the same amount. To correct for it, subtract the delay from that gateway's `fine_ns` before solving.

## Coverage cells (`[h3]` extra)

`lorawan_locate.coverage` bins signal samples into H3 resolution-9 hexagons, about 200 m across. It merges them exactly into parent cells for zoomed-out views: counts and sums add, maxima take the maximum, and the time span widens. The cells come out as GeoJSON.

```python
from lorawan_locate import coverage

cell = coverage.cell_of(45.0123, 10.0456)
rows = [{'gateway_id': 'gw1', 'cell': cell, 'samples': 3, 'rssi_sum': -312, 'rssi_max': -101}]
geojson = coverage.collection(coverage.merge(rows, 7))   # best gateway per cell, at resolution 7
```

Each row needs `gateway_id`, `cell`, `samples`, `rssi_sum` and `rssi_max`. It can also carry `snr_sum`, `snr_samples`, `snr_max`, `first_at` and `last_at`. `merge()` carries any other keys along.

`coverage.plausible(lat, lon, gateway_latlon)` rejects a position that is Null Island, out of range, or more than 100 km from the gateway that heard it. `coverage.grade(rssi, snr)` gives `good`, `fair` or `weak`.

## Model limits

- The Earth is treated as flat around the gateways, with an equirectangular projection. Over tens of km this is a few metres off.
- Altitude is ignored.
- TDOA range differences share a reference gateway but are treated as independent.
- Each TDOA fix uses one uplink.
- The solver is pure Python. A solve takes milliseconds, which is fine for a device every few minutes and slow for millions an hour.

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md). Real-world accuracy results are the most useful contribution.

## License

Apache-2.0. See [LICENSE](LICENSE).

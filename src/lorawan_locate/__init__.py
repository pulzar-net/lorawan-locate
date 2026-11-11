"""Where a LoRaWAN device probably is, from the gateways that heard it.

Pure geometry, standard library only: no database, no network. The caller
decides which gateways to believe and hands their positions in.

Three methods, tried in this order:

  tdoa      Time difference of arrival. Three or more gateways that stamped
            the same uplink with a fine timestamp (fineTimeSinceGpsEpoch --
            GPS-disciplined concentrators only) place it by the differences
            between those times: each pair puts the device on a hyperbola,
            and the hyperbolas cross at the device. Solved by Gauss-Newton on
            the range differences against the earliest gateway, from several
            starting points so a local minimum does not win.
  rssi      Three or more gateways with RSSI. A log-distance path-loss model,
            RSSI = A - 10 n log10(d), solved for the position and A together
            (A, the RSSI a metre away, depends on the device's power and
            antenna, which nobody told us), with a weak prior on A so three
            gateways are enough. Shadowing makes RSSI ranging coarse, and the
            accuracy radius says so.
  centroid  Two gateways, or a failed RSSI solve. The RSSI-weighted midpoint,
            with a radius as wide as the honest answer: somewhere around
            these two.

One gateway places nothing.

Simplifications, said once here rather than discovered: a flat Earth around
the gateways (an equirectangular projection, metres off over tens of km),
gateway altitude ignored, and the range differences treated as independent
although they share the reference gateway.

The accuracy radius is the root of the trace of the position covariance:
about a 1-sigma radial error, inside which the device is most of the time and
not always. tests/test_accuracy.py measures how often.
"""

import math

__version__ = '0.1.0'

EARTH_RADIUS_M = 6_371_008.8
SPEED_OF_LIGHT = 299_792_458.0

# Urban to suburban at 915 and 868 MHz. Not fitted per site: a fit needs
# devices of known position. Pass your own to solve() once you have them.
PATH_LOSS_EXPONENT = 3.0
SHADOWING_DB = 6.0
# RSSI one metre away for a 14 dBm end device with unity-gain antennas, and
# how far a real one strays from it.
REFERENCE_RSSI = -15.0
REFERENCE_SIGMA_DB = 15.0
# What a fine timestamp is good to. Semtech quotes tens of nanoseconds for
# the SX1303's; 30 ns is about 9 m of range.
TIMESTAMP_SIGMA_NS = 30.0

MIN_ACCURACY_M = 25.0
# A TDOA fit whose residuals average more than this many sigma squared has a
# timestamp in it that does not belong: a clock off, or a reflection.
OUTLIER_CHI2 = 9.0
# Gateways that stray from a straight line by less than this share of its
# length are on it, for the mirror ambiguity that brings.
# ponytail: a fixed share; a covariance test if bent lines ever matter.
COLLINEAR = 0.05
# Gateways closer together than this are one place.
_SAME_PLACE_M = 10.0
# A solve further than this from every gateway it used is a broken clock or
# a lie, not a device.
MAX_RANGE_M = 50_000.0
_MIN_RANGE_M = 10.0


class Projection:
    """Metres east and north of a reference point."""

    def __init__(self, lat0, lon0):
        self.lat0, self.lon0 = lat0, lon0
        self.k = math.cos(math.radians(lat0))

    def xy(self, lat, lon):
        return (math.radians(lon - self.lon0) * EARTH_RADIUS_M * self.k,
                math.radians(lat - self.lat0) * EARTH_RADIUS_M)

    def latlon(self, x, y):
        return (self.lat0 + math.degrees(y / EARTH_RADIUS_M),
                self.lon0 + math.degrees(x / (EARTH_RADIUS_M * self.k)))


def distance_m(lat1, lon1, lat2, lon2):
    """Great-circle distance (haversine)."""
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dp, dl = p2 - p1, math.radians(lon2 - lon1)
    a = math.sin(dp / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2
    return 2 * EARTH_RADIUS_M * math.asin(min(1.0, math.sqrt(a)))


# ── small dense linear algebra ─────────────────────────────────────────────


def _solve(matrix, vector):
    """Gaussian elimination with partial pivoting. None when singular."""
    n = len(vector)
    a = [list(row) + [vector[i]] for i, row in enumerate(matrix)]
    for col in range(n):
        pivot = max(range(col, n), key=lambda r: abs(a[r][col]))
        if abs(a[pivot][col]) < 1e-12:
            return None
        a[col], a[pivot] = a[pivot], a[col]
        for r in range(n):
            if r != col:
                f = a[r][col] / a[col][col]
                for c in range(col, n + 1):
                    a[r][c] -= f * a[col][c]
    return [a[i][n] / a[i][i] for i in range(n)]


def _inverse(matrix):
    n = len(matrix)
    columns = []
    for i in range(n):
        col = _solve(matrix, [1.0 if j == i else 0.0 for j in range(n)])
        if col is None:
            return None
        columns.append(col)
    return [[columns[j][i] for j in range(n)] for i in range(n)]


def _normal(jacobian, residuals):
    """JᵀJ and Jᵀr."""
    k = len(jacobian[0])
    jtj = [[sum(row[i] * row[j] for row in jacobian) for j in range(k)] for i in range(k)]
    jtr = [sum(row[i] * r for row, r in zip(jacobian, residuals)) for i in range(k)]
    return jtj, jtr


def _least_squares(model, start, iterations=60):
    """Levenberg-Marquardt. `model(theta)` -> (residuals, jacobian), both
    already divided by each measurement's sigma. Returns (theta, cost, jtj)."""
    theta = list(start)
    residuals, jacobian = model(theta)
    cost = sum(r * r for r in residuals)
    damping = 1e-3
    for _ in range(iterations):
        jtj, jtr = _normal(jacobian, residuals)
        damped = [[v + (damping * jtj[i][i] if i == j else 0.0) for j, v in enumerate(row)]
                  for i, row in enumerate(jtj)]
        step = _solve(damped, [-g for g in jtr])
        if step is None:
            break
        trial = [t + s for t, s in zip(theta, step)]
        trial_residuals, trial_jacobian = model(trial)
        trial_cost = sum(r * r for r in trial_residuals)
        if trial_cost < cost:
            theta, residuals, jacobian = trial, trial_residuals, trial_jacobian
            improvement, cost = cost - trial_cost, trial_cost
            damping = max(damping / 3, 1e-9)
            if improvement < 1e-9 * (1 + cost) and max(abs(s) for s in step) < 1e-3:
                break
        else:
            damping *= 4
            if damping > 1e9:
                break
    jtj, _ = _normal(jacobian, residuals)
    return theta, cost, jtj


def _radius(covariance):
    if covariance is None:
        return None
    trace = covariance[0][0] + covariance[1][1]
    return math.sqrt(trace) if trace > 0 else None


def _starts(points, weights=None):
    """Where to start a solve: the (weighted) centroid, and a point just
    inside each gateway, so a minimum near any of them can be found."""
    weights = weights or [1.0] * len(points)
    total = sum(weights)
    cx = sum(p[0] * w for p, w in zip(points, weights)) / total
    cy = sum(p[1] * w for p, w in zip(points, weights)) / total
    starts = [(cx, cy)]
    for px, py in points:
        starts.append((px + 0.2 * (cx - px) + 1.0, py + 0.2 * (cy - py) + 1.0))
    return starts


# ── the methods, in metres on a plane ──────────────────────────────────────


def tdoa(points, times_ns, sigma_ns=TIMESTAMP_SIGMA_NS):
    """(x, y, accuracy_m) from >= 3 points and their arrival times, or None."""
    solved = _tdoa(points, times_ns, sigma_ns)
    return solved and solved[:3]


def _tdoa(points, times_ns, sigma_ns):
    """tdoa() and how well it fit: (x, y, accuracy_m, cost per degree of
    freedom, 0 when there is none to spare)."""
    if len(points) < 3:
        return None
    ref = min(range(len(points)), key=lambda i: times_ns[i])
    # Differences beyond what light covers across the furthest pair are clocks
    # that do not agree, not geometry.
    spread = max(math.dist(a, b) for a in points for b in points)
    if (max(times_ns) - min(times_ns)) * 1e-9 * SPEED_OF_LIGHT > spread + 1000:
        return None
    others = [i for i in range(len(points)) if i != ref]
    sigma_m = sigma_ns * 1e-9 * SPEED_OF_LIGHT
    rx, ry = points[ref]

    def model(theta):
        x, y = theta
        d_ref = max(math.hypot(x - rx, y - ry), 1e-6)
        residuals, jacobian = [], []
        for i in others:
            gx, gy = points[i]
            d_i = max(math.hypot(x - gx, y - gy), 1e-6)
            measured = (times_ns[i] - times_ns[ref]) * 1e-9 * SPEED_OF_LIGHT
            residuals.append(((d_i - d_ref) - measured) / sigma_m)
            jacobian.append([((x - gx) / d_i - (x - rx) / d_ref) / sigma_m,
                             ((y - gy) / d_i - (y - ry) / d_ref) / sigma_m])
        return residuals, jacobian

    # Starts outside the gateways too: hyperbolas can cross twice, and the
    # second crossing is usually out there.
    cx = sum(p[0] for p in points) / len(points)
    cy = sum(p[1] for p in points) / len(points)
    starts = _starts(points) + [(cx + k * (px - cx) + 1.0, cy + k * (py - cy) + 1.0)
                                for px, py in points for k in (2, 4)]
    fits = [_least_squares(model, start) for start in starts]
    theta, cost, jtj = min(fits, key=lambda f: f[1])
    # Residuals larger than the clocks allow say the timestamps disagree; the
    # residual's own size then sets the error, not the optimistic sigma.
    dof = len(others) - 2
    scale = max(1.0, cost / dof) if dof > 0 else 1.0
    covariance = _inverse(jtj)
    if covariance is None:
        return None
    radius = _radius([[v * scale for v in row] for row in covariance])
    if radius is None:
        return None
    x, y = theta
    # Another place that fits as well (within one sigma squared) and is not
    # this one: the timestamps cannot choose, so neither does the answer. It
    # is the midpoint, with a radius that reaches both.
    rivals = [f[0] for f in fits if f[1] <= cost + 1.0 and math.dist(f[0], theta) > 3 * radius]
    if rivals:
        other = max(rivals, key=lambda t: math.dist(t, theta))
        half = math.dist(other, theta) / 2
        x, y, radius = (x + other[0]) / 2, (y + other[1]) / 2, half + radius
    return x, y, radius, cost / dof if dof > 0 else 0.0


def rssi_multilateration(points, rssi, n=PATH_LOSS_EXPONENT, shadowing=SHADOWING_DB,
                         reference=REFERENCE_RSSI):
    """(x, y, accuracy_m) from >= 3 points and their RSSI, or None."""
    if len(points) < 3:
        return None
    factor = 10 * n / math.log(10)

    def model(theta):
        x, y, a = theta
        residuals, jacobian = [], []
        for (gx, gy), measured in zip(points, rssi):
            d2 = max((x - gx) ** 2 + (y - gy) ** 2, _MIN_RANGE_M ** 2)
            predicted = a - 10 * n * math.log10(math.sqrt(d2))
            residuals.append((predicted - measured) / shadowing)
            jacobian.append([-factor * (x - gx) / d2 / shadowing,
                             -factor * (y - gy) / d2 / shadowing,
                             1.0 / shadowing])
        residuals.append((a - reference) / REFERENCE_SIGMA_DB)
        jacobian.append([0.0, 0.0, 1.0 / REFERENCE_SIGMA_DB])
        return residuals, jacobian

    weights = [10 ** (r / (10 * n)) for r in rssi]
    best = None
    for sx, sy in _starts(points, weights):
        theta, cost, jtj = _least_squares(model, (sx, sy, reference))
        if best is None or cost < best[1]:
            best = (theta, cost, jtj)
    theta, cost, jtj = best
    covariance = _inverse(jtj)
    if covariance is None:
        return None
    dof = len(points) - 3
    scale = max(1.0, cost / dof) if dof > 0 else 1.0
    radius = _radius([[v * scale for v in row[:2]] for row in covariance[:2]])
    if radius is None:
        return None
    return theta[0], theta[1], radius


def centroid(points, rssi, n=PATH_LOSS_EXPONENT, reference=REFERENCE_RSSI):
    """(x, y, accuracy_m): the RSSI-weighted centroid. Weighted by the inverse
    of each gateway's modelled range, so a gateway 10 dB louder pulls about
    twice as hard at n = 3."""
    if not points:
        return None
    weights = [10 ** ((r - reference) / (10 * n)) for r in rssi]
    total = sum(weights)
    x = sum(p[0] * w for p, w in zip(points, weights)) / total
    y = sum(p[1] * w for p, w in zip(points, weights)) / total
    # Two gateways say "somewhere around here", not "on the line between".
    reach = max(math.hypot(p[0] - x, p[1] - y) for p in points)
    span = max(math.dist(a, b) for a in points for b in points)
    return x, y, max(reach, span / 2)


# ── from receptions ─────────────────────────────────────────────────────────


def solve(observations, positions, *, path_loss_exponent=PATH_LOSS_EXPONENT,
          shadowing_db=SHADOWING_DB, reference_rssi=REFERENCE_RSSI,
          timestamp_sigma_ns=TIMESTAMP_SIGMA_NS):
    """A position from receptions, or None.

    `observations`: uplinks of one device, newest first, each a list of
    receptions {gateway_id, rssi, snr, fine_ns}. `positions`: {gateway_id:
    (lat, lon)} for the gateways that may be used; receptions from any other
    gateway are ignored. Leaving a gateway out of `positions` is how a caller
    says it does not believe it.

    TDOA uses the newest uplink heard by three or more timestamping gateways.
    RSSI averages each gateway over every uplink given, in dB: that averages
    out fast fading, not shadowing, which belongs to the place.

    The keyword arguments are the site's calibration. Subtract a gateway's
    known cable and RF-chain delay from its fine_ns before calling.
    """
    positions = {g: p for g, p in positions.items() if _valid_position(p)}
    usable = [[r for r in uplink if r.get('gateway_id') in positions and _finite(r.get('rssi'))]
              for uplink in observations]
    heard = {r['gateway_id'] for uplink in usable for r in uplink}
    if len(heard) < 2:
        return None
    lat0 = sum(positions[g][0] for g in heard) / len(heard)
    lon0 = sum(positions[g][1] for g in heard) / len(heard)
    proj = Projection(lat0, lon0)
    xy = {g: proj.xy(*positions[g]) for g in heard}
    if _places([xy[g] for g in heard]) < 2:
        return None

    result = None
    for uplink in usable:
        timed = {}
        for r in uplink:
            if _finite(r.get('fine_ns')):
                timed.setdefault(r['gateway_id'], r['fine_ns'])
        if len(timed) >= 3:
            solved, ids = _tdoa_robust(sorted(timed), xy, timed, timestamp_sigma_ns)
            if solved:
                result = (*solved[:3], 'tdoa', ids, 1)
                break

    if result is None:
        sums = {}
        for uplink in usable:
            for r in uplink:
                total, count = sums.get(r['gateway_id'], (0.0, 0))
                sums[r['gateway_id']] = (total + r['rssi'], count + 1)
        ids = sorted(sums)
        mean = [sums[g][0] / sums[g][1] for g in ids]
        points = [xy[g] for g in ids]
        solved = rssi_multilateration(points, mean, path_loss_exponent, shadowing_db,
                                      reference_rssi) if _places(points) >= 3 else None
        method = 'rssi'
        if solved is None or not _plausible(solved, points):
            solved, method = centroid(points, mean, path_loss_exponent, reference_rssi), 'centroid'
        result = (*solved, method, ids, sum(1 for u in usable if u))

    x, y, accuracy, method, ids, uplinks = result
    if method != 'centroid':
        x, y, accuracy = _unfold((x, y, accuracy), [xy[g] for g in ids])
    lat, lon = proj.latlon(x, y)
    return {
        'latitude': round(lat, 7),
        'longitude': round(lon, 7),
        'accuracyM': round(max(accuracy, MIN_ACCURACY_M), 1),
        'method': method,
        'gatewayIds': ids,
        'uplinks': uplinks,
    }


def _tdoa_robust(ids, xy, timed, sigma_ns):
    """(TDOA solve, the gateways it used) or (None, ids). When the fit over
    every gateway is impossible or poor, leaving out one gateway is tried:
    a broken clock or a reflection is usually one gateway's. The smaller fit
    is used when only one such fit works, or, with five or more gateways, the
    one that fits best -- with four, any three fit exactly, so a poor fit
    cannot say which gateway is wrong and stands, with its widened radius."""
    def fit(subset):
        points = [xy[g] for g in subset]
        if _places(points) < 3:
            return None
        solved = _tdoa(points, [timed[g] for g in subset], sigma_ns)
        return solved if solved and _plausible(solved[:3], points) else None

    full = fit(ids)
    if len(ids) < 4 or (full is not None and full[3] <= OUTLIER_CHI2):
        return full, ids
    fits = []
    for drop in ids:
        subset = [g for g in ids if g != drop]
        solved = fit(subset)
        if solved is not None and solved[3] <= OUTLIER_CHI2:
            fits.append((solved, subset))
    if len(fits) == 1 or (fits and len(ids) >= 5):
        return min(fits, key=lambda f: f[0][3])
    return full, ids


def _unfold(solved, points):
    """Gateways in a line cannot tell a device from its mirror image across
    the line, so a solve from them is moved onto the line with a radius that
    reaches both sides."""
    x, y, accuracy = solved
    n = len(points)
    cx, cy = sum(p[0] for p in points) / n, sum(p[1] for p in points) / n
    sxx = sum((p[0] - cx) ** 2 for p in points)
    syy = sum((p[1] - cy) ** 2 for p in points)
    sxy = sum((p[0] - cx) * (p[1] - cy) for p in points)
    angle = 0.5 * math.atan2(2 * sxy, sxx - syy)
    ux, uy = math.cos(angle), math.sin(angle)
    along = [(p[0] - cx) * ux + (p[1] - cy) * uy for p in points]
    across = [(p[1] - cy) * ux - (p[0] - cx) * uy for p in points]
    if max(map(abs, across)) > COLLINEAR * (max(along) - min(along)):
        return solved
    offset = (y - cy) * ux - (x - cx) * uy
    return x + offset * uy, y - offset * ux, accuracy + abs(offset)


def _places(points):
    """How many distinct places these gateways are at."""
    # ponytail: a grid, so two gateways either side of a line 1 m apart count
    # as two; fine for "are these all one mast".
    return len({(round(x / _SAME_PLACE_M), round(y / _SAME_PLACE_M)) for x, y in points})


def _finite(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def _valid_position(position):
    try:
        lat, lon = position
    except (TypeError, ValueError):
        return False
    return (_finite(lat) and _finite(lon) and -90 <= lat <= 90 and -180 <= lon <= 180
            and (abs(lat) > 1e-6 or abs(lon) > 1e-6))


def _plausible(solved, points):
    x, y, accuracy = solved
    if not all(math.isfinite(v) for v in (x, y, accuracy)):
        return False
    return min(math.hypot(x - px, y - py) for px, py in points) <= MAX_RANGE_M


# ── ChirpStack v4 ──────────────────────────────────────────────────────────


def gps_time_ns(value):
    """A protobuf Duration as ChirpStack renders it in JSON --
    "1380901234.000341892s", or {"seconds", "nanos"} -- in integer
    nanoseconds, or None. Parsed from the text, not through float: a float
    holds about 16 digits and this has 19, so the nanoseconds that make a fine
    timestamp worth having would be rounded away."""
    if isinstance(value, dict):
        try:
            return int(value.get('seconds') or 0) * 1_000_000_000 + int(value.get('nanos') or 0)
        except (TypeError, ValueError):
            return None
    if not isinstance(value, str) or not value.endswith('s'):
        return None
    whole, _, fraction = value[:-1].partition('.')
    if not whole.isdigit() or (fraction and not fraction.isdigit()) or len(fraction) > 9:
        return None
    return int(whole) * 1_000_000_000 + int((fraction or '0').ljust(9, '0'))


def receptions(rx_info):
    """ChirpStack `rxInfo[]` as receptions for solve(): one per gateway,
    lower-cased, keeping the loudest antenna's RSSI and SNR -- two antennas of
    one gateway are not two places the signal was heard -- and any antenna's
    fine timestamp."""
    best = {}
    for rx in rx_info or []:
        gateway = (rx.get('gatewayId') or '').lower()
        if not gateway or rx.get('rssi') is None:
            continue
        fine = gps_time_ns(rx.get('fineTimeSinceGpsEpoch'))
        kept = best.get(gateway)
        if kept is not None and kept['rssi'] >= rx['rssi']:
            if kept['fine_ns'] is None:
                kept['fine_ns'] = fine
            continue
        best[gateway] = {'gateway_id': gateway, 'rssi': rx['rssi'], 'snr': rx.get('snr'),
                         'fine_ns': fine if fine is not None else kept and kept['fine_ns']}
    return list(best.values())


def gateway_positions(events):
    """{gateway_id: (lat, lon)} from the `location` each rxInfo carries. That
    is the position the gateway was configured with, which whoever runs it
    chose: pass your own to solve_events() for gateways you do not trust."""
    positions = {}
    for event in events:
        for rx in event.get('rxInfo') or []:
            loc = rx.get('location') or {}
            lat, lon = loc.get('latitude'), loc.get('longitude')
            gateway = (rx.get('gatewayId') or '').lower()
            if gateway and lat is not None and lon is not None and (abs(lat) > 1e-6 or abs(lon) > 1e-6):
                positions[gateway] = (lat, lon)
    return positions


def solve_events(events, positions=None, **calibration):
    """solve() from ChirpStack v4 uplink events (the JSON of an `up` event)
    of one device. Sorted newest first by `time`; positions default to the
    ones the events carry."""
    # ponytail: sorts `time` as text, right for one LNS's format; parse it if
    # events from mixed sources ever disagree on offsets or digits.
    events = sorted(events, key=lambda e: e.get('time') or '', reverse=True)
    if positions is None:
        positions = gateway_positions(events)
    return solve([receptions(e.get('rxInfo')) for e in events], positions, **calibration)

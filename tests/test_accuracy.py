"""How far off the solver is, measured on synthetic geometry.

A device at a known point, gateways at known points, and what each gateway
would report generated from them: arrival times with Gaussian clock noise;
RSSI from the log-distance model with shadowing that belongs to the link (the
same on every uplink, so averaging cannot remove it) plus fast fading that
does not. Then solve() has to find the device again, with its default
calibration, which is not told the noise.

`python tests/test_accuracy.py` prints the table the README quotes. The test
holds each row to bounds a little looser than what was measured, so a change
that makes the solver worse, or its accuracy radius less honest, fails.
"""

import math
import random

import lorawan_locate as ll

C = ll.SPEED_OF_LIGHT
SQUARE = [(0, 0), (3000, 0), (0, 3000), (3000, 3000)]
FIVE = SQUARE + [(1500, -1000)]
PAIR = [(0, 0), (3000, 0)]
ROAD = [(0, 0), (1500, 60), (3000, 0)]
TRIALS = 300


def inside(rng):
    return rng.uniform(300, 2700), rng.uniform(300, 2700)


def outside(rng):
    """1 to 3 km east of the square: the geometry TDOA is worst at."""
    return rng.uniform(4000, 6000), rng.uniform(0, 3000)


def between(rng):
    return rng.uniform(0, 3000), rng.uniform(-1000, 1000)


def tdoa_uplinks(rng, tx, gateways, sigma_ns):
    t0 = 1_400_000_000_000_000_000
    return [[{'gateway_id': f'g{i}', 'rssi': -100,
              'fine_ns': t0 + round(math.dist(tx, g) / C * 1e9 + rng.gauss(0, sigma_ns))}
             for i, g in enumerate(gateways)]]


def nlos_uplinks(rng, tx, gateways, sigma_ns, late_ns=500):
    """One gateway, at random, hears a reflection: up to `late_ns` late, as
    multipath in a built-up area does. Never early -- nothing beats the
    straight line."""
    uplinks = tdoa_uplinks(rng, tx, gateways, sigma_ns)
    rng.choice(uplinks[0])['fine_ns'] += round(rng.uniform(0, late_ns))
    return uplinks


def bad_clock_uplinks(rng, tx, gateways, sigma_ns, off_ns=2000):
    """One gateway, at random, with its clock `off_ns` out."""
    uplinks = tdoa_uplinks(rng, tx, gateways, sigma_ns)
    rng.choice(uplinks[0])['fine_ns'] += off_ns
    return uplinks


def rssi_uplinks(rng, tx, gateways, shadowing, n=3.0, fading=3.0, uplinks=10):
    power = -12 + rng.gauss(0, 3)
    links = [rng.gauss(0, shadowing) for _ in gateways]
    return [[{'gateway_id': f'g{i}',
              'rssi': power - 10 * n * math.log10(max(math.dist(tx, g), 10)) + s + rng.gauss(0, fading),
              'fine_ns': None}
             for i, (g, s) in enumerate(zip(gateways, links))]
            for _ in range(uplinks)]


SCENARIOS = [
    # (name, gateways, where, uplinks(rng, tx, gateways))
    ('TDOA, 4 gw, 10 ns', SQUARE, inside, lambda r, t, g: tdoa_uplinks(r, t, g, 10)),
    ('TDOA, 4 gw, 30 ns', SQUARE, inside, lambda r, t, g: tdoa_uplinks(r, t, g, 30)),
    ('TDOA, 4 gw, 100 ns', SQUARE, inside, lambda r, t, g: tdoa_uplinks(r, t, g, 100)),
    ('TDOA, 3 gw (triangle), 30 ns', SQUARE[:3], lambda r: (r.uniform(300, 1800), r.uniform(300, 1800)),
     lambda r, t, g: tdoa_uplinks(r, t, g, 30)),
    ('TDOA, 4 gw, 30 ns, one link up to 500 ns late', SQUARE, inside,
     lambda r, t, g: nlos_uplinks(r, t, g, 30)),
    ('TDOA, 5 gw, 30 ns, one link up to 500 ns late', FIVE, inside,
     lambda r, t, g: nlos_uplinks(r, t, g, 30)),
    ('TDOA, 5 gw, 30 ns, one clock 2 us off', FIVE, inside,
     lambda r, t, g: bad_clock_uplinks(r, t, g, 30)),
    ('TDOA, 3 gw along a road, 30 ns', ROAD, lambda r: (r.uniform(0, 3000), r.uniform(-1500, 1500)),
     lambda r, t, g: tdoa_uplinks(r, t, g, 30)),
    ('TDOA, 3 gw (triangle), 30 ns, outside', SQUARE[:3], outside,
     lambda r, t, g: tdoa_uplinks(r, t, g, 30)),
    ('TDOA, 4 gw, 30 ns, outside', SQUARE, outside, lambda r, t, g: tdoa_uplinks(r, t, g, 30)),
    ('RSSI, 5 gw, 4 dB', FIVE, inside, lambda r, t, g: rssi_uplinks(r, t, g, 4)),
    ('RSSI, 5 gw, 6 dB', FIVE, inside, lambda r, t, g: rssi_uplinks(r, t, g, 6)),
    ('RSSI, 5 gw, 8 dB', FIVE, inside, lambda r, t, g: rssi_uplinks(r, t, g, 8)),
    ('RSSI, 3 gw, 6 dB', SQUARE[:3], inside, lambda r, t, g: rssi_uplinks(r, t, g, 6)),
    ('RSSI, 5 gw, 6 dB, n=2.5 actual', FIVE, inside, lambda r, t, g: rssi_uplinks(r, t, g, 6, n=2.5)),
    ('RSSI, 5 gw, 6 dB, n=3.5 actual', FIVE, inside, lambda r, t, g: rssi_uplinks(r, t, g, 6, n=3.5)),
    ('RSSI, 5 gw, 6 dB, outside', FIVE, outside, lambda r, t, g: rssi_uplinks(r, t, g, 6)),
    ('RSSI, 3 gw along a road, 6 dB', ROAD, lambda r: (r.uniform(0, 3000), r.uniform(-1500, 1500)),
     lambda r, t, g: rssi_uplinks(r, t, g, 6)),
    ('Centroid, 2 gw, 6 dB', PAIR, between, lambda r, t, g: rssi_uplinks(r, t, g, 6)),
]

# name: (median m, 90th percentile m, share within accuracyM, share within 2x)
# -- the bounds the test holds, set a little looser than measured.
BOUNDS = {
    'TDOA, 4 gw, 10 ns': (5, 10, 0.95, 0.97),
    'TDOA, 4 gw, 30 ns': (10, 20, 0.95, 0.97),
    'TDOA, 4 gw, 100 ns': (35, 65, 0.52, 0.9),
    'TDOA, 3 gw (triangle), 30 ns': (15, 25, 0.93, 0.97),
    'TDOA, 4 gw, 30 ns, one link up to 500 ns late': (50, 100, 0.43, 0.86),
    'TDOA, 5 gw, 30 ns, one link up to 500 ns late': (15, 55, 0.81, 0.9),
    'TDOA, 5 gw, 30 ns, one clock 2 us off': (15, 25, 0.93, 0.97),
    'TDOA, 3 gw along a road, 30 ns': (1050, 1750, 0.87, 0.97),
    'TDOA, 3 gw (triangle), 30 ns, outside': (120, 2650, 0.54, 0.85),
    'TDOA, 4 gw, 30 ns, outside': (60, 200, 0.67, 0.91),
    'RSSI, 5 gw, 4 dB': (700, 1450, 0.88, 0.97),
    'RSSI, 5 gw, 6 dB': (1000, 2250, 0.8, 0.96),
    'RSSI, 5 gw, 8 dB': (1300, 2650, 0.77, 0.96),
    'RSSI, 3 gw, 6 dB': (1250, 2450, 0.87, 0.96),
    'RSSI, 5 gw, 6 dB, n=2.5 actual': (900, 1650, 0.85, 0.96),
    'RSSI, 5 gw, 6 dB, n=3.5 actual': (1750, 4850, 0.76, 0.97),
    # Outside the gateways the RSSI radius is overconfident: the README says so.
    'RSSI, 5 gw, 6 dB, outside': (2950, 4850, 0.42, 0.57),
    'RSSI, 3 gw along a road, 6 dB': (1100, 2200, 0.89, 0.97),
    'Centroid, 2 gw, 6 dB': (800, 1250, 0.94, 0.97),
}


def measure(name, gateways, where, uplinks, trials=TRIALS, seed=11):
    rng = random.Random(seed)
    proj = ll.Projection(45.0, 10.0)
    positions = {f'g{i}': proj.latlon(*g) for i, g in enumerate(gateways)}
    errors, within, within2, methods = [], 0, 0, {}
    for _ in range(trials):
        tx = where(rng)
        solved = ll.solve(uplinks(rng, tx, gateways), positions)
        error = ll.distance_m(solved['latitude'], solved['longitude'], *proj.latlon(*tx))
        errors.append(error)
        within += error <= solved['accuracyM']
        within2 += error <= 2 * solved['accuracyM']
        methods[solved['method']] = methods.get(solved['method'], 0) + 1
    errors.sort()
    return {'name': name, 'median': errors[len(errors) // 2], 'p90': errors[int(0.9 * len(errors))],
            'within': within / trials, 'within2': within2 / trials,
            'methods': methods}


def test_accuracy():
    for scenario in SCENARIOS:
        row = measure(*scenario)
        median, p90, within, within2 = BOUNDS[row['name']]
        assert row['median'] <= median, row
        assert row['p90'] <= p90, row
        assert row['within'] >= within, row
        assert row['within2'] >= within2, row


def table():
    lines = ['| Scenario | Median error | 90th percentile | Within accuracyM | Within 2x | Method used |',
             '|---|---|---|---|---|---|']
    for scenario in SCENARIOS:
        r = measure(*scenario)
        methods = ', '.join(f'{m} {c}' for m, c in sorted(r['methods'].items()))
        lines.append(f"| {r['name']} | {r['median']:.0f} m | {r['p90']:.0f} m | "
                     f"{r['within']:.0%} | {r['within2']:.0%} | {methods} |")
    return '\n'.join(lines)


if __name__ == '__main__':
    print(table())

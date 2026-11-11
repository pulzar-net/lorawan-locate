"""The solver's behaviour, and reading ChirpStack events.

Pinned:
  * through solve(): RSSI averages uplinks, TDOA wins when the newest uplink
    has three timestamping gateways, timestamps that disagree with physics
    fall back to RSSI, two gateways give a centroid, one gives nothing, and
    gateways left out of `positions` are ignored even when they lie loudest;
  * fine timestamps parse to the nanosecond, a gateway heard on two antennas
    is one reception, and an end-to-end solve from rxInfo works;
  * the CLI prints one line per device.
"""

import io
import json
import math
import random
import sys

import lorawan_locate as ll
from lorawan_locate import __main__ as cli

C = ll.SPEED_OF_LIGHT
GATEWAYS = [(0, 0), (3000, 0), (0, 3000), (3000, 3000), (1500, -1000)]
PROJ = ll.Projection(45.0, 10.0)
IDS = [f'aa555a000000000{i}' for i in range(5)]
POSITIONS = {gid: PROJ.latlon(*xy) for gid, xy in zip(IDS, GATEWAYS)}
TRUTH = (1200.0, 1700.0)


def rssi_at(tx, gw, power, rng, shadowing=4.0):
    d = max(math.hypot(tx[0] - gw[0], tx[1] - gw[1]), 10)
    return power - 10 * ll.PATH_LOSS_EXPONENT * math.log10(d) + rng.gauss(0, shadowing)


def arrival_ns(tx, gw):
    return 1_000_000_000_000 + round(math.dist(tx, gw) / C * 1e9)


def error_m(solved):
    return ll.distance_m(solved['latitude'], solved['longitude'], *PROJ.latlon(*TRUTH))


def rssi_uplinks(rng, count=10):
    return [[{'gateway_id': gid, 'rssi': rssi_at(TRUTH, xy, -10, rng), 'snr': 5.0, 'fine_ns': None}
             for gid, xy in zip(IDS, GATEWAYS)] for _ in range(count)]


def test_rssi_averages_uplinks():
    solved = ll.solve(rssi_uplinks(random.Random(7)), POSITIONS)
    assert solved['method'] == 'rssi', solved
    assert error_m(solved) <= 2 * solved['accuracyM'], solved
    assert solved['uplinks'] == 10 and solved['gatewayIds'] == IDS


def test_tdoa_wins_and_broken_clocks_fall_back():
    uplinks = rssi_uplinks(random.Random(7))
    timed = [{'gateway_id': gid, 'rssi': -100, 'snr': 1.0, 'fine_ns': arrival_ns(TRUTH, xy)}
             for gid, xy in zip(IDS[:4], GATEWAYS[:4])]
    solved = ll.solve([timed] + uplinks, POSITIONS)
    assert solved['method'] == 'tdoa' and solved['uplinks'] == 1
    assert error_m(solved) < 10, solved
    # 50 microseconds apart cannot be one frame heard across 3 km. With four
    # gateways the broken one is the one whose removal makes physics work...
    broken = [dict(r, fine_ns=r['fine_ns'] + (50_000 if i == 0 else 0)) for i, r in enumerate(timed)]
    solved = ll.solve([broken] + uplinks, POSITIONS)
    assert solved['method'] == 'tdoa' and solved['gatewayIds'] == IDS[1:4], solved
    assert error_m(solved) < 10, solved
    # ...with three there is no telling, and RSSI it is.
    assert ll.solve([broken[:3]] + uplinks, POSITIONS)['method'] == 'rssi'


def test_one_bad_clock_of_five_is_left_out():
    gateways = dict(zip(IDS, GATEWAYS))
    timed = [{'gateway_id': g, 'rssi': -100, 'fine_ns': arrival_ns(TRUTH, xy) + (2000 if g == IDS[2] else 0)}
             for g, xy in gateways.items()]
    solved = ll.solve([timed], POSITIONS)
    assert solved['method'] == 'tdoa' and IDS[2] not in solved['gatewayIds'], solved
    assert error_m(solved) < 10, solved
    # Four gateways, one 2 us off: any three fit, so nobody is left out and
    # the radius widens to admit the doubt.
    solved = ll.solve([timed[:4]], POSITIONS)
    assert solved['gatewayIds'] == IDS[:4] and error_m(solved) <= 2 * solved['accuracyM'], solved


def test_gateways_in_a_line():
    """Three gateways along a road: the device and its mirror image across
    the road fit equally well, so the answer is on the road, wide enough for
    both sides."""
    road = [(0, 0), (1500, 0), (3000, 0)]
    positions = {g: PROJ.latlon(*xy) for g, xy in zip(IDS, road)}
    for side in (800, -800):
        tx = (1000, side)
        timed = [{'gateway_id': g, 'rssi': rssi_at(tx, xy, -12, random.Random(1), 0),
                  'fine_ns': arrival_ns(tx, xy)} for g, xy in zip(IDS, road)]
        for uplink in (timed, [dict(r, fine_ns=None) for r in timed]):
            solved = ll.solve([uplink], positions)
            err = ll.distance_m(solved['latitude'], solved['longitude'], *PROJ.latlon(*tx))
            assert abs(PROJ.xy(solved['latitude'], solved['longitude'])[1]) < 1, solved
            assert err <= solved['accuracyM'], (solved, err)


def test_bad_input():
    pair = [[{'gateway_id': 'a', 'rssi': -80}, {'gateway_id': 'b', 'rssi': -90}]]
    good = PROJ.latlon(0, 0)
    for bad in [(float('nan'), 1), (None, 1), (200, 5), (0, 0), 'x', (1,)]:
        assert ll.solve(pair, {'a': bad, 'b': good}) is None, bad
    for rssi in ['x', None, float('inf'), True]:
        assert ll.solve([[{'gateway_id': 'a', 'rssi': rssi}, *pair[0][1:]]],
                        {'a': PROJ.latlon(3000, 0), 'b': good}) is None, rssi
    assert ll.solve([], {}) is None and ll.solve([[]], {'a': good}) is None
    # Three gateways on one mast are one place: nothing to solve from.
    mast = {g: good for g in 'abc'}
    assert ll.solve([[{'gateway_id': g, 'rssi': -90, 'fine_ns': 1} for g in 'abc']], mast) is None
    # Two on the mast and one across town: a centroid, not a confident TDOA.
    mast['c'] = PROJ.latlon(3000, 0)
    assert ll.solve([[{'gateway_id': g, 'rssi': -90, 'fine_ns': 1} for g in 'abc']], mast)['method'] == 'centroid'


def test_two_gateways_centroid_one_nothing():
    pair = [[{'gateway_id': IDS[0], 'rssi': -80}, {'gateway_id': IDS[1], 'rssi': -110}]]
    solved = ll.solve(pair, POSITIONS)
    assert solved['method'] == 'centroid'
    to_loud = ll.distance_m(solved['latitude'], solved['longitude'], *POSITIONS[IDS[0]])
    to_quiet = ll.distance_m(solved['latitude'], solved['longitude'], *POSITIONS[IDS[1]])
    assert to_loud < to_quiet and solved['accuracyM'] >= 1500 - 1
    assert ll.solve([[{'gateway_id': IDS[0], 'rssi': -80}]], POSITIONS) is None
    assert ll.solve(rssi_uplinks(random.Random(1)), {IDS[0]: POSITIONS[IDS[0]]}) is None


def test_untrusted_gateways_are_ignored():
    """A gateway the caller left out, claiming -30 dBm and timestamps that
    would pull a TDOA solve to its own point, moves nothing."""
    rng = random.Random(3)
    honest = dict(zip(IDS[:3], GATEWAYS[:3]))
    liars = {f'bad{i}': PROJ.latlon(*xy) for i, xy in enumerate([(5500, 1500), (5600, 1600), (5400, 1700)])}
    uplinks = []
    for _ in range(5):
        uplink = [{'gateway_id': g, 'rssi': rssi_at(TRUTH, xy, -12, rng), 'fine_ns': None}
                  for g, xy in honest.items()]
        uplink += [{'gateway_id': g, 'rssi': -30, 'fine_ns': i} for i, g in enumerate(liars)]
        uplinks.append(uplink)
    solved = ll.solve(uplinks, {g: POSITIONS[g] for g in honest})
    assert solved['method'] == 'rssi' and solved['gatewayIds'] == sorted(honest)
    assert error_m(solved) <= 2 * solved['accuracyM']
    # The same receptions, believing only the liars: their TDOA, their point.
    assert ll.solve(uplinks, liars)['method'] == 'tdoa'


def test_calibration_is_passed_through():
    uplinks = rssi_uplinks(random.Random(7))
    default = ll.solve(uplinks, POSITIONS)
    tuned = ll.solve(uplinks, POSITIONS, path_loss_exponent=2.0, shadowing_db=10.0)
    assert (tuned['latitude'], tuned['accuracyM']) != (default['latitude'], default['accuracyM'])


def test_gps_time():
    assert ll.gps_time_ns('1373000000.000000123s') == 1_373_000_000_000_000_123
    assert ll.gps_time_ns('1373000000.5s') == 1_373_000_000_500_000_000
    assert ll.gps_time_ns('1373000000s') == 1_373_000_000_000_000_000
    assert ll.gps_time_ns({'seconds': '12', 'nanos': 34}) == 12_000_000_034
    for bad in (None, 1373000000.5, '1373000000', '-1s', '1.0000000001s', 'xs'):
        assert ll.gps_time_ns(bad) is None, bad


def test_receptions_one_per_gateway():
    rows = ll.receptions([
        {'gatewayId': 'AA555A0000000001', 'rssi': -101, 'snr': 2.5,
         'fineTimeSinceGpsEpoch': '1373000000.000000123s'},
        {'gatewayId': 'aa555a0000000001', 'rssi': -97, 'snr': 3.0, 'antenna': 1},
        {'gatewayId': 'aa555a0000000002', 'rssi': -120, 'snr': -9.0},
        {'gatewayId': 'aa555a0000000003'},
        {'rssi': -50},
    ])
    assert rows == [
        # The louder antenna's RSSI, the other antenna's timestamp.
        {'gateway_id': 'aa555a0000000001', 'rssi': -97, 'snr': 3.0, 'fine_ns': 1_373_000_000_000_000_123},
        {'gateway_id': 'aa555a0000000002', 'rssi': -120, 'snr': -9.0, 'fine_ns': None},
    ]


def event(dev_eui, time, receptions):
    return {'time': time, 'deviceInfo': {'devEui': dev_eui}, 'rxInfo': receptions}


def chirpstack_events():
    """Four GPS gateways time-stamping one uplink of device 01, and an older
    uplink without; device 02 heard by one gateway."""
    rx = []
    for gid, xy in zip(IDS[:4], GATEWAYS[:4]):
        ns = arrival_ns(TRUTH, xy)
        lat, lon = POSITIONS[gid]
        rx.append({'gatewayId': gid.upper(), 'rssi': -100, 'snr': 2.0,
                   'fineTimeSinceGpsEpoch': f'{ns // 10**9}.{ns % 10**9:09d}s',
                   'location': {'latitude': lat, 'longitude': lon, 'altitude': 30}})
    old = [{k: v for k, v in r.items() if k != 'fineTimeSinceGpsEpoch'} for r in rx]
    lone = [{'gatewayId': IDS[0], 'rssi': -90, 'location': {'latitude': 0, 'longitude': 0}}]
    return [event('01', '2026-11-11T10:00:00Z', old), event('01', '2026-11-11T10:05:00Z', rx),
            event('02', '2026-11-11T10:05:00Z', lone)]


def test_solve_events():
    events = chirpstack_events()
    assert ll.gateway_positions(events) == {g: POSITIONS[g] for g in IDS[:4]}, 'Null Island skipped'
    solved = ll.solve_events(events[:2])
    assert solved['method'] == 'tdoa' and error_m(solved) < 10, solved
    # Positions given: only those gateways, wherever the events say they are.
    assert ll.solve_events(events[:2], {g: POSITIONS[g] for g in IDS[:2]})['method'] == 'centroid'


def test_cli(tmp_path, monkeypatch, capsys):
    path = tmp_path / 'events.jsonl'
    path.write_text('\n'.join(json.dumps(e) for e in chirpstack_events()))
    cli.main([str(path)])
    lines = [json.loads(line) for line in capsys.readouterr().out.splitlines()]
    assert [line['devEui'] for line in lines] == ['01', '02']
    assert lines[0]['location']['method'] == 'tdoa' and lines[1]['location'] is None

    gateways = tmp_path / 'gateways.json'
    gateways.write_text(json.dumps({IDS[0].upper(): POSITIONS[IDS[0]], IDS[1]: POSITIONS[IDS[1]]}))
    monkeypatch.setattr(sys, 'stdin', io.StringIO(json.dumps(chirpstack_events())))
    cli.main(['--gateways', str(gateways)])
    assert json.loads(capsys.readouterr().out.splitlines()[0])['location']['method'] == 'centroid'

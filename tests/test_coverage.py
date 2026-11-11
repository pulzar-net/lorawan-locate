"""Coverage cells (the [h3] extra): binning, exact merging, GeoJSON."""

from datetime import datetime, timedelta, timezone

import pytest

h3 = pytest.importorskip('h3')
from lorawan_locate import coverage  # noqa: E402

T0 = datetime(2026, 10, 13, 12, 0, tzinfo=timezone.utc)
LAT, LON = 45.0123, 10.0456


def cell_row(gateway_id, cell, samples, rssi_mean, rssi_max, at=T0):
    return {'gateway_id': gateway_id, 'cell': cell, 'samples': samples,
            'rssi_sum': rssi_mean * samples, 'rssi_max': rssi_max,
            'snr_sum': 5.0 * samples, 'snr_samples': samples, 'snr_max': 9.0,
            'first_at': at, 'last_at': at}


def test_binning():
    cell = coverage.cell_of(LAT, LON)
    assert h3.get_resolution(cell) == 9 and coverage.valid_cell(cell)
    assert coverage.cell_of(LAT + 0.00002, LON) == cell
    assert coverage.cell_of(LAT + 0.01, LON) != cell
    lat, lon = coverage.centre(cell)
    assert abs(lat - LAT) < 0.003 and abs(lon - LON) < 0.003
    ring = coverage.boundary(cell)
    assert ring[0] == ring[-1] and len(ring) == 7, 'a closed hexagon'
    assert abs(ring[0][0] - LON) < 0.01 and abs(ring[0][1] - LAT) < 0.01, '[lon, lat]'
    assert not coverage.valid_cell('nonsense') and not coverage.valid_cell(None)


def test_plausible():
    assert coverage.plausible(LAT, LON)
    assert not coverage.plausible(0, 0), 'Null Island is a receiver with no fix'
    assert not coverage.plausible(None, LON) and not coverage.plausible(91, LON)
    assert coverage.plausible(LAT, LON, (LAT + 0.5, LON))
    assert not coverage.plausible(LAT, LON, (LAT + 3.0, LON)), '330 km from its gateway'


def test_merge_is_exact():
    parent = h3.cell_to_parent(coverage.cell_of(LAT, LON), 8)
    siblings = list(h3.cell_to_children(parent, 9))[:3]
    rows = [
        cell_row('gw1', siblings[0], 4, -100.0, -90, at=T0),
        cell_row('gw1', siblings[1], 6, -110.0, -95, at=T0 + timedelta(hours=1)),
        cell_row('gw2', siblings[2], 2, -80.0, -78, at=T0 - timedelta(hours=1)),
    ]
    assert len(coverage.merge(rows, 9)) == 3
    coarse = coverage.merge(rows, 8)
    assert set(coarse) == {('gw1', parent), ('gw2', parent)}
    gw1 = coarse[('gw1', parent)]
    assert gw1['samples'] == 10 and gw1['rssi_sum'] == -400.0 - 660.0
    assert gw1['rssi_max'] == -90 and gw1['snr_samples'] == 10
    assert gw1['first_at'] == T0 and gw1['last_at'] == T0 + timedelta(hours=1)
    assert coverage.summary(gw1)['rssi'] == -106.0
    assert rows[0]['samples'] == 4, 'merge must not mutate its input'

    top, every = coverage.best(coarse)[parent]
    assert top['gatewayId'] == 'gw2' and [g['gatewayId'] for g in every] == ['gw2', 'gw1']

    collection = coverage.collection(coarse)
    feature = collection['features'][0]
    assert feature['id'] == parent and feature['geometry']['type'] == 'Polygon'
    assert feature['properties']['gatewayId'] == 'gw2' and feature['properties']['gateways'] == 2
    only = coverage.collection(coarse, gateway_id='gw1')['features'][0]['properties']
    assert only['gatewayId'] == 'gw1' and only['rssi'] == -106.0 and only['grade'] == 'fair'
    assert coverage.collection(coarse, gateway_id='nobody')['features'] == []


def test_grade():
    assert coverage.grade(-90, 5) == 'good' and coverage.grade(-90, -3) == 'fair'
    assert coverage.grade(-110, None) == 'fair' and coverage.grade(-120, -10) == 'weak'
    assert coverage.grade(None, 3) is None

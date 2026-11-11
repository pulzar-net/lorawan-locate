"""Coverage cells: where each gateway was heard from, binned and summarised.

H3 rather than a geohash grid. A geohash cell is a rectangle whose width
shrinks with latitude, so the same "cell" is a different area in Austin and in
Anchorage and a neighbour is one of eight at two different distances. H3's
hexagons are close to equal area everywhere and every neighbour is equally far,
which is what a coverage map is drawn and compared on. The Python binding is
one wheel with no dependencies (`pip install lorawan-locate[h3]`); a map needs
no H3 at all, because it is sent each cell's outline as GeoJSON.

Samples are binned at resolution 9 (edges of about 200 m, a short walk) and
merged into parents for views further out. Merging is exact for counts, sums
and maxima, which is why a cell stores sums rather than means.

Pure functions; storing the cells is the caller's.

A cell row is a dict: gateway_id, cell, samples, rssi_sum, rssi_max, and
optionally snr_sum, snr_samples, snr_max, first_at, last_at (datetimes or ISO
strings), plus anything else the caller keeps, which merge() carries along.
"""

import h3

from . import distance_m

RESOLUTION = 9
RESOLUTIONS = (6, 7, 8, 9)

# A gateway heard from further than this is a GNSS fix or a phone gone wrong,
# or a sample somebody made up; LoRa's records are longer, but not from a
# sensor that stays put.
MAX_SAMPLE_RANGE_M = 100_000.0
# A phone fix less precise than this places a sample in the wrong cell.
MAX_FIX_ACCURACY_M = 100.0


def cell_of(latitude, longitude, resolution=RESOLUTION):
    return h3.latlng_to_cell(latitude, longitude, resolution)


def valid_cell(cell):
    try:
        return isinstance(cell, str) and h3.is_valid_cell(cell)
    except (TypeError, ValueError):
        return False


def centre(cell):
    return h3.cell_to_latlng(cell)


def plausible(latitude, longitude, gateway=None):
    """Whether a sample's position is worth keeping.

    Null Island is a receiver with no fix. `gateway`, the (lat, lon) of the
    gateway that heard it when known, bounds how far away the device can have
    been.
    """
    if latitude is None or longitude is None:
        return False
    if not (-90 <= latitude <= 90 and -180 <= longitude <= 180):
        return False
    if abs(latitude) < 1e-6 and abs(longitude) < 1e-6:
        return False
    if gateway is not None:
        return distance_m(latitude, longitude, *gateway) <= MAX_SAMPLE_RANGE_M
    return True


def merge(cells, resolution):
    """Cells keyed (gateway_id, cell) at `resolution`, from finer cells.

    Each input is a coverage_cells row as a dict. Counts and sums add,
    maxima take the max, and first/last widen.
    """
    out = {}
    for row in cells:
        cell = row['cell']
        if h3.get_resolution(cell) > resolution:
            cell = h3.cell_to_parent(cell, resolution)
        key = (row['gateway_id'], cell)
        into = out.get(key)
        if into is None:
            out[key] = {**row, 'cell': cell}
            continue
        into['samples'] += row['samples']
        into['rssi_sum'] += row['rssi_sum']
        into['rssi_max'] = max(into['rssi_max'], row['rssi_max'])
        into['snr_sum'] = (into.get('snr_sum') or 0) + (row.get('snr_sum') or 0)
        into['snr_samples'] = (into.get('snr_samples') or 0) + (row.get('snr_samples') or 0)
        if row.get('snr_max') is not None:
            into['snr_max'] = row['snr_max'] if into.get('snr_max') is None else max(into['snr_max'], row['snr_max'])
        if row.get('first_at') is not None:
            into['first_at'] = min(filter(None, (into.get('first_at'), row['first_at'])))
        if row.get('last_at') is not None:
            into['last_at'] = max(filter(None, (into.get('last_at'), row['last_at'])))
    return out


def summary(row):
    """One (gateway, cell) aggregate as the API reports it."""
    snr_samples = row.get('snr_samples') or 0
    return {
        'gatewayId': row['gateway_id'],
        'samples': row['samples'],
        'rssi': round(row['rssi_sum'] / row['samples'], 1),
        'rssiMax': row['rssi_max'],
        'snr': round(row['snr_sum'] / snr_samples, 1) if snr_samples else None,
        'snrMax': row.get('snr_max'),
        'firstAt': _iso(row.get('first_at')),
        'lastAt': _iso(row.get('last_at')),
    }


def best(merged):
    """{cell: (summary of the best gateway, every gateway's summary)}.

    Best is the highest mean RSSI: the gateway a device in that cell would
    most likely be heard by, which is the question "is there coverage here"
    asks.
    """
    by_cell = {}
    for (_gateway, cell), row in merged.items():
        by_cell.setdefault(cell, []).append(summary(row))
    out = {}
    for cell, gateways in by_cell.items():
        gateways.sort(key=lambda g: (-g['rssi'], -g['samples'], g['gatewayId']))
        out[cell] = (gateways[0], gateways)
    return out


def boundary(cell):
    """The cell's outline as a closed GeoJSON ring, [lon, lat]."""
    ring = [[round(lon, 6), round(lat, 6)] for lat, lon in h3.cell_to_boundary(cell)]
    return ring + [ring[0]]


def feature(cell, properties):
    return {
        'type': 'Feature',
        'id': cell,
        'geometry': {'type': 'Polygon', 'coordinates': [boundary(cell)]},
        'properties': {'cell': cell, **properties},
    }


def collection(merged, gateway_id=None):
    """A GeoJSON FeatureCollection of merge()'s cells: the best gateway's
    summary per cell, or one gateway's with `gateway_id`."""
    features = []
    for cell, (top, gateways) in best(merged).items():
        if gateway_id:
            top = next((x for x in gateways if x['gatewayId'] == gateway_id), None)
            if top is None:
                continue
        features.append(feature(cell, {**top, 'grade': grade(top['rssi'], top['snr']),
                                       'gateways': len(gateways)}))
    features.sort(key=lambda f: f['id'])
    return {'type': 'FeatureCollection', 'features': features}


def grade(rssi, snr):
    """A word for the signal, on the thresholds a technician already uses
    for LoRa: above -100 dBm is good, below -115 is marginal."""
    if rssi is None:
        return None
    if rssi > -100 and (snr is None or snr > 0):
        return 'good'
    if rssi > -115:
        return 'fair'
    return 'weak'


def _iso(value):
    return value.isoformat() if hasattr(value, 'isoformat') else value

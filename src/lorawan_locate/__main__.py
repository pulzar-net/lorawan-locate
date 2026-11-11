"""lorawan-locate [EVENTS] [--gateways FILE]

Reads ChirpStack v4 uplink events -- one JSON object, an array, or one per
line -- from EVENTS or stdin, and prints one JSON line per device: its solved
position, or null when the gateways that heard it cannot place it.
"""

import argparse
import json
import sys

from . import solve_events


def _events(text):
    text = text.strip()
    if not text:
        return []
    try:
        parsed = json.loads(text)
    except json.JSONDecodeError:
        return [json.loads(line) for line in text.splitlines() if line.strip()]
    return parsed if isinstance(parsed, list) else [parsed]


def main(argv=None):
    parser = argparse.ArgumentParser(prog='lorawan-locate', description=__doc__.split('\n\n')[1])
    parser.add_argument('events', nargs='?', type=argparse.FileType('r'), default=sys.stdin)
    parser.add_argument('--gateways', type=argparse.FileType('r'),
                        help='JSON {gatewayId: [lat, lon]}: use only these gateways, at these '
                             'positions, instead of the locations in rxInfo')
    args = parser.parse_args(argv)
    positions = None
    if args.gateways:
        positions = {g.lower(): tuple(p) for g, p in json.load(args.gateways).items()}
    by_device = {}
    for event in _events(args.events.read()):
        dev_eui = (event.get('deviceInfo') or {}).get('devEui')
        by_device.setdefault(dev_eui, []).append(event)
    for dev_eui, events in by_device.items():
        print(json.dumps({'devEui': dev_eui, 'location': solve_events(events, positions)}))


if __name__ == '__main__':
    main()

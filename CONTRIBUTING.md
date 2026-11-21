# Contributing

Issues and pull requests are welcome.

## Setup

```bash
python -m venv .venv
.venv/bin/pip install -e '.[h3]' pytest
.venv/bin/python -m pytest
```

## Ground rules

- **Core stays standard library only.** Anything that needs a dependency goes behind an extra, the way H3 is.
- **Accuracy claims need evidence.** A change that affects accuracy should come with the before and after output of `python tests/test_accuracy.py`. If it changes the numbers, update the README table and `BOUNDS` in the same pull request. Loosening a bound needs a reason.
- **Prefer an honest radius to a small one.** A solve that is wrong should say so through `accuracyM`.
- **Every behaviour change needs a test,** in `tests/test_locate.py` for the solver and input handling, or in `tests/test_coverage.py`.

## Real data

The most useful thing you can contribute is a measurement against real receptions: devices at surveyed or GNSS-known positions, heard by gateways at known positions. Please open an issue with:

- how many gateways there were, their models, and whether they have GPS fine timestamps;
- how the true positions were known, and how precisely;
- the error you measured (median, 90th percentile, and how often it was within `accuracyM`) and the calibration you used.

Please don't post device or gateway identifiers or exact coordinates you aren't free to share. Anonymised offsets in metres are enough.

## Releases

Bump `__version__` in `src/lorawan_locate/__init__.py`, then push an annotated tag
whose message is the release notes (`git tag -a vX.Y.Z --cleanup=whitespace -F notes.md`).
There is no changelog file: what changed is in the
[GitHub releases](https://github.com/pulzar-net/lorawan-locate/releases).

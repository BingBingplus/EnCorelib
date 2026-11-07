"""Reader for the box files written by EnCorelib (E0.txt, E1.txt, E_0.txt, ...).

A box file has one box per line:

    Box 0: [lo, hi] x [lo, hi] [x [lo, hi]]

Numbers may be in scientific notation.  `read_boxes` returns a list of
boxes; each box is a list of (lo, hi) pairs, one per coordinate.
"""

import re

_NUM = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?"
_PAIR = re.compile(r"\[\s*(" + _NUM + r")\s*,\s*(" + _NUM + r")\s*\]")


def read_boxes(path):
    boxes = []
    with open(path, "r") as fh:
        for line in fh:
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            pairs = _PAIR.findall(line)
            if pairs:
                boxes.append([(float(a), float(b)) for a, b in pairs])
    return boxes


def bounds(boxes, i):
    """(min lo, max hi) of coordinate i over a list of boxes."""
    lo = min(b[i][0] for b in boxes)
    hi = max(b[i][1] for b in boxes)
    return lo, hi


def pad(lo, hi, frac=0.05):
    span = hi - lo
    if span <= 0:
        span = 1.0
    return lo - frac * span, hi + frac * span

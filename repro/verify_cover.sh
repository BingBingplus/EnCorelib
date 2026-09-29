#!/usr/bin/env bash
# ===========================================================================
#  repro/verify_cover.sh  --  independent check of the computed cover
#
#  The paper guarantees  End(B0,H) subset Union(C).  This script tests that
#  guarantee from the outside: it takes a grid of sample points of B0,
#  integrates each of them separately with ./Encoretraj (a plain validated
#  integration), and
#  checks that every resulting enclosure lies inside some box of C.
#
#
#  Usage
#      bash repro/verify_cover.sh                 all examples, all methods
#      bash repro/verify_cover.sh eg1 eg5
#      GRID=7 EPS=0.1 bash repro/verify_cover.sh eg2
#
#  Requires:  make          (./EnCorelib)
#             make Encoretraj  (./Encoretraj)
#             python3
# ===========================================================================
set -u

cd "$(dirname "$0")/.." || exit 1
BIN=./EnCorelib;   [ -x "$BIN" ]  || BIN=./EnCorelib.exe
[ -x "$BIN" ]   || BIN=./bin/EnCorelib.exe
TRAJ=./Encoretraj; [ -x "$TRAJ" ] || TRAJ=./Encoretraj.exe
[ -x "$TRAJ" ]  || TRAJ=./bin/Encoretraj.exe
PY=${PYTHON:-python3}
command -v "$PY" >/dev/null 2>&1 || PY=python
for f in "$BIN" "$TRAJ"; do
  [ -x "$f" ] || { echo "build first (or use bin/):  make && make Encoretraj" >&2; exit 1; }
done

GRID=${GRID:-4}          # samples per coordinate
EPS=${EPS:-0.3}
METHODS=${METHODS:-"0 1 2"}
OUT=repro/out/verify
mkdir -p "$OUT"

#      tag | n | vars | rhs | box | H
eg1="eg1|2|x y|2*x-2*x*y -y+x*y|0.9 1.1 2.9 3.1|4"
eg2="eg2|2|x y|y (1-x^2)*y-x|-3.1 -2.9 2.9 3.1|1"
eg3="eg3|2|x y|y x^2|0.95 1.05 -1.05 -0.95|1"
eg4="eg4|2|x y|x-x^3/3-y+0.5 0.08*(x+0.7-0.8*y)|0.9 1.1 -0.1 0.1|1"
eg5="eg5|3|x y z|10*(y-x) x*(28-z)-y x*y-8*z/3|14.999 15.001 14.999 15.001 35.999 36.001|1"
eg6="eg6|3|x y z|-y-z x+0.2*y 0.2+z*(x-5.7)|0.9 1.1 1.9 2.1 2.9 3.1|1"

WANT=${*:-"eg1 eg2 eg3 eg4 eg5 eg6"}
rc=0

for key in $WANT; do
  spec=$(eval "printf '%s' \"\${$key:-}\"")
  [ -n "$spec" ] || { echo "unknown example: $key" >&2; continue; }
  IFS='|' read -r tag dim vars rhs box H <<EOF
$spec
EOF

  # --- validated images of a grid of sample points of B0 -----------------
  pts=$OUT/$tag.pts
  : > "$pts"
  "$PY" - "$dim" "$GRID" $box <<'PYEOF' | while read -r line; do
import sys, itertools
dim  = int(sys.argv[1]); grid = int(sys.argv[2])
b    = [float(v) for v in sys.argv[3:]]
axes = []
for i in range(dim):
    lo, hi = b[2*i], b[2*i+1]
    axes.append([lo] if hi == lo else
                [lo + (hi-lo)*k/(grid-1) for k in range(grid)])
for p in itertools.product(*axes):
    print(" ".join(repr(v) for v in p))
PYEOF
    args=""
    for v in $line; do args="$args $v $v"; done
    "$TRAJ" "$dim" $vars $rhs 20 "$H" 0.05 $args 2>/dev/null | tail -1 >> "$pts"
  done

  for M in $METHODS; do
    ( cd "$OUT" && "../../../$BIN" 3 0 "$M" 0 0 "$dim" $vars $rhs \
        "$EPS" 20 "$H" 0 $box >/dev/null 2>&1 ) || true
    mv -f "$OUT/E1.txt" "$OUT/$tag.E1.m$M" 2>/dev/null || true

    res=$("$PY" - "$OUT/$tag.E1.m$M" "$pts" <<'PYEOF'
import re, sys
NUM  = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?"
PAIR = re.compile(r"\[\s*(" + NUM + r")\s*,\s*(" + NUM + r")\s*\]")
boxes = []
for line in open(sys.argv[1]):
    p = PAIR.findall(line)
    if p: boxes.append([(float(a), float(b)) for a, b in p])
pts = []
for line in open(sys.argv[2]):
    v = [float(t) for t in line.split()]
    pts.append([(v[1+2*i], v[2+2*i]) for i in range((len(v)-1)//2)])
bad = 0
for q in pts:
    if not any(all(b[i][0] <= q[i][0] and q[i][1] <= b[i][1] for i in range(len(q)))
               for b in boxes):
        bad += 1
print("%d %d %d" % (len(boxes), bad, len(pts)))
PYEOF
)
    set -- $res
    if [ "$2" = "0" ]; then
      printf 'ok    %-4s method=%s  #(C)=%-5s  all %s sample images covered\n' "$tag" "$M" "$1" "$3"
    else
      printf 'FAIL  %-4s method=%s  #(C)=%-5s  %s of %s sample images NOT covered\n' \
             "$tag" "$M" "$1" "$2" "$3"
      rc=1
    fi
  done
done

exit $rc

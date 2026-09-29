#!/usr/bin/env bash
# ===========================================================================
#  repro/figure4.sh  --  data for Figure 4 ("Geometry of EndCovers")
#
#  For each example of Table 1 this runs EndCover at output level 3, which
#  writes
#       E0.txt   the initial sub-boxes that were integrated
#       E1.txt   the eps-cover C at time H          (one box per E0 box)
#  into repro/out/fig4/<tag>/.
#
#  Draw them with
#       python tools/plot_boxes.py repro/out/fig4/eg1 --title "Eg1, H=4"
#
#  In 2D the script produces the two panels of Figure 4 (subdivision of B0
#  on top, the end-cover below); in 3D it produces the xy-, xz- and
#  yz-projections.
#
#  Usage
#      bash repro/figure4.sh                 all six examples
#      bash repro/figure4.sh eg1 eg5         only those
#      EPS=0.05 bash repro/figure4.sh eg1
# ===========================================================================
set -u

cd "$(dirname "$0")/.." || exit 1
BIN=./EnCorelib
[ -x "$BIN" ] || BIN=./EnCorelib.exe
[ -x "$BIN" ] || BIN=./bin/EnCorelib.exe
if [ ! -x "$BIN" ]; then echo "build first (or use bin/):  make" >&2; exit 1; fi

EPS=${EPS:-0.1}
METHOD=${METHOD:-0}

#      tag | n | vars | rhs | box | H
eg1="eg1|2|x y|2*x-2*x*y -y+x*y|0.9 1.1 2.9 3.1|4"
eg2="eg2|2|x y|y (1-x^2)*y-x|-3.1 -2.9 2.9 3.1|2"
eg3="eg3|2|x y|y x^2|0.95 1.05 -1.05 -0.95|4"
eg4="eg4|2|x y|x-x^3/3-y+0.5 0.08*(x+0.7-0.8*y)|0.9 1.1 -0.1 0.1|1"
eg5="eg5|3|x y z|10*(y-x) x*(28-z)-y x*y-8*z/3|14.999 15.001 14.999 15.001 35.999 36.001|4"
eg6="eg6|3|x y z|-y-z x+0.2*y 0.2+z*(x-5.7)|0.9 1.1 1.9 2.1 2.9 3.1|1"

WANT=${*:-"eg1 eg2 eg3 eg4 eg5 eg6"}

for key in $WANT; do
  spec=$(eval "printf '%s' \"\${$key:-}\"")
  if [ -z "$spec" ]; then echo "unknown example: $key" >&2; continue; fi
  IFS='|' read -r tag dim vars rhs box H <<EOF
$spec
EOF
  dir=repro/out/fig4/$tag
  mkdir -p "$dir"
  echo "== $tag  (n=$dim, H=$H, eps=$EPS) =="
  ( cd "$dir" && "../../../../$BIN" 3 0 "$METHOD" 0 0 "$dim" $vars $rhs \
       "$EPS" 20 "$H" 0 $box ) || echo "  FAILED"
  printf 'example=%s n=%s H=%s eps=%s method=%s\nvars=%s\nrhs=%s\nB0=%s\n' \
         "$tag" "$dim" "$H" "$EPS" "$METHOD" "$vars" "$rhs" "$box" > "$dir/run.info"
done

echo
echo "data written under repro/out/fig4/"
echo "plot with:  python tools/plot_boxes.py repro/out/fig4/eg1"

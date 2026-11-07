#!/usr/bin/env bash
# ===========================================================================
#  src/doc/compare_options.sh  --  the measurements behind section 6 of the user
#                               guide ("Choosing the options")
#
#      bash src/doc/compare_options.sh              one run at a time
#      JOBS=6 bash src/doc/compare_options.sh       six runs at a time
#      DRY=1 bash src/doc/compare_options.sh        list the runs, run nothing
#      ONLY=combo bash src/doc/compare_options.sh   only the runs whose tag starts so
#                                               (the guide's last table, two at a time)
#      LIMIT=1200 bash src/doc/compare_options.sh   stop each run after 1200 s (default 600)
#      python3 src/doc/make_tables.py               the LaTeX tables of section 6
#
#  Each of the six problems of Table 1 is run at its longer horizon with one
#  option changed at a time from the defaults (mode=0 method=0 stepB=0
#  tubedegree=0 order=20).  Lorenz is also run to T=5 with every stepB, and
#  every problem once more with method=2 and stepB=3 together.
#  Every run writes doc/out/options/<tag>.out (standard output) and .err.
#  A run still going after LIMIT seconds is stopped, and its .out file ends
#  with the line "stopped after LIMIT s".
# ===========================================================================
set -u

cd "$(dirname "$0")/../.." || exit 1
make --no-print-directory -s EndCorelib >/dev/null || { echo "build first: make" >&2; exit 1; }
BIN=$PWD/EndCorelib
OUT=doc/out/options
mkdir -p "$OUT"

#  problem | n vars f_1 .. f_n | initial box | horizon
#  (the right-hand sides multiplied out, see src/examples/_template.txt)
declare -A RHS BOX HOR
RHS[eg1]="2 x y 2*x-2*x*y -y+x*y";                         BOX[eg1]="0.9 1.1 2.9 3.1";                           HOR[eg1]=5.5
RHS[eg2]="2 x y y -x+y-x^2*y";                             BOX[eg2]="-3.1 -2.9 2.9 3.1";                         HOR[eg2]=2
RHS[eg3]="2 x y y x^2";                                    BOX[eg3]="0.95 1.05 -1.05 -0.95";                     HOR[eg3]=4
RHS[eg4]="2 x y 0.5+x-y+(-1/3)*x^3 0.056+0.08*x-0.064*y";  BOX[eg4]="0.9 1.1 -0.1 0.1";                          HOR[eg4]=4
RHS[eg5]="3 x y z -10*x+10*y 28*x-y-x*z (-8/3)*z+x*y";     BOX[eg5]="14.999 15.001 14.999 15.001 35.999 36.001"; HOR[eg5]=4
RHS[eg6]="3 x y z -y-z x+0.2*y 0.2-5.7*z+x*z";             BOX[eg6]="0.9 1.1 1.9 2.1 2.9 3.1";                   HOR[eg6]=4

# job TAG MODE METHOD STEPB TUBEDEGREE PROBLEM EPS ORDER [HORIZON]
job() { echo "$1|2 $2 $3 $4 $5 ${RHS[$6]} $7 $8 ${9:-${HOR[$6]}} 0 ${BOX[$6]}"; }

jobs_list() {
  for e in eg1 eg2 eg3 eg4 eg5 eg6; do
    for s in 0 2; do job "eps_${e}_e1.0_s$s" 0 0 $s 0 $e 1.0 20; done
    for eps in 0.5 0.1; do
      for s in 0 1 2 3; do job "stepB_${e}_e${eps}_s$s" 0 0 $s 0 $e $eps 20; done
      for m in 1 2;     do job "method_${e}_e${eps}_m$m" 0 $m 0 0 $e $eps 20; done
    done
    job "mode1_${e}_e0.1" 1 0 0 0 $e 0.1 20
  done
  for o in 10 15 30; do
    job "order_eg1_e0.5_o$o" 0 0 0 0 eg1 0.5 $o
    job "order_eg5_e0.5_o$o" 0 0 0 0 eg5 0.5 $o
    job "order_eg6_e0.5_o$o" 0 0 0 0 eg6 0.5 $o
  done
  for s in 0 1 2 3; do job "lorenzT5_s$s" 0 0 $s 0 eg5 1.0 20 5; done
  # the two faster settings together, method=2 and stepB=3, next to the
  # settings of the example files; run two at a time for the guide (JOBS=2)
  for e in eg1 eg2 eg3 eg4 eg5 eg6; do
    job "combo0_${e}_e0.1" 0 0 0 0 $e 0.1 20
    job "combo_${e}_e0.1"  0 2 3 0 $e 0.1 20
  done
  job "combo0_lorenzT5" 0 0 0 0 eg5 1.0 20 5
  job "combo_lorenzT5"  0 2 3 0 eg5 1.0 20 5
}

# ONLY=prefix runs just the runs whose tag starts with prefix, e.g. ONLY=combo
selected() { if [ -n "${ONLY:-}" ]; then grep "^${ONLY}"; else cat; fi; }

if [ "${DRY:-0}" = "1" ]; then jobs_list | selected; exit 0; fi

LIMIT=${LIMIT:-600}

run_one() {
  set -f
  local tag=${1%%|*} args=${1#*|} rc
  timeout "$LIMIT" "$BIN" $args > "$OUT/$tag.out" 2> "$OUT/$tag.err"
  rc=$?
  [ "$rc" -eq 124 ] && echo "stopped after $LIMIT s" >> "$OUT/$tag.out"
  echo "  $tag  exit=$rc"
}
export -f run_one
export BIN OUT LIMIT

jobs_list | selected | tr '\n' '\0' | xargs -0 -P "${JOBS:-1}" -I{} bash -c 'run_one "$1"' _ {}
echo "done; now:  python3 src/doc/make_tables.py"

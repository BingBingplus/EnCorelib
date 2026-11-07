#!/usr/bin/env bash
# ===========================================================================
#  doc/run_examples.sh  --  regenerate the program output shown in the user
#                            guide (doc/userguide.tex)
#
#      bash doc/run_examples.sh            the short examples (about 5 min)
#      LONG=1 bash doc/run_examples.sh     also the Lorenz runs (5 more min)
#
#  Every example writes doc/out/<name>.txt: the command as it is typed,
#  after "$ ", followed by what it printed.  The guide includes these files
#  verbatim, so running this script after a change to the code brings the
#  examples of the guide up to date; rebuild the PDF afterwards.
#
#  The examples run in the top directory like any other run, so they leave
#  E0.txt, E1.txt, ... there and append rows to out.txt.
# ===========================================================================
set -u
export LC_ALL=C           # messages of make and the shell in English

cd "$(dirname "$0")/.." || exit 1
make --no-print-directory -s EnCorelib >/dev/null || { echo "build first: make" >&2; exit 1; }

OUT=doc/out
mkdir -p "$OUT"

# make reports which CAPD and SymEngine it found; those lines name
# directories of this machine and are left out of the examples.
quiet_make() { grep -v -e '^CAPD source tree:' -e '^CAPD via ' -e '^SymEngine:'; }

# ex  NAME 'COMMAND'    keep standard output
# ex2 NAME 'COMMAND'    keep standard output and standard error, interleaved
ex()  { printf '$ %s\n' "$2" > "$OUT/$1.txt"; eval "$2" 2>/dev/null | quiet_make >> "$OUT/$1.txt"; echo "  $OUT/$1.txt"; }
ex2() { printf '$ %s\n' "$2" > "$OUT/$1.txt"; eval "$2" 2>&1      | quiet_make >> "$OUT/$1.txt"; echo "  $OUT/$1.txt"; }

LV='2 0 0 0 0 2 x y 2*x-2*x*y -y+x*y'        # Eg1, all options at their defaults
LVBOX='0.9 1.1  2.9 3.1'

echo "short examples"
ex  check          'make check'
ex2 usage          './EnCorelib'
ex  quickstart     'make run-eg1Volterra'
ex  showargs       'make show-args EG=src/examples/eg1Volterra.mk'
ex  cli_eg1        "./EnCorelib $LV 1.0 20 4 0  $LVBOX"
ex  lorenz_T2      'make run EG=src/examples/eg5Lorenz.mk T=2'
ex  boundary       'make run-boundary-eg1Volterra'

# a problem of one's own (section 4 of the guide)
ex  pendulum       'make run EG=src/examples/dampedPendulum.mk'
ex2 pendulum_debug 'make run EG=src/examples/dampedPendulum.mk eps=0.05 debug=1'
ex  sci            './EnCorelib 1 0 0 0 0 2 x y y "-sin(x)-2e-1*y" 0.1 20 10 0  0.9 1.1  -0.1 0.1'
ex2 parse_error    './EnCorelib 1 0 0 0 0 2 x y y "-sin(x)-0.2*Y" 0.1 20 10 0  0.9 1.1  -0.1 0.1'

# progress reports (section 5)
ex2 debug1         "./EnCorelib $LV 0.5 20 4 1  $LVBOX"
ex2 debug2         "./EnCorelib $LV 0.5 20 4 2  $LVBOX"
ex2 debug3         "./EnCorelib $LV 0.5 20 4 3  $LVBOX"

# output files (section 3) and the data of the two figures
ex  files          "./EnCorelib 3 0 0 0 0 2 x y 2*x-2*x*y -y+x*y 0.1 20 4 0  $LVBOX"
mkdir -p "$OUT/eg1_eps0.1" && cp E0.txt E1.txt convex_hull.txt "$OUT/eg1_eps0.1/"
ex  files_head     'head -n 3 E1.txt'
ex  files_log      'tail -n 2 out.txt'
ex  eijgenraam     'make run EG=src/examples/fig1Eijgenraam.mk'
mkdir -p "$OUT/eijgenraam" && cp E0.txt E1.txt "$OUT/eijgenraam/"

# the box-by-box check (section 8)
ex  audit          'make audit'

if [ "${LONG:-0}" = "1" ]; then
  echo "long examples"
  ex2 lorenz_T5        'make run EG=src/examples/eg5Lorenz.mk T=5 debug=1'
  ex2 lorenz_T5_stepB3 'make run EG=src/examples/eg5Lorenz.mk T=5 stepB=3 debug=1'
  ex2 lorenz_T100      'make run EG=src/examples/eg5Lorenz.mk T=100 stepB=3 debug=1'
fi

echo "figures:  python3 doc/make_figures.py"

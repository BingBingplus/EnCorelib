#!/usr/bin/env bash
# ===========================================================================
#  repro/check.sh  --  fast self-test (a few seconds)
#
#  Runs a handful of configurations and compares B1 (rounded to two decimals,
#  as printed in Table 2) and #(C) against known-good values.  Use it after
#  building, or after changing the code, to confirm nothing has drifted.
#
#      bash repro/check.sh          (or:  make check)
#
#  BASELINE=1 bash repro/check.sh   prints what the code currently produces,
#                                   in the form the `check` calls below take,
#                                   so the expectations can be re-baselined
#                                   after an intentional change.
# ===========================================================================
set -u

cd "$(dirname "$0")/.." || exit 1
BIN=./EnCorelib
[ -x "$BIN" ] || BIN=./EnCorelib.exe
[ -x "$BIN" ] || BIN=./bin/EnCorelib.exe
if [ ! -x "$BIN" ]; then echo "build first (or use bin/):  make" >&2; exit 1; fi

BASELINE=${BASELINE:-0}
pass=0; fail=0

fmt_awk='
function r2(v,   s) { s = sprintf("%.2f", v); return (s == "-0.00") ? "0.00" : s }
{ half=""; mid="";
  for (i = 1; i <= NF; i += 2) {
    lo = $i + 0; hi = $(i+1) + 0;
    mid  = mid  (i > 1 ? ", " : "") r2((lo+hi)/2);
    half = half (i > 1 ? ", " : "") r2((hi-lo)/2);
  }
  printf "(%s) +- (%s)", mid, half }'

# check <label> <expected B1> <expected #C> -- <EnCorelib arguments...>
check() {
  local label="$1" wantB1="$2" wantC="$3"; shift 4
  local out hull gotB1 gotC
  out=$("$BIN" "$@" 2>&1)
  hull=$(printf '%s' "$out" | sed -n 's/^Hull(T=[^)]*)=//p')
  gotC=$(printf '%s' "$out" | sed -n 's/^#(C)=\([0-9]*\).*/\1/p')
  if [ -z "$hull" ]; then
    printf 'FAIL  %-36s  no output\n' "$label"; fail=$((fail+1)); return
  fi
  gotB1=$(printf '%s' "$hull" | tr -d '{}[]' | awk -F, "$fmt_awk")
  if [ "$BASELINE" = "1" ]; then
    printf '%-36s  "%s"  %s\n' "$label" "$gotB1" "$gotC"; return
  fi
  if [ "$gotB1" = "$wantB1" ] && [ "$gotC" = "$wantC" ]; then
    printf 'ok    %-36s  %s  #(C)=%s\n' "$label" "$gotB1" "$gotC"; pass=$((pass+1))
  else
    printf 'FAIL  %-36s\n        expected  %s  #(C)=%s\n        got       %s  #(C)=%s\n' \
           "$label" "$wantB1" "$wantC" "$gotB1" "$gotC"; fail=$((fail+1))
  fi
}

echo "the default: mode=EndCover, method=0, stepB=2, order=20"
check "Eg1 Volterra  H=4  eps=1.0"  "(1.48, 0.19) +- (0.36, 0.04)"   5  -- \
      2 0 0 2 0 2 x y 2*x-2*x*y -y+x*y 1.0 20 4 0  0.9 1.1 2.9 3.1
check "Eg2 VanDerPol H=1  eps=1.0"  "(-2.14, 0.57) +- (0.27, 0.22)"   1  -- \
      2 0 0 2 0 2 x y y "(1-x^2)*y-x" 1.0 20 1 0  -3.1 -2.9 2.9 3.1
check "Eg3 Quadratic H=1  eps=0.1"  "(0.28, -0.58) +- (0.15, 0.15)"   21  -- \
      2 0 0 2 0 2 x y y "x^2" 0.1 20 1 0  0.95 1.05 -1.05 -0.95
check "Eg4 FitzHugh  H=1  eps=0.1"  "(1.76, 0.17) +- (0.08, 0.10)"   20  -- \
      2 0 0 2 0 2 x y "x-x^3/3-y+0.5" "0.08*(x+0.7-0.8*y)" 0.1 20 1 0  0.9 1.1 -0.1 0.1
check "Eg5 Lorenz    H=1  eps=1.0"  "(-6.95, 3.00, 35.14) +- (0.03, 0.01, 0.04)"   1  -- \
      2 0 0 2 0 3 x y z "10*(y-x)" "x*(28-z)-y" "x*y-8*z/3" 1.0 20 1 0 \
      14.999 15.001 14.999 15.001 35.999 36.001
check "Eg6 Roessler  H=1  eps=1.0"  "(-1.74, 1.86, 0.03) +- (0.15, 0.18, 0.00)"   1  -- \
      2 0 0 2 0 3 x y z -y-z "x+0.2*y" "0.2+z*(x-5.7)" 1.0 20 1 0 \
      0.9 1.1 1.9 2.1 2.9 3.1

echo
echo "method=2 (Split only), stepB=0"
check "Eg1 Volterra  H=4  eps=1.0" "(1.47, 0.19) +- (0.38, 0.05)"   5  -- \
      2 0 2 0 0 2 x y 2*x-2*x*y -y+x*y 1.0 20 4 0  0.9 1.1 2.9 3.1
check "Eg1 Volterra  H=4  eps=0.1" "(1.46, 0.19) +- (0.26, 0.03)"  85  -- \
      2 0 2 0 0 2 x y 2*x-2*x*y -y+x*y 0.1 20 4 0  0.9 1.1 2.9 3.1
check "Eg2 VanDerPol H=1  eps=1.0" "(-2.14, 0.57) +- (0.29, 0.29)"  1  -- \
      2 0 2 0 0 2 x y y "(1-x^2)*y-x" 1.0 20 1 0  -3.1 -2.9 2.9 3.1
check "Eg3 Quadratic H=1  eps=0.1" "(0.28, -0.58) +- (0.15, 0.15)" 21  -- \
      2 0 2 0 0 2 x y y "x^2" 0.1 20 1 0  0.95 1.05 -1.05 -0.95
check "Eg5 Lorenz    H=1  eps=1.0" \
      "(-6.95, 3.00, 35.14) +- (0.03, 0.01, 0.04)" 1 -- \
      2 0 2 0 0 3 x y z "10*(y-x)" "x*(28-z)-y" "x*y-8*z/3" 1.0 20 1 0 \
      14.999 15.001 14.999 15.001 35.999 36.001
check "Eg6 Roessler  H=1  eps=1.0" \
      "(-1.74, 1.86, 0.03) +- (0.16, 0.18, 0.00)" 1 -- \
      2 0 2 0 0 3 x y z -y-z "x+0.2*y" "0.2+z*(x-5.7)" 1.0 20 1 0 \
      0.9 1.1 1.9 2.1 2.9 3.1

echo
echo "the other refinement strategies, StepB variants and tube degrees"
check "Eg2 H=1 eps=0.5  method=1"      "(-2.13, 0.58) +- (0.20, 0.09)"  5  -- \
      2 0 1 2 0 2 x y y "(1-x^2)*y-x" 0.5 20 1 0  -3.1 -2.9 2.9 3.1
check "Eg2 H=1 eps=0.5  stepB=0"       "(-2.13, 0.58) +- (0.21, 0.10)" 5 -- \
      2 0 0 0 0 2 x y y "(1-x^2)*y-x" 0.5 20 1 0  -3.1 -2.9 2.9 3.1
check "Eg2 H=1 eps=0.5  stepB=3"       "(-2.13, 0.58) +- (0.20, 0.09)" 5 -- \
      2 0 0 3 0 2 x y y "(1-x^2)*y-x" 0.5 20 1 0  -3.1 -2.9 2.9 3.1
check "Eg2 H=1 eps=0.5  tubedegree=1"  "(-2.13, 0.58) +- (0.20, 0.09)" 5 -- \
      2 0 0 2 1 2 x y y "(1-x^2)*y-x" 0.5 20 1 0  -3.1 -2.9 2.9 3.1

echo
echo "Boundary mode"
check "Eg1 boundary 2D  H=1"           "(0.08, 1.46) +- (0.01, 0.07)" 4 -- \
      2 1 0 2 0 2 x y 2*x-2*x*y -y+x*y 1.0 20 1 0  0.9 1.1 2.9 3.1
check "Eg6 boundary 3D  H=1"           "(-1.74, 1.86, 0.03) +- (0.15, 0.18, 0.00)" 6 -- \
      2 1 0 2 0 3 x y z -y-z "x+0.2*y" "0.2+z*(x-5.7)" 1.0 20 1 0 \
      0.9 1.1 1.9 2.1 2.9 3.1

if [ "$BASELINE" = "1" ]; then exit 0; fi
echo
echo "$pass passed, $fail failed"
[ "$fail" -eq 0 ]

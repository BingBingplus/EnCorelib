#!/usr/bin/env bash
# ===========================================================================
#  src/repro/audit_cover.sh  --  box-by-box check of the computed covers
#
#  For each problem this runs EndCorelib at output level 3 and hands E0.txt
#  and E1.txt to ./EndCoreaudit, which checks that every box of the cover C
#  meets the validated images of the centre and the corners of its own
#  sub-box of B0 (see src/tools/audit_cover.cpp).  src/repro/verify_cover.sh checks
#  that the sample images lie in some box of C; this check also catches a
#  box of C that does not enclose the image of its sub-box.
#
#  Usage
#      bash src/repro/audit_cover.sh                all problems, methods 0 1 2
#      bash src/repro/audit_cover.sh eg1 eij
#      EPS=0.05 METHODS="0 2" STEPB=3 bash src/repro/audit_cover.sh eg2
#
#  Requires:  make             (./EndCorelib)
#             make EndCoreaudit (./EndCoreaudit)
# ===========================================================================
set -u

cd "$(dirname "$0")/../.." || exit 1
BIN=./EndCorelib;    [ -x "$BIN" ]  || BIN=./EndCorelib.exe
[ -x "$BIN" ]       || BIN=./bin/EndCorelib.exe
AUD=./EndCoreaudit;  [ -x "$AUD" ]  || AUD=./EndCoreaudit.exe
[ -x "$AUD" ]       || AUD=./bin/EndCoreaudit.exe
for f in "$BIN" "$AUD"; do
  [ -x "$f" ] || { echo "build first:  make && make EndCoreaudit" >&2; exit 1; }
done

EPS=${EPS:-0.1}
METHODS=${METHODS:-"0 1 2"}
STEPB=${STEPB:-0}
OUT=repro/out/audit
mkdir -p "$OUT"

#      tag | n vars rhs | box | H      (right-hand sides multiplied out)
eij="eij|2 x y 0 x^2|-1 1 0 0|1"
eg1="eg1|2 x y 2*x-2*x*y -y+x*y|0.9 1.1 2.9 3.1|4"
eg2="eg2|2 x y y -x+y-x^2*y|-3.1 -2.9 2.9 3.1|1"
eg3="eg3|2 x y y x^2|0.95 1.05 -1.05 -0.95|1"
eg4="eg4|2 x y 0.5+x-y+(-1/3)*x^3 0.056+0.08*x-0.064*y|0.9 1.1 -0.1 0.1|1"
eg5="eg5|3 x y z -10*x+10*y 28*x-y-x*z (-8/3)*z+x*y|14.999 15.001 14.999 15.001 35.999 36.001|1"
eg6="eg6|3 x y z -y-z x+0.2*y 0.2-5.7*z+x*z|0.9 1.1 1.9 2.1 2.9 3.1|1"

WANT=${*:-"eij eg1 eg2 eg3 eg4 eg5 eg6"}
rc=0
TOP=$PWD
echo "eps=$EPS stepB=$STEPB order=20"

for key in $WANT; do
  spec=$(eval "printf '%s' \"\${$key:-}\"")
  [ -n "$spec" ] || { echo "unknown problem: $key" >&2; continue; }
  IFS='|' read -r tag prob box H <<EOF
$spec
EOF
  for m in $METHODS; do
    d=$OUT/${tag}_m$m
    mkdir -p "$d"
    ( set -f; cd "$d" && "$TOP/$BIN" 3 0 "$m" "$STEPB" 0 $prob "$EPS" 20 "$H" 0 $box \
        > run.txt 2>&1 )
    printf '%-4s H=%-3s method=%d  ' "$tag" "$H" "$m"
    ( set -f; "$AUD" $prob 20 "$H" "$d/E0.txt" "$d/E1.txt" ) || rc=1
  done
done
exit $rc

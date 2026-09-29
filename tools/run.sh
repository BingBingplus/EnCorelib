#!/usr/bin/env bash
# ===========================================================================
#  tools/run.sh  --  invoke EnCorelib from the Makefile
#
#  The Makefile keeps the variable names and the right-hand sides as single
#  comma-separated strings (`var`, `ff`).  Expanding them straight into a
#  recipe would let the shell interpret characters such as '(' in
#  "x*(28-z)-y", so the splitting is done here instead, with proper quoting.
#
#  Usage
#      bash tools/run.sh BIN output_mode mode method stepB tubedegree n \
#                        VAR_CSV FF_CSV eps order T debug  lo1 hi1 ...
#
#  A right-hand side may not itself contain a comma (commas separate the
#  f_i); write `2*x-2*x*y` rather than a function call with two arguments.
# ===========================================================================
set -u

if [ "$#" -lt 13 ]; then
  echo "tools/run.sh: too few arguments" >&2
  exit 2
fi

BIN=$1;          shift
OUTPUT_MODE=$1;  shift
MODE=$1;         shift
METHOD=$1;       shift
STEPB=$1;        shift
TUBEDEGREE=$1;   shift
N=$1;            shift
VAR_CSV=$1;      shift
FF_CSV=$1;       shift
EPS=$1;          shift
ORDER=$1;        shift
T=$1;            shift
DEBUG=$1;        shift
# whatever is left is the initial box: lo1 hi1 lo2 hi2 ...

split_csv() {           # split_csv <csv> -> fills the global array SPLIT
  local csv=$1 old=$IFS
  SPLIT=()
  IFS=','
  for part in $csv; do SPLIT+=("$part"); done
  IFS=$old
}

split_csv "$VAR_CSV"; VARS=("${SPLIT[@]}")
split_csv "$FF_CSV";  FFS=("${SPLIT[@]}")

if [ "${#VARS[@]}" -ne "$N" ] || [ "${#FFS[@]}" -ne "$N" ]; then
  echo "tools/run.sh: n=$N but got ${#VARS[@]} variables and ${#FFS[@]} functions" >&2
  exit 2
fi

if [ "${SHOW_ONLY:-0}" = "1" ]; then
  printf '%s ' "$BIN" "$OUTPUT_MODE" "$MODE" "$METHOD" "$STEPB" "$TUBEDEGREE" "$N"
  printf "'%s' " "${VARS[@]}" "${FFS[@]}"
  printf '%s ' "$EPS" "$ORDER" "$T" "$DEBUG" "$@"
  printf '\n'
  exit 0
fi

exec "$BIN" "$OUTPUT_MODE" "$MODE" "$METHOD" "$STEPB" "$TUBEDEGREE" "$N" \
     "${VARS[@]}" "${FFS[@]}" "$EPS" "$ORDER" "$T" "$DEBUG" "$@"

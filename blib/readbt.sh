#!/bin/bash
base=$(mktemp XXXXX)
rm -f $base
mkdir $base
csplit -s -z -f $base/xx - '%===== BACKTRACE START ==========%1' '/===== BACKTRACE END ==========/'
sed '/\[.*\]/ {s/.*\[\(.*\)\].*/\1/}' ${base}/xx00 | addr2line -C -e $1 -f -i
rm -rf ${base}

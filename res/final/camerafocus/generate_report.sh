#! /bin/sh
SCRIPTDIR="$(dirname $0)"
./generate_report.py $SCRIPTDIR/cleanedup/ --calibration=$SCRIPTDIR/../../../calibrations/20161219172120574/ --camera=7 --extract-num

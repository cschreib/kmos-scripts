#!/bin/bash

CALIBS="01 02 03 04"
RAW_DIR="/data2/slow_data/programs/kmos-cluster"
GRATING="KKK"

for CAL in ${CALIBS}; do
    mkdir -p calib-$CAL
    cd calib-$CAL
    ../reduce calib $RAW_DIR"/calib-$CAL/" grating=$GRATING
    cd ..
done


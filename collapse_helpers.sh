#!/bin/bash

GRATING="KKK"
HELPERS="[129994]"
SCIS="01 02 03 04 05 06"

for SCI in ${SCIS}; do
    mkdir -p sci-$SCI/helpers
    cd sci-$SCI/helpers
    ../../reduce helpers ../ grating=$GRATING helpers=$HELPERS
    ./reduce.sh
    cd ../../
done

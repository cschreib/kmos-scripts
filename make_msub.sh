#!/bin/bash

SCIS="01 02 03 04 05 06 07 08 09 10"
for SCI in $SCIS; do
    mkdir -p sci-msub-$SCI
    cd sci-msub-$SCI/
    ../median_sub ../sci-$SCI
    cd ../
done

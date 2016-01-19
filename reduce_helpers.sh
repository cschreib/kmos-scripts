#!/bin/bash

DIRS="01 02 03 04 05 06"
for i in $DIRS; do
    cd sci-$i/helpers/
    ./reduce.sh
    cd ../../
done

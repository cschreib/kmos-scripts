#!/bin/bash

RAW_DIR="/data2/slow_data/programs/kmos-cluster"
GRATING="KKK"

reduce_wrapper() {
    mkdir -p stdstar-$1
    cd stdstar-$1
    ../reduce stdstar $RAW_DIR/calib-std-$1/ $2 grating=$GRATING
    cd ..
}

reduce_wrapper 01 calib=[../calib-01,../calib-03] # 2015-12-26T06:50:53.933
reduce_wrapper 02 calib=[../calib-01,../calib-03] # 2015-12-26T08:29:13.527
reduce_wrapper 03 calib=[../calib-02,../calib-03] # 2015-12-27T06:14:04.526
reduce_wrapper 04 calib=[../calib-04,../calib-03] # 2016-01-14T06:00:57.680
reduce_wrapper 05 calib=[../calib-04,../calib-03] # 2016-01-14T08:29:05.593

#!/bin/bash

RAW_DIR="/data2/slow_data/programs/kmos-cluster"
GRATING="KKK"

reduce_wrapper() {
    mkdir -p $1
    cd $1
    ../reduce sci $RAW_DIR/$1/ $2 $3 grating=$GRATING
    cd ..
}

reduce_wrapper sci-01 calib=[../calib-01,../calib-03] stdstar=../stdstar-01 # 2015-12-26T06:04:36.909
reduce_wrapper sci-02 calib=[../calib-01,../calib-03] stdstar=../stdstar-02 # 2015-12-26T07:46:54.425
reduce_wrapper sci-03 calib=[../calib-02,../calib-03] stdstar=../stdstar-03 # 2015-12-27T06:53:24.551
reduce_wrapper sci-04 calib=[../calib-04,../calib-03] stdstar=../stdstar-04 # 2016-01-14T05:23:19.716
reduce_wrapper sci-05 calib=[../calib-04,../calib-03] stdstar=../stdstar-04 # 2016-01-14T06:41:46.130
reduce_wrapper sci-06 calib=[../calib-04,../calib-03] stdstar=../stdstar-05 # 2016-01-14T07:41:30.778


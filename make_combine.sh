#!/bin/bash

reduce_wrapper() {
    mkdir -p $1-master
    cd $1-master
    ../reduce combine ../$1
    cd ..
}

reduce_wrapper sci



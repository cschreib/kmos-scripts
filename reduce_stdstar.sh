DIRS="01 02 03 04 05"
for i in $DIRS; do
    cd stdstar-$i/
    ./reduce.sh
    cd ..
done

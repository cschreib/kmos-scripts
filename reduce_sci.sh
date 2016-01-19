DIRS="01 02 03 04 05 06"
for i in $DIRS; do
    cd sci-$i/
    ./reduce.sh
    cd ..
done

DIRS="01 02 03 04"
for i in $DIRS; do
	cd calib-$i/
	./reduce.sh
	cd ..
done

#!/bin/sh

# Clean all and create directories for output
rm -rf *.dat *~

for (( k=1; k<5; k++ ))
do
	if [ $k -lt 10 ];
	then
		rm -rf result/test0$k
		mkdir result/test0$k
	else
		rm -rf result/test$k
		mkdir result/test$k
	fi	
done

# Benchmark case
./examples/image_segmentation/img_seg test/01_test_CV.dat
mv result/*.dat result/test01

# Adaptivity test
./examples/image_segmentation/img_seg test/02_test_LE.dat
mv result/*.dat result/test02

# Noisy image
./examples/image_segmentation/img_seg test/03_test_noise.dat
mv result/*.dat result/test03

# Pattern image
./examples/image_segmentation/img_seg test/04_test_tigers.dat
mv result/*.dat result/test04

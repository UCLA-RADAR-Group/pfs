#!/bin/bash
#run as : ./validate.sh

fft_param_1=0
down_param_1=0

# test tone data

# Test 1: fft

pfs_fft -m 32 -r 2.98023223876953125 -n 1 -f 3.125 -s 1000,11000 -b -o result.fftb test_tone.bin 

# note that you can get an ASCII version of this FFT by executing
# pfs_fft -m 32 -r 2.98023223876953125 -n 1 -f 3.125 -s 1000,11000 -o result.fft test_tone.bin 
# the result shows a strong narrow spike at +1 MHz

od -f correct.fftb | head -100 > correct.tmp1
od -f result.fftb  | head -100 > result.tmp1

cat correct.tmp1 | awk '{print $2}' > correct_hd.cmp
cat result.tmp1 | awk '{print $2}' > result_hd.cmp

paste correct_hd.cmp result_hd.cmp > fft.cmp

cat fft.cmp | awk '{if(($1-$2)/$1>0.0001 || ($1-$2)/$1<-0.0001) print ($1-$2)}' > err

if [ $(du -k err | cut -f1) -eq "0" ];then # test passed because comparison shows no meaningful differences
    fft_param_1=1; else fft_param_1=0;
fi

if [ $(du -k result_hd.cmp | cut -f1) -eq "0" ];then # test failed because output file size is 0
    fft_param_1=0;
fi

# Test 2: downsampling 

pfs_downsample -m 32 -d 1000 -s 1000 -I 0.0354 -Q -0.0148 -o result.d1000 test_tone.bin 

od -f correct.d1000 | head -100 > correct.tmp2
od -f result.d1000  | head -100 > result.tmp2

cat correct.tmp2 | awk '{print $2}' > correct_hd.cmp
cat result.tmp2 | awk '{print $2}' > result_hd.cmp

paste correct_hd.cmp result_hd.cmp > dwnsmpl.cmp

cat dwnsmpl.cmp | awk '{if(($1-$2)/$1>0.0001 || ($1-$2)/$1<-0.0001) print ($1-$2)}' > err

if [ $(du -k err | cut -f1) -eq "0" ];then # test passed because comparison shows no meaningful differences
    down_param_1=1; else down_param_1=0;
fi


#=====================================================


echo "================================================="
echo "=                    RESULTS                    ="
echo "================================================="
if [ $fft_param_1 -eq 1 ]; then echo " FFT test PASSED "; fi
if [ $down_param_1 -eq 1 ]; then echo " Downsampling test PASSED "; fi

if [ $fft_param_1 -eq 0 ]; then echo " FFT test FAILED "; fi
if [ $down_param_1 -eq 0 ]; then echo " Downsampling test FAILED "; fi

# clean up 
rm *tmp1 *tmp2 *cmp err
rm result*

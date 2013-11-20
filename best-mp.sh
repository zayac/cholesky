#!/bin/bash
#
# This script runs a series of benchmarks for S-Net Cholesky decomposition.
# Change the number of cores, block sizes, and matrix dimensions below.
#

#
# A list of CPU core numbers according to MACHINE-SPECIFIC!!!
# information from /sys/devices/system/cpu/...
# This is a list of core numbers such that cores are part of the
# same L3 cache and processor socket as much as possible.
#
#CPUS=0,4,8,12,16,20,1,5,9,13,17,21,2,6,10,14,18,22,3,7,11,15,19,23,24,28,32,36,40,44,25,29,33,37,41,45,26,30,34,38,42,46,27,31,35,39,43,47
#CPUS=0,4,1,5,2,6,3,7
CPUS=0,1,2,3,4,5,6,7,8,9,10,11
#CPUS=0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47
#
# Execute a single run of the benchmarked program.
#
function run () {
    # number of CPUs
    N=$1
    # executable
    X=$2
    # matrix dimension
    D=$3
    # block size
    B=$4
    # input file
    I="a-$D"
    # CPU set
    C=`echo "$CPUS" | grep -Eo "^([,]?[0-9]+){$N}"`
    # array file
    A=a-$D

    if [ ! -f $A ]; then
        echo "### Generating matrix data for $D ..." >&2
        ./symposdef $D -o $A || exit 1
    fi
    if [ ! -f $I ]; then
        echo "### Generating input image $D $B ..." >&2
        ./gen $D $B $A /dev/null > $I || exit 1
    fi

    log=/dev/shm/$USER-$$-$RANDOM.log
    err=/dev/shm/$USER-$$-$RANDOM.err
    trap 'rm -f -- $log $err; exit 1' 1 2 15
    echo "### $X $D $B -i $I " >&2
    $X $D $B -i $I > $log 2>$err
    grep -v -e 'Warning: Destroying' $err
    elapsed=`grep '^Time for' $log | grep -Eo '[0-9]+\.[-+e0-9]+'`
    thr=${X##*/}
    thr=${thr##-*}
    printf "%s %d %d %d %.3f\n" $thr $N $D $B $elapsed

    rm -f $log $err
    trap 1 2 15
}

#
# Loop over all of the needed parameters.
#
for cpu in 12
do
    for exe in ./openmp
    do
        for dim in 4096
        do
            for bs in 128 256 64 512 32 1024 2048 4096
            do
#                if [ $bs -ge 4096 ]
#                then
                    run $cpu $exe $dim $bs
#                fi
            done
        done
    done
done


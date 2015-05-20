#!/bin/bash

# sudo sh check_snipper_raid_disks.sh
# Check the health of the snipper raid disk partitions.

# Check the health of megaraid disks 0-1 on /dev/sda
# smartctl -a -d megaraid,0 /dev/sda | grep 'SMART Health'
# smartctl -a -d megaraid,1 /dev/sda | grep 'SMART Health'
for i in {0..1}
do 
        echo "smartctl -a -d megaraid,$i /dev/sda | grep 'SMART Health'"
        smartctl -a -d megaraid,$i /dev/sda | grep 'SMART Health' 
done

# Check the health of megaraid disks 2-7 on /dev/sdb
# smartctl -a -d sat+megaraid,2 /dev/sdb | grep 'SMART overall'
# smartctl -a -d sat+megaraid,3 /dev/sdb | grep 'SMART overall'
# smartctl -a -d sat+megaraid,4 /dev/sdb | grep 'SMART overall'
# smartctl -a -d sat+megaraid,5 /dev/sdb | grep 'SMART overall'
# smartctl -a -d sat+megaraid,6 /dev/sdb | grep 'SMART overall'
# smartctl -a -d sat+megaraid,7 /dev/sdb | grep 'SMART overall'
# smartctl -a -d sat+megaraid,8 /dev/sdb | grep 'SMART overall'
for i in {2..7}
do 
        echo "smartctl -a -d sat+megaraid,$i /dev/sdb | grep 'SMART overall-health'"
        smartctl -a -d sat+megaraid,$i /dev/sdb | grep 'SMART overall-health'
done

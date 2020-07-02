#!/bin/bash -l

config=$1
numa_bal=`cat $config | grep NUMA_BAL | awk '{print $3}'`
thp=`cat $config | grep THP | awk '{print $3}'`
nps=`cat $config | grep NPS | awk '{print $3}'`

numa_res=$(cat /proc/sys/kernel/numa_balancing)
if [[ "$numa_res" != "$numa_bal" ]]; then
    echo "Caution : NUMA balancing is $numa_res"
    exit 1
fi

thp_res=$(cat /sys/kernel/mm/transparent_hugepage/enabled | grep -Eo "\[.*\]" | sed -e "s@\[@@g" |  sed -e "s@\]@@g")
if [[ "$thp_res" != "$thp" ]]; then
    echo "Caution : THP is $thp_res"
    exit 1
fi

likwid-topology > topology.txt
nnuma=$(grep "NUMA domains:" topology.txt | cut -d":" -f 2)
nsockets=$(grep "Sockets:" topology.txt | cut -d":" -f 2)
NPS_float=$(echo "$nnuma/$nsockets" | bc -l)
NPS=${NPS_float%.*}
if (( $(echo "$NPS != $nps" |bc -l) )); then
    echo "Caution : NPS is $NPS, not $nps"
    exit 1
fi

rm topology.txt

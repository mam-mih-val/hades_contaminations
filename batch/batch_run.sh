#!/bin/bash

format='+%Y/%m/%d-%H:%M:%S'

date $format

job_num=$(($SLURM_ARRAY_TASK_ID))

filelist=$lists_dir/$job_num.list

cd $output_dir
mkdir -p $job_num
cd $job_num

while read line; do
    echo $line >> list.txt
done < $filelist
echo >> list.txt

source /etc/profile.d/modules.sh
module use /cvmfs/it.gsi.de/modulefiles/
module load compiler/gcc/9

echo "loading " $ownroot
source $ownroot

echo "executing $build_dir/analyse -i list.txt -o AuAu@1.23AGeV_PT2.root -p 2"
$build_dir/analyse -i list.txt -o AuAu@1.23AGeV_PT2.root -p 2

echo "executing $build_dir/analyse -i list.txt -o AuAu@1.23AGeV_PT2.root -p 2 -s"
$build_dir/analyse -i list.txt -o AuX@1.23AGeV_PT2.root -p 2 -s

echo "executing $build_dir/analyse -i list.txt -o AuAu@1.23AGeV_PT3.root -p 3"
$build_dir/analyse -i list.txt -o AuAu@1.23AGeV_PT3.root -p 3

echo "executing $build_dir/analyse -i list.txt -o AuX@1.23AGeV_PT3.root -p 3 -s"
$build_dir/analyse -i list.txt -o AuX@1.23AGeV_PT3.root -p 3 -s

echo JOB FINISHED!
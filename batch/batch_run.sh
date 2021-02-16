#!/bin/bash

format='+%Y/%m/%d-%H:%M:%S'

date $format

job_num=$(($SLURM_ARRAY_TASK_ID))

filelist=$lists_dir/$(ls $lists_dir | sed "${job_num}q;d")

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

echo "executing $build_dir/analyse -i list.txt -o output.root p 2"
$build_dir/analyse -i list.txt -o au_au.root -p 2 -e /lustre/nyx/hades/user/mmamaev/hades_contaminations/efficiency/efficiency_protons_auau123.root

echo "executing $build_dir/analyse -i list.txt -o output.root p 2"
$build_dir/analyse -i list.txt -o au_x.root -p 2 -e /lustre/nyx/hades/user/mmamaev/hades_contaminations/efficiency/efficiency_protons_auau123.root -s

echo JOB FINISHED!
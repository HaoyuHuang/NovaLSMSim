#!/bin/bash
home_dir="/proj/bg-PG0/haoyu"
script_dir="$home_dir/scripts"
results="/tmp/results"
exp_results_dir="$home_dir/nova-tinyrange-sim-java"
nova_dir="$home_dir/nova"
dryrun="false"

rm -rf $results
mkdir -p $results
mkdir -p $exp_results_dir

max_jobs=8
current_jobs=0
dbsize="10000000"
iterations="100000000"
enable_minor="true"

current_jobs=$(ps aux | grep nova_subrange_sim | grep java | grep -cv grep)
echo $current_jobs

function run_bench() {
	log_name="nova-$dist-$num_memtables-$subranges-$sampling-$enable_minor-$seed"
	echo "!!!!!! Running experiment $log_name"
	cmd="java -ea -jar $nova_dir/nova_subrange_sim.jar $dist $num_memtables $dbsize $subranges $iterations $sampling $enable_minor $seed"
	echo $cmd

	if [[ $dryrun == "true" ]]; then
		return
	fi
	nohup $cmd >& $results/$log_name &

	# Wait for all jobs to complete.
	while [ "$current_jobs" -gt "$max_jobs" ]
	do
		echo "Waiting jobs to complete. Number of running jobs: $current_jobs"
		sleep 10
		current_jobs=$(ps aux | grep nova_subrange_sim | grep java | grep -cv grep)
		echo "Waiting jobs to complete. Number of running jobs: $current_jobs"
		cp $results/nova-* $exp_results_dir/
	done
}

# cp $nova_dir/zipfian /tmp/zipfian

for seed in "0" "1" "2" "3" "4"
do
for enable_minor in "true" #"false"
do
for sampling in "1" #"1"
do
for dist in "zipfian" #"uniform"
do
for num_memtables in "256" "128" "4" "12" "64"
do
subranges=$((num_memtables/2))
current_jobs=$((current_jobs+1))
run_bench
done
done
done
done
done

# Wait for all jobs to complete.
while [ $current_jobs -gt 0 ]
do
  echo "Waiting jobs to complete. Number of running jobs: $current_jobs"
  sleep 10
  current_jobs=$(ps aux | grep nova_subrange_sim | grep java | grep -cv grep)
  echo "Waiting jobs to complete. Number of running jobs: $current_jobs"
  cp $results/nova-* $exp_results_dir/
done

# Save files
cp $results/nova-* $exp_results_dir/
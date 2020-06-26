# NovaLSMSimulator
NovaLSMSimulator to evaluate subranges. 

# Build
* nova_compute_ref_count.jar: Generated from https://github.com/HaoyuHuang/NovaLSMSim/blob/master/tst/edu/usc/nova/ComputeRefCounts.java
* nova_subrange_sim.jar: Generated from https://github.com/HaoyuHuang/NovaLSMSim/blob/master/db/edu/usc/nova/TinyRangeSim.java
# Run
Generate a file that stores the references count of each key with 10,000,000 keys database. The file is written in /tmp/zipfian. This is used as a reference the compute load imbalance. 
```bash
java -jar nova_compute_ref_count.jar 10000000
```
Run the following to reproduce the results shown in the paper.
```bash
bash scripts/nova_subrange_sim.sh
```



#!/bin/bash
#
# ****************************************
# Missingness Simulation Study
#
# BluePebble Automation Script
# This bash scripts automates updating to the most recent
# code version from GitHub and submitting the job to BP
#
# Emma Tarmey
#
# Started:          06/10/2025
# Most Recent Edit: 10/11/2025
# ****************************************

echo ""

# delete older version
rm -f -r Missingness-Sim-Study

# clone most recent version
git clone https://github.com/RaspberryEmma/Missingness-Sim-Study

# change wd
cd Missingness-Sim-Study
cd bash

# import python
module load languages/python/3.12.3



# submit simulation to BP HPC
for i in 1 2;
do
	echo   "Submitting job: launch_BP_run_"$i".sh"
	sbatch "launch_BP_run_"$i".sh"
done


# check jobs submitted correctly
sleep 5.0
echo ""
sacct -X
echo ""


#!/bin/bash
#SBATCH --mail-user=indrani.nayak@nationwidechildrens.org
#SBATCH --mail-type=FAIL
#SBATCH --cpus-per-task=1
#SBATCH --job-name="bootstapped_gamma"
#SBATCH --time=1-0
set -e 
ml purge
ml load bionetgen
ml load Python/3.8
export PYTHONPATH=~/.local/lib/python3.8/site-packages/:$PYTHONPATH
#python ca_signal_fuction.py analysis5
for((i=0; i<50; i++)); do 
    #python ca_signal_fuction.py analysis${i}
    # Ali suggestion
    sbatch --time=2-0 --cpus-per-task=1 --wrap="python ca_signal_function.py analysis${i}"

done
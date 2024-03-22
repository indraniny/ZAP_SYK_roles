#!/bin/bash
#SBATCH --mail-user=Indrani.Nayak@nationwidechildrens.org
#SBATCH --mail-type=FAIL
#SBATCH --cpus-per-task=1
#SBATCH --job-name="hCD3z_CD3z_R2"
#SBATCH --time=5-00:00:00
set -e 
ml purge
#ml load bionetgen
ml load Python/3.8
export PYTHONPATH=~/.local/lib/python3.8/site-packages/:$PYTHONPATH
python ca_signal_fuction.py
#python3 -m pdb ca_signal_fuction.py


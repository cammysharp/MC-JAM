#!/bin/bash
#
#SBATCH --job-name=specbin
#SBATCH --output=out_%j.txt
#
#SBATCH --ntasks=1
#SBATCH --time=60:00
#SBATCH --mem=20G

ml python-scientific/3.10.4-foss-2022a

source /fred/oz059/cammy/JAM/bin/activate
python MC_JAM_singlethread.py
deactivate


#!/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=1-0
sibeliaz -k 15 -b 300 -t 32 -a 2240 -n -o sibeliaz_out for_sibeliaz_2_contigs.fna

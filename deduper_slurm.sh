#!/bin/bash
#SBATCH -A bgmp                ###
#SBATCH --partition=bgmp       ### Partition
#SBATCH --job-name=Deduper        ### Job Name
#SBATCH --output=Deduper_»j.out ### file to store output
#SBATCH --error=Deduper-&j.err ### File to store job error messages
#SBATCH --time=12:00:00        ### WallTime
#SBATCH --nodes=1              ### Number of Nodes
#SBATCH --ntasks=1             ### Number of tasks per array job
#SBATCH --cpus-per-task=8            ### Number of CPU cores per task

#Pathways to files:
conda activate base 

# -f "/projects/bgmp/tonib/bioinfo/Bi624/Dedupe/deduper/Deduper/C1_SE_uniqAlign.sam"
# -o = "/projects/bgmp/tonib/bioinfo/Bi624/Dedupe/deduper/Deduper/C1_SE_uniqAlign_sbatched.sam"
# -u = "/projects/bgmp/tonib/bioinfo/Bi624/Dedupe/deduper/Deduper/STL96.txt"

./brooks_deduper.py -f /projects/bgmp/tonib/bioinfo/Bi624/Dedupe/deduper/Deduper/C1_SE_uniqAlign.sam -u /projects/bgmp/tonib/bioinfo/Bi624/Dedupe/deduper/Deduper/STL96.txt -o C1_SE_uniqAlign_sbatched.sam

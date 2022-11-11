#!/bin/bash
#SBATCH -A bgmp                ###
#SBATCH --partition=bgmp       ### Partition
#SBATCH --job-name=Deduper        ### Job Name
#SBATCH --output=Deduper_»j.out ### file to store output
#SBATCH --error=Deduper-&j.err ### File in which to store job error messages
#SBATCH --time=12:00:00        ### WallTime
#SBATCH --nodes=1              ### Number of Nodes
#SBATCH --ntasks=1             ### Number of tasks per array job
#SBATCH --cpus-per-task=8            ### Number of CPU cores per task, same as saying 8 cores

#Pathways to files:
conda activate base 

# file= "/projects/bgmp/tonib/bioinfo/Bi624/Dedupe/deduper/Deduper/C1_SE_uniqAlign.sam"
# out_file= "/projects/bgmp/tonib/bioinfo/Bi624/Dedupe/deduper/Deduper/C1_SE_uniqAlign_sbatched.sam"
# umi= "/projects/bgmp/tonib/bioinfo/Bi624/Dedupe/deduper/Deduper/STL96.txt"

# ./projects/bgmp/tonib/bioinfo/Bi624/Dedupe/deduper/scripts/brooks_deduper.py -f $file -u $umi -o $out_file

./brooks_deduper.py -f /projects/bgmp/tonib/bioinfo/Bi624/Dedupe/deduper/Deduper/C1_SE_uniqAlign.sam -u /projects/bgmp/tonib/bioinfo/Bi624/Dedupe/deduper/Deduper/STL96.txt -o C1_SE_uniqAlign_sbatched.sam
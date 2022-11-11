#!/usr/bin/env python3.10

import argparse
import re
import math

def get_args():
    parser = argparse.ArgumentParser(description="A program to check and remove any duplicate reads after aligning to a reference genome")
    #parser.add_argument("-qs", "--quality-score", help="quality score cutoff value", type=int)
    parser.add_argument("-f", "--file", help="designates absolute file path to sorted sam file", required=True, type=str)
    parser.add_argument("-o", "--outfile", help="designates absolute file path to sorted sam file", required=False, type=str)
    parser.add_argument("-u", "--umi", help="designates file containing the list of UMIs", required=True, type=str)
    # parser.add_argument("-h", "--help", help="prints a USEFUL help message and any assumptions your code makes", required=False)
    return parser.parse_args()
args=get_args()

#output=open(args.outfile, 'w')

#empty dictionaries
known_umi_dict = {}
unique_record= {}

#set my counters
unknown_umi_count = 0
header_count = 0
unique_record_count = 0
duplicate_record_count =0

#to run script: ./brooks_deduper.py -u /projects/bgmp/tonib/bioinfo/Bi624/Dedupe/deduper/Deduper/STL96.txt -f /projects/bgmp/tonib/bioinfo/Bi624/Dedupe/deduper/Deduper/test.sam
#to run script: ./brooks_deduper.py -u /projects/bgmp/tonib/bioinfo/Bi624/Dedupe/deduper/Deduper/STL96.txt -f /projects/bgmp/tonib/bioinfo/Bi624/Dedupe/deduper/Deduper/test_case -o out_deduper.sam

# umi_file= /projects/bgmp/tonib/bioinfo/Bi624/Dedupe/deduper/Deduper/STL96.txt
# sam_file= /projects/bgmp/tonib/bioinfo/Bi624/Dedupe/deduper/Deduper/test.sam

# opening known umi file and storing umis as keys in a dictionary 
with open(args.umi, 'r') as umi_file:
    for line in umi_file:
        umi= line.strip() #strip new line character
        known_umi_dict[umi]=()
#        print(known_umi_dict)

#for cigar string operators
leftmost_S = 0
rightmost_S = 0
int_M = 0
int_N = 0
int_D = 0


#opening sam file and specifiying wanted header positions
with open(args.file, 'r') as sam_file, open(args.outfile, 'w') as out_sam:
#with open(args.file, 'r') as sam_file:
    for line in sam_file: #read lines of sam file
        if line.startswith("@"): #if the line does not start with @, continue
            header_count += 1
            #print(header_count)
        else:
            header= line.split("\t") #split header by tabs
            umi_sam= header[0].split(":")[-1]  #position of umi dictionary value
            #print(umi_sam)

            if umi_sam not in known_umi_dict:
                unknown_umi_count += 1
                continue

            chromosome= header[2] #position of chromosome
            position= int(header[3]) #mapping position
            flag= int(header[1])
            cigar= header[5]
        
            if ((flag & 16) == 16):
                direction = 1 #if the direction is equal to 1 then this is the reverse direction
            else:
                direction = 0 #if the direction is equal to 0 then this is the forward direction
      
            tru_cigar= re.findall(r'(\d+)([A-Z]{1})', cigar) #figured out from help w/ Peter to create a list that seperates the number and letter in the cigar string
            # tru_cigar=re.findall(r'([0-9+])([MIDNSHPX=]+)', cigar) # initially tried to use this but it only returned the last digit in the list. 
            #print(tru_cigar)
            if direction == 0: #if we are on the forward strand
                for i in tru_cigar:
                    S= re.search(r'^\d+S', cigar)
                    if S in i:
                        S= S[0]
                        leftmost_S = int(S)
                tru_pos = position - leftmost_S 
                #print(tru_pos)
                 
            else: #we are on the reverse strand 
                for i in tru_cigar:
                    S= re.search(r'\d+S$', cigar)
                    if S in i:
                        S= S[0]
                        rightmost_S = int(S)

                    M = re.search(r'(\d+)[M]', cigar)               
                    if M in i:
                        M=M[0]
                        int_M = int(M)
                        
                    N = re.search(r'(\d+)[N]', cigar)
                    if N in i:
                        N=N[0]
                        int_N = int(N)
                        
                    D = re.search(r'(\d+)[D]', cigar)
                    if D in i:
                        D=D[0]
                        int_D = int(D)

                    tru_pos= (int(position)+ int_M + int_N + int_D + rightmost_S)       
                    #print("this is adjusted position:", tru_pos)
            
            unique= (umi_sam, chromosome, tru_pos, direction) 
            #print(unique)

            if unique not in unique_record:
                unique_record[unique] = 1
                unique_record_count +=1
                out_sam.write(line)
            else:
                duplicate_record_count += 1


print(f"Number of unique records:", unique_record_count)
print(f"Number of headers:", header_count)
print(f"Number of unknown umis:", unknown_umi_count)
print(f"Number of duplicates:", duplicate_record_count)


# #Final output:
# Number of unique records: 13937544
# Number of headers: 64
# Number of unknown umis: 0
# Number of duplicates: 4248866

#sort command
#cat C1_SE_uniqAlign_sbatched.sam | grep -v "@" | cut -f 3 | sort | uniq -c > unique_reads_per_chrom.txt
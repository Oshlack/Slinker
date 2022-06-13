import pandas as pd
import pysam as py
import sys, argparse

# Handle CLI

sam_path = sys.argv[1]
region = sys.argv[2]
filter = sys.argv[3]
out = sys.argv[4]

if not sam_path or not region:
    print("Please add both a sam and region")
    sys.exit()
  
# Same script in Canvas tool

mapq=255

multi = {}
unique = {}
    
sam = py.AlignmentFile(sam_path, "rb")

# Removed the loop here
offset = 1

""" Taken from Alex Dobin's excellent STAR aligner SJ script - reimplemented in Python."""
count = 0
readSJs = {}

for read in sam.fetch(str(region)): # Changed this to just need the superTranscript id.
    cigar = read.cigartuples
    t = 1
    g = read.reference_start + 1 + offset # as 0 based - SAM format is 1 based.

    if read.query_name not in readSJs.keys():
        readSJs[read.query_name] = {}

    for operation in cigar:
        op_type = operation[0]
        op_length = operation[1]
        if (op_type == 4) or (op_type == 1):
            t += op_length
        elif op_type == 2:
            g += op_length
        elif op_type == 3:
            sj1 = read.reference_name + ":" + str(g) + "-" + str(g + op_length - 1)
            try:
                readSJs[read.query_name][sj1] += 1
            except KeyError:
                readSJs[read.query_name][sj1] = 1

            if readSJs[read.query_name][sj1] == 1:
                if read.mapping_quality >= mapq:
                    unique[sj1] = unique.get(sj1, 0) + 1
                else:
                    multi[sj1] = multi.get(sj1, 0) + 1
            else:
                count += 1 

            g += op_length
        else: # M operation
            g += op_length
            t += op_length

sj = pd.DataFrame([unique, multi]).transpose().fillna(0)
sj.columns = ["Unique", "Multi"]
sj = sj.astype(int)

if filter:
    sj = sj[sj["Unique"] > int(filter)]

if out:
    print(sj)

print(count)

sj.to_csv("slinker_junctions.txt")
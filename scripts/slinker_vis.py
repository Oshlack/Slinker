#=======================================================================================================================
#
#   CASE vs. CONTROL TEST - Build a simple case vs. control plot with some annotation
#	Output is both a html file and a png of the resulting plot.
#
#   Author: Breon Schmidt
#   License: MIT
#
#=======================================================================================================================

''' --------------------------------------------------------------------------------------------------------------------
Imports
---------------------------------------------------------------------------------------------------------------------'''

import sys
import Slinker as sl

''' --------------------------------------------------------------------------------------------------------------------
R U N   T E S T
---------------------------------------------------------------------------------------------------------------------'''

# Setup
gene = sys.argv[1]
case_bam = sys.argv[2]
fmt = sys.argv[3]
resources = sys.argv[4]
output = sys.argv[5]
width = int(sys.argv[6])
junctions = int(sys.argv[7])
log = sys.argv[8]

if log == "true":
	log = True
else:
	log = False

# To get the case id we want to get all values on the left and the right of %, then get that string.
fmt_left = fmt.split("%")[0].split("*")[-1] # There should only be one %
fmt_right = fmt.split("%")[1].split("*")[-1]

case_id = case_bam.split(fmt_left)[1].split(fmt_right)[0]

print(case_id)

# Slinkerfy
slinker = sl.Slinker(gene=gene, case_id=case_id, resources=resources, padding=100, min_junctions=junctions)

# Plot
slinker.plot(log=log, cpm=False, min_junctions=junctions, width=width,
             title="<b>Slinker</b> - " + case_id + " (" + gene + ")",
             save=output+"/"+case_id+"_"+gene)


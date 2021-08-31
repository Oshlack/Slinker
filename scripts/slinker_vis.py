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
case_id = sys.argv[2]
resources = sys.argv[3]
output = sys.argv[4]
width = int(sys.argv[5])
junctions = int(sys.argv[6])
log = sys.argv[7]

if log == "true":
	log = True
else:
	log = False

case_id = case_id.split("_Aligned.sortedByCoord.out.bam")[0].split("/")[-1]

# Slinkerfy
slinker = sl.Slinker(gene=gene, case_id=case_id, resources=resources, padding=100)

# Plot
slinker.plot(log=log, cpm=False, min_junctions=junctions, width=width,
             title="<b>Slinker</b> - " + case_id + " (" + gene + ")",
             save=output+"/"+case_id+"_"+gene)


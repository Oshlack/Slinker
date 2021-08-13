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
genome = sys.argv[4]
output = sys.argv[5]

case_id = case_id.split("_Aligned.sortedByCoord.out.bam")[0].split("/")[-1]

# Slinkerfy
slinker = sl.Slinker(gene=gene, case_id=case_id, genome=genome, resources=resources, padding=100)

# Plot
slinker.plot(log=False, cpm=False, min_junctions=10, width=1000,
             title="<b>Slinker</b> - " + case_id + " (" + gene + ")",
             save=output+"/"+case_id+"_"+gene)


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

import Slinker as sl

''' --------------------------------------------------------------------------------------------------------------------
R U N   T E S T
---------------------------------------------------------------------------------------------------------------------'''

# Setup
gene = "ETV6"
case_id = "unhealthy"
resources = "resources"
genome = 19

# Slinkerfy
slinker = sl.Slinker(gene=gene, case_id=case_id, genome=genome, resources=resources, padding=100)

# Plot
slinker.plot(log=False, cpm=False, min_junctions=10, width=1000,
             title="<b>Slinker</b> - " + case_id + " (" + gene + ")",
             save=False)

# case_id+"_"+gene
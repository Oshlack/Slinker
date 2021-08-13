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

import Canvas as cv

''' --------------------------------------------------------------------------------------------------------------------
R U N   T E S T
---------------------------------------------------------------------------------------------------------------------'''

''' Set plot variables'''
#region = {"chr": 12, "start": 11802788, "end": 12048325} # Whole gene
region = {"chr": 12, "start": 11976000, "end": 11995000}
height = 1000
width = 1000

''' Set junctions variables'''
min_junctions = 10

''' Load the samples '''
bam_dir = "source"
samples = cv.load_samples(bam=bam_dir)

''' Then construct the plot layout. We simply need to create a numbered dictionary object. '''

layout = {}

layout[1] = {'type': 'axis', 'size': 1}

layout[2] = {'title': "Case",
			 'type': 'coverage',
			 'data': samples[0],
			 'size': 3,
			 'title_bgcolor': 'rgba(87, 22, 162, 1)',
			 'bgcolor': 'rgba(243, 232, 255, 1)',
			 'fill': 'rgba(137, 58, 228, 0.5)',
			 'line': 'rgba(87, 22, 162, 0.5)',
			 'log': False,
			 'cpm': False}

layout[3] = {'type': 'junctions',
			 'data': samples[0],
			 'size': 1,
			 'title_bgcolor': 'rgba(137, 58, 228, 1)',
			 'line': 'rgba(137, 58, 228, 0.5)',
			 'bgcolor': 'rgba(243, 232, 255, 1)',
			 'support': min_junctions}

layout[4] = {'type': 'gene',
			 'title_bgcolor': "rgba(255, 177, 51, 1)",
			 'bgcolor': "#fcf2d4",
			 'form': "gene",
			 'path': "etv6.gtf",
			 'title': "Transcripts",
			 'size': 4}

layout[5] = {'title': "Control",
			 'type': 'coverage',
			 'data': samples[1],
			 'size': 3,
			 'title_bgcolor': 'rgba(3, 181, 170, 1)',
			 'bgcolor': 'rgba(128, 161, 212, 1)',
			 'fill': 'rgba(23, 195, 178, 0.5)',
			 'line': 'rgba(3, 181, 170, 0.3)',
			 'log': False,
			 'cpm': False}

layout[6] = {'type': 'junctions',
			 'data': samples[1],
			 'size': 1,
			 'title_bgcolor': 'rgba(23, 195, 178, 1)',
			 'line': 'rgba(23, 195, 178, 1)',
			 'bgcolor': 'rgba(204, 252, 248, 1)',
			 'support': min_junctions}

highlights = [(11978578, 11979490, "rgba(249, 233, 0, 0.2)")]

''' The create the plot '''
test_plot = cv.Plot(layout, region, highlights=highlights,
					title="Example", height=height, width=width)










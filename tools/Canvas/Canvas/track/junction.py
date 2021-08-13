#=======================================================================================================================
#
#   CANVAS
#   Author: Breon Schmidt
#   License: MIT
#
#=======================================================================================================================

''' --------------------------------------------------------------------------------------------------------------------
Imports
---------------------------------------------------------------------------------------------------------------------'''

''' Internal '''


''' External '''
import plotly.graph_objects as go
import numpy as np
import pandas as pd

''' --------------------------------------------------------------------------------------------------------------------
Classes
---------------------------------------------------------------------------------------------------------------------'''

class Junctions():

	"""
	A class that represents a bam file and associated information.

	...

	Methods
	-------

	"""

	def __init__(self, sample):
		self.sample = sample

	def generate(self, region, mapq=255):

		multi = {}
		unique = {}

		chrom = region["chr"] if isinstance(region["chr"], list) else [region["chr"]]
		for i in range(0, len(chrom)):
			offset = region["offset"][i]

			""" Taken from Alex Dobin's excellent STAR aligner SJ script - reimplemented in Python."""
			for read in self.sample.sam.fetch(str(chrom[i]), start=region["start"], end=region["end"]):
				cigar = read.cigartuples
				t = 1
				g = read.reference_start + 1 + offset # as 0 based - SAM format is 1 based.

				for operation in cigar:
					op_type = operation[0]
					op_length = operation[1]
					if (op_type == 4) or (op_type == 1):
						t += op_length
					elif op_type == 4:
						g += op_length
					elif op_type == 3:
						sj1 = read.reference_name + ":" + str(g) + "-" + str(g + op_length - 1)
						if read.mapping_quality >= mapq:
							unique[sj1] = unique.get(sj1, 0) + 1
						else:
							multi[sj1] = multi.get(sj1, 0) + 1
						g += op_length
					else: # M operation
						g += op_length
						t += op_length

		self.sj = pd.DataFrame([unique, multi]).transpose().fillna(0)
		self.sj.columns = ["Unique", "Multi"]
		self.sj = self.sj.astype(int)

	def filter(self, min_support, by="Unique"):

		if by == "both":
			sj = self.sj[(self.sj["Unique"] > min_support) and (self.sj["Multi"] > min_support)]
		elif by == "Unique":
			sj = self.sj[self.sj["Unique"] > min_support]
		elif by == "Multi":
			sj = self.sj[self.sj["Multi"] > min_support]
		else:
			sj = self.sj

		return sj

	def get(self, region, min_support):
		self.generate(region)
		sj_filtered = self.filter(min_support)

		return sj_filtered

	def print(self, canvas, region, min_support=0, line="rgba(87, 22, 162, 0.5)",
                    bgcolor=False, row=1, col=1):

		splice_junctions = self.get(region, min_support)
		distances = [-15, 10, -30, 25]
		distance = 0

		x_coords = []
		y_coords = []

		for breakpoints, size in splice_junctions.iterrows():
			bp_info = breakpoints.split(":")[1].split("-")
			bp_start = int(bp_info[0])
			bp_end = int(bp_info[1])
			offset = -6

			if bp_end < region["start"] or bp_start > region["end"]:
				continue

			mid_point = (bp_end + bp_start)/2
			x_coords += [bp_start - 1, mid_point, bp_end + 1, bp_end + 1]
			y_coords += [0, offset, 0, None]
			support = str(size["Unique"])

			canvas.add_annotation(x=mid_point, y=offset, text=support, row=row, col=col, ay=distances[distance])
			distance = 0 if distance == len(distances) -1 else distance + 1

		canvas.append_trace(go.Scatter(x=x_coords, y=y_coords, line_shape='spline', mode='lines',
									   line=dict(width=2, smoothing=1, color=line), showlegend=False),
							row=row, col=col)

		canvas.update_yaxes(visible=False, row=row, col=col, tickcolor='white',
							tickfont=dict(color='white', size=12), range=[-10, 0], autorange=False)
		canvas.update_annotations(dict(xref="x"+str(col), yref="y"+str(row), showarrow=True,
									   arrowhead=7, ax=0, borderwidth=0, bgcolor="#000000", borderpad=0,
									font=dict(size=8, color="#ffffff")),
								  row=row, col=col)


		if bgcolor:
			canvas.add_shape(dict(type="rect", layer="below", x0=0, y0=-10, x1=1, y1=0, fillcolor=bgcolor),
							 line=dict(width=0), row=row, col=col)
			canvas.layout.shapes[-1]['xref']='paper'


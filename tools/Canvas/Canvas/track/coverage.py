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

''' External '''
import plotly.graph_objects as go
import numpy as np

''' --------------------------------------------------------------------------------------------------------------------
Classes
---------------------------------------------------------------------------------------------------------------------'''

class Coverage():

	"""
	A class that represents a bam file and associated information.

	...

	Methods
	-------

	"""

	def __init__(self, sample, cov_range=False):
		self.sample = sample
		self.cov_range = cov_range

	def _assign_block(self, block, coverage_array, region, reads, offset):

		"""Build up coverage by adding 1 across each, gapless block."""

		start = block[0] - region["start"] + offset
		end = block[1] - region["start"] + offset
		size = region["end"] - region["start"]

		start_in = (0 <= start <= size)
		end_in = (0 <= end <= size)

		if start_in and end_in:
			coverage_array[start: end] += 1
			reads += 1
		elif start_in:  # 5' region of read in region
			coverage_array[start: size] += 1
			reads += 1
		elif end_in:  # 3' region of read in region
			coverage_array[0: end] += 1
			reads += 1

		return coverage_array, reads


	def get(self, coverage_array, region, offset=False, log=False, cpm=False):

		reads = 0
		count = 0

		for read in self.sample.sam.fetch(str(region["chr"]), start=region["start"], end=region["end"]):
			blocks = read.get_blocks()
			count += 1
			for block in blocks:
				coverage_array, reads = self._assign_block(block, coverage_array, region, reads, offset)

		if reads == 0:
			print("Warning: No reads within supplied region for", self.sample.name)

		return coverage_array

	def print(self, canvas, region, log=False, cpm=True, row=1, col=1,
			  hover_template='<b>Coverage</b>: %{y}' + '<br/><b>Coord</b>: %{x}', coord_map=False,
			  line="rgba(87, 22, 162, 0.5)", fill="rgba(115, 29, 216, 0.5)"):

		''' Print a coverage track onto the supplied canvas plot. '''

		''' Get coverage histogram '''

		x = np.array(range(region["start"], region["end"]))
		y = np.zeros(region["end"]-region["start"])

		chrom = region["chr"] if isinstance(region["chr"], list) else [region["chr"]]
		for i in range(0, len(chrom)):
			y = self.get(y, {"chr": chrom[i],
							 "start": region["start"],
							 "end": region["end"]},
						 offset = region["offset"][i],
						 log=log, cpm=cpm)

		if cpm and not log:
			y = y / self.sample.lib_size * 1000000
		elif cpm and log:
			y = np.log2(y + 1) / self.sample.lib_size * 1000000
		elif log:
			y = np.log2(y + 1)

		self.cov_range = max(y)

		if not isinstance(coord_map, bool):
			alt_x = []
			for st_pos in x:
				alt_x.append(coord_map[coord_map == st_pos].index.values[0])
		else:
			alt_x = []

		''' Print Coverage '''
		canvas.append_trace(go.Scatter(x=x, y=y, text=alt_x, fill='tozeroy', showlegend=False,
									   hovertemplate=hover_template, fillcolor=fill, line_color=line, line_shape='hv'),
							row=row, col=col)

		canvas.update_yaxes(title_text=self.sample.name, row=row, col=col,
							title_font=dict(size=12, color='white'), tickcolor='white',
							tickfont=dict(color='white', size=12), range=(0, max(y)), autorange=False)

		canvas.update_xaxes(visible=False, row=row, col=col)



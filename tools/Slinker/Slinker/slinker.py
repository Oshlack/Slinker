#=======================================================================================================================
#
#   CANVAS
#   Author: Breon Schmidt
#   License: MIT
#
#=======================================================================================================================

''' ---------------------------pip----------------------------------------------------------------------------------------
Imports
---------------------------------------------------------------------------------------------------------------------'''

''' Internal '''
from Slinker.errors import error
import Slinker.build as build

''' External '''
import Canvas as cv

''' --------------------------------------------------------------------------------------------------------------------
Classes
---------------------------------------------------------------------------------------------------------------------'''

'''
Input into slinker vis (outside of plotting interface stuff)

gene id
case id
directory with bams
reference annotation
genome
build path
'''

class Slinker():

	def __init__(self, gene=False, case_id=False, resources="resources", padding=10):

		self.gene = gene
		self.case_id = case_id
		self.resources = resources
		self._error_check()

		'''Plot'''
		self.padding = padding

		'''Samples'''
		self.samples = cv.load_samples(resources)
		self.case, self.controls = self._split_samples()

		'''Assembly'''
		self.color_event = self._get_colors()
		self.assembly = build.Assembly(gene, resources, colors=self.color_event)
		skipped_exons = self.assembly.skipped_exon(self.samples, self.assembly.st.st_region, self.case_id, min_support=10)
		self.assembly.novel_regions["se"]["pos"] = skipped_exons


	def _get_colors(self):

		color_event = {
			"ne": 'rgba(244, 178, 102, 1)',
			"ee": 'rgba(216, 17, 89, 1)',
			"ri": 'rgba(249, 233, 0, 1)',
			"te": 'rgba(66, 107, 179, 1)',
			"se": 'rgba(89, 52, 79, 1)'
		}

		return color_event

	def _error_check(self):
		if not self.gene:
			error(100, "gene")
		elif not self.case_id:
			error(100, "case id")

	def _split_samples(self):

		'''From a list of samples, segment the nominated case samples out from the lot.'''

		case_sample = [sample for sample in self.samples if sample.name == self.case_id][0]
		control_samples = [sample for sample in self.samples if sample.name != self.case_id]

		return case_sample, control_samples

	def plot(self, height=None, width=1000, log=False, cpm=True,
			 min_junctions=False, title="Slinker Figure", save=False):

		'''Progressively build a Slinker layout using the CANVAS paradigm.
		   Note, dictionary id must reflect row no. Nested will overlap.'''

		if not height:
			self.assembly.gene.add_transcripts(self.assembly.assembly.table)
			no_transcripts = self.assembly.gene.count_transcripts(self.assembly.gene.region)
			height = 400 + no_transcripts*100 + len(self.controls)*300

		hover_template = '<b>Coverage</b>: %{y}'+'<br><b>ST Coord</b>: %{x}<br>'+ \
						 '<b>Gen Coord</b>: %{text}</b><extra></extra>'

		layout = {}

		layout[1] = {'type': 'axis', 'size': 100}

		layout[2] = {'title': self.case_id,
					 'type': 'coverage',
					 'data': self.case,
					 'size': 200,
					 'title_bgcolor': 'rgba(87, 22, 162, 1)',
					 'bgcolor': 'rgba(243, 232, 255, 1)',
					 'fill': 'rgba(137, 58, 228, 0.5)',
					 'line': 'rgba(87, 22, 162, 0.5)',
					 'coord_map': self.assembly.st.st_map,
					 'hover_template': hover_template,
					 'log': log,
					 'cpm': cpm}

		layout[3] = {'type': 'junctions',
					 'data': self.case,
					 'size': 100,
					 'title_bgcolor': 'rgba(137, 58, 228, 1)',
					 'line': 'rgba(137, 58, 228, 0.5)',
					 'bgcolor': 'rgba(243, 232, 255, 1)',
					 'support': min_junctions}

		layout[4] = {'type': 'gene',
					 'title_bgcolor': "rgba(247, 146, 86, 1)",
					 'bgcolor': "rgba(253, 236, 216, 1)",
					 'form': "gene",
					 'gtf': self.assembly.st.ref_st_gtf,
					 'title': "Block",
					 'size': 100 + 50}


		'''A decision has to be made here, is it big? Or not so big?'''
		layout[5] = {'type': 'gene',
					 'title_bgcolor': "rgba(255, 177, 51, 1)",
					 'bgcolor': "#fcf2d4",
					 'form': "gene",
					 'colors': self.assembly.novel_regions,
					 'gtf': self.assembly.st.st_gtf,
					 'title': "superTranscripts",
					 'size': no_transcripts*100 + 50}

		c = 6

		for control in self.controls:

			layout[c] = {'title': self.case_id,
						 'type': 'coverage',
						 'data': control,
						 'size': 200,
						 'title_bgcolor': 'rgba(3, 181, 170, 1)',
						 'bgcolor': 'rgba(128, 161, 212, 1)',
						 'fill': 'rgba(23, 195, 178, 0.5)',
						 'line': 'rgba(3, 181, 170, 0.3)',
						 'coord_map': self.assembly.st.st_map,
						 'hover_template': hover_template,
						 'log': log,
						 'cpm': cpm}

			layout[c + 1] = {'type': 'junctions',
						 'data': control,
						 'size': 100,
						 'title_bgcolor': 'rgba(23, 195, 178, 1)',
						 'line': 'rgba(23, 195, 178, 1)',
						 'bgcolor': 'rgba(204, 252, 248, 1)',
						 'support': min_junctions}

			c += 2

		'''Highlight regions of interest'''
		for event in self.assembly.novel_regions.keys():
			for position in self.assembly.novel_regions[event]["pos"]:
				color = self.color_event[event].split(",")
				color = color[0] + ", " + color[1] + ", " + color[2] + ", " + "0.2)"
				if "highlights" in locals():
					highlights += [(position.left, position.right, color)]
				else:
					highlights = [(position.left, position.right, color)]

		if "highlights" not in locals():
			highlights = False

		''' The create the plot '''
		plot = cv.Plot(layout, self.assembly.st.st_region, highlights=highlights, title=title,
					   height=height, width=width, padding=self.padding)

		if save:
			plot.canvas.write_html(save + ".html")
			plot.canvas.write_image(save + ".png")

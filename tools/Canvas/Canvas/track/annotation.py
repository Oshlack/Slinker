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
import re
import pandas as pd
import plotly.graph_objects as go

''' --------------------------------------------------------------------------------------------------------------------
Classes
---------------------------------------------------------------------------------------------------------------------'''


class Annotation():

	"""


	"""

	def __init__(self, path=False, gtf=False, form=None, colors=None):
		self.path = path
		self.gtf = gtf
		self.form = form
		self.colors = colors

	def get(self):

		if self.path.split(".")[-1] == "gtf":
			gtf = GTF(self.path)

		return gtf

	def _gene_format(self, region, annotation):

		genes = {}

		''' First, given the region we want to know which genes are within it'''
		filtered_genes = annotation.table[annotation.table["feature"] == "gene"]
		#filtered_genes = filtered_genes[filtered_genes["seqname"] == region["chr"]]

		''' Now, create new gene objects for our lucky contestants'''
		for index, gene in filtered_genes.iterrows():
			genes[gene["gene_name"]] = Gene(gene["gene_name"], gene["start"], gene["end"], gene["strand"])
			transcript_annotations = annotation.table[annotation.table["gene_name"] == gene["gene_name"]]
			genes[gene["gene_name"]].add_transcripts(transcript_annotations)

		return genes

	def process_for_print(self, gene, region, features=False, padding=100, height=100, arrow_steps=20, transcript_no=0):

		if not features:
			features = {
				"utr": {"x": [], "y": []},
				"colors": {
					"light": {"x": [], "y": [], "c": "rgba(99, 99, 99, 1)"},
					"dark": {"x": [], "y": [], "c": "rgba(55, 55, 55, 1)"}},
				"lines": {"x":[], "y": [], "arrows": {"x": [], "y": []}},
				"annotation" : {"features": {
									"text": [], "x": [], "y": []},
								"transcripts": {
									"text": [], "x": [], "y": []}}
			}

		for transcript in gene.transcripts.values():

			if not transcript.in_viewport(region):
				continue
			else:
				transcript_no += 1

			exons, exon_no = transcript.get_exons(region)
			utrs = transcript.get_utr(region)
			y0 = height*transcript_no - height*(2/5) + padding

			'''First assign exons to correct color arrays'''
			if len(exons) == 0 and len(utrs) == 0:  # Can't figure out a way to generate these lists without iterating
				continue

			for pos in range(0, len(exons)):
				event = "light" if exon_no.iloc[pos] % 2 == 0 else "dark"
				text = str(exon_no.iloc[pos])

				'''Need to perform a colour check'''
				if self.colors is not None:
					for n_event in self.colors.keys():
						for n_event_pos in self.colors[n_event]["pos"]:
							n_left = n_event_pos.left
							n_right = n_event_pos.right
							left = exons[pos].left
							right = exons[pos].right

							if n_left == left and n_right == right:
								event = n_event
								text = n_event.upper()
								if event not in features["colors"].keys():
									features["colors"][event] = {"x": [], "y": [], "c": self.colors[n_event]["c"]}

								break

				offset = height*(1.25/5)
				text_offset = offset if exon_no.iloc[pos] % 2 == 0 else -offset

				'''Is this exon partly an UTR'''
				utr_overlap = utrs.overlaps(exons[pos])
				utr = utrs[utr_overlap]

				if len(utr) > 0:
					if utr.right == exons[pos].right:
						exon_end = utr.left[0]
						exon_start = exons[pos].left
					elif utr.left == exons[pos].left:
						exon_start = utr.right[0]
						exon_end = exons[pos].right

					features["utr"]["x"] += [utr.left[0], utr.right[0], utr.right[0]]
					features["utr"]["y"] += [y0, y0, None]

				else:
					exon_end = exons[pos].right
					exon_start = exons[pos].left

				features["colors"][event]["x"] += [exon_start, exon_end, exon_end]
				features["colors"][event]["y"] += [y0, y0, None]
				features["annotation"]["features"]["text"] += [text]
				features["annotation"]["features"]["x"] += [(exon_start + exon_end) / 2]
				features["annotation"]["features"]["y"] += [y0 + text_offset]

			'''Now let's calculate some lines'''
			transcript_vp_start = transcript.start if transcript.start > region["start"] else region["start"]
			transcript_vp_end = region["end"] if transcript.end > region["end"] else transcript.end
			features["lines"]["x"] += [transcript_vp_start, transcript_vp_end, transcript_vp_end]
			features["lines"]["y"] += [y0, y0, None]

			'''Arrows are relative to each line, every X bases'''
			distance = region["end"] - region["start"]
			first = True
			for i in range(transcript_vp_start, transcript_vp_end, int(distance/arrow_steps)):
				if first:
					first = False
					continue

				direction_offset = i-distance*0.005 if transcript.strand == "+" else i+distance*0.005
				features["lines"]["arrows"]["x"] += [direction_offset, i, direction_offset, i]
				features["lines"]["arrows"]["y"] += [y0+3, y0, y0-3, None]

			'''And add relevant annotations'''
			name_offset = height*(2/5)
			mid_point = (transcript_vp_start + transcript_vp_end) / 2
			features["annotation"]["transcripts"]["text"] += [transcript.name]
			features["annotation"]["transcripts"]["x"] += [mid_point]
			features["annotation"]["transcripts"]["y"] += [y0 - name_offset]

		return features, transcript_no

	def print(self, canvas, region, row=1, padding=25, height=100, bgcolor="#fcf2d4", title=False):

		if self.path:
			annotation = self.get()
		else:
			annotation = self.gtf

		if self.form == "gene":
			formatted = self._gene_format(region, annotation)
		else:
			formatted = None

		''' Create dictionaries depicting the colours of the exons/features and the X/Y coordinates. 
			This is preferable to simply printing each exon as there is a slowdown with each track added to the plot'''
		transcript_no = 0
		for gene in formatted.values():

			features, transcript_no = self.process_for_print(gene, region, padding=padding, height=height,
															 transcript_no=transcript_no)

			# Add transcript names
			for no in range(0, len(features["annotation"]["transcripts"]["x"])):
				canvas.add_annotation(x=features["annotation"]["transcripts"]["x"][no], xref="x",
									  y=features["annotation"]["transcripts"]["y"][no], yref="y",
									  text=features["annotation"]["transcripts"]["text"][no],
									  showarrow=False, row=row, col=1,
									  font={"color": "#000000", "size": 11})

			# Transcript lines
			canvas.add_trace(go.Scatter(x=features["lines"]["x"], y=features["lines"]["y"],
										mode='lines', hoverinfo='skip', marker_color="#999999",
										line_width=1, line_dash="dot"),
							row=row, col=1)

			# Transcript direction
			canvas.add_trace(go.Scatter(x=features["lines"]["arrows"]["x"], y=features["lines"]["arrows"]["y"],
										mode='lines', hoverinfo='skip', marker_color="#333333", line_width=1),
							row=row, col=1)

			max_width = 15
			lw = height*(2/5) if height*(2/5) < max_width else max_width

			# Features
			for color in features["colors"].keys():
				canvas.add_trace(go.Scatter(x=features["colors"][color]["x"], y=features["colors"][color]["y"],
											 mode='lines', hoverinfo='skip', line_width=lw,
											 marker_color=features["colors"][color]["c"]),
									  row=row, col=1)

			# UTR
			canvas.add_trace(go.Scatter(x=features["utr"]["x"], y=features["utr"]["y"],
										 mode='lines', hoverinfo='skip', line_width=int(lw/2),
										 marker_color="#000000"),
								  row=row, col=1)

			# Feature names
			canvas.add_trace(go.Scatter(x=features["annotation"]["features"]["x"],
										y=features["annotation"]["features"]["y"],
										text=features["annotation"]["features"]["text"],
										 mode='text+markers', marker_color="rgba(00,00,00,0)", hoverinfo='skip',
										 textposition="middle center"),
								  row=row, col=1)


		# Plot specifics
		track_height = height*transcript_no + 2*padding
		canvas.update_traces(textfont_size=9)
		canvas.update_layout(showlegend=False)
		canvas.add_shape(dict(type="rect", layer="below", x0=0, y0=0, x1=1, y1=track_height,
								   fillcolor=bgcolor), line=dict(width=0), row=row, col=1)
		canvas.layout.shapes[-1]['xref'] = 'paper'
		canvas.update_xaxes(visible=False, row=row, col=1)
		canvas.update_yaxes(title_text=title, showticklabels=False, showgrid=False,
							row=row, col=1, range=[0, track_height], autorange=False)

class GTF:

	def __init__(self, path):
		self.path = path
		self.table = self.load()

	def _get_attribute(self, attributes, attribute):
		found = re.search(attribute + r' \"(.*?)\";', attributes)
		if found:
			return re.search(attribute + r' \"(.*?)\";', attributes).group(1)
		else:
			return "None"

	def load(self):

		table = pd.read_csv(self.path, sep="\t", header=None)
		table.rename(columns={table.columns[-1]: "attribute"}, inplace=True)

		# Format columns
		for attribute in ["gene_id", "transcript_id", "reference_id", "gene_name"]:
			attribute_formatted = table["attribute"].apply(self._get_attribute, attribute=attribute)
			table[attribute] = attribute_formatted

		table.drop("attribute", axis=1, inplace=True)
		table.columns = ["seqname", "source", "feature", "start", "end",
								  "score", "strand", "frame", "gene_id",
								  "transcript_id", "reference_id", "gene_name"]

		return table

class Gene:

	def __init__(self, name, start, end, strand="NA"):
		self.name = name
		self.start = start
		self.end = end
		self.strand = strand
		self.transcripts = {}
		self.no_transcripts = 0

	def format_id(self, name):
		fname = name if "(" not in name else name.split(" ")[0]
		return fname

	def add_transcripts(self, annotation):

		transcript_annotations = annotation[annotation["feature"] == "transcript"]

		for index, transcript_annotation in transcript_annotations.iterrows():

			transcript_id = transcript_annotation["transcript_id"]
			transcript_name = transcript_annotation["reference_id"]
			transcript_start = transcript_annotation["start"]
			transcript_end = transcript_annotation["end"]

			'''Break down transcripts into regions'''
			transcript_annotation = annotation[annotation["transcript_id"] == transcript_id]
			exons = transcript_annotation[transcript_annotation["feature"] == "exon"]
			cds = transcript_annotation[transcript_annotation["feature"] == "CDS"]
			utr = transcript_annotation[transcript_annotation["feature"].isin(["five_prime_utr", "three_prime_utr"])]

			exon_list = list(zip(exons["start"], exons["end"]))
			cds_list = list(zip(cds["start"], cds["end"]))
			utr_list = list(zip(utr["start"], utr["end"]))

			self.transcripts[transcript_id] = Transcript(transcript_id, transcript_start, transcript_end,
														 self.strand, exons=exon_list,  cds=cds_list,
														 utr=utr_list)

	def count_transcripts(self, region):

		no_transcripts = 0
		for transcript in self.transcripts.keys():
			if self.transcripts[transcript].in_viewport(region):
				no_transcripts += 1

		return no_transcripts


class Transcript:

	def __init__(self, name, start, end, strand, exons=False, cds=False, utr=False):
		self.name = name
		self.start = start
		self.end = end
		self.strand = strand
		self.viewport = False

		if exons:
			self.exons = pd.arrays.IntervalArray.from_tuples(exons)
			self.cds = pd.arrays.IntervalArray.from_tuples(cds)
			self.utr = pd.arrays.IntervalArray.from_tuples(utr)

	def in_viewport(self, region):

		viewport_exons = (region["start"] <= self.exons.left) & (self.exons.right <= region["end"])
		exons = self.exons[viewport_exons]

		if exons.shape[0] > 0:
			return True
		else:
			False

	def get_exons(self, region):
		viewport_exons = (region["start"] <= self.exons.left) & (self.exons.right <= region["end"])
		exon_range_max = len(self.exons) + 1
		exon_no = range(1, exon_range_max) if self.strand == "+" else range(exon_range_max - 1, 0, -1)
		exon_no = pd.Series(exon_no)[viewport_exons]
		exons = self.exons[viewport_exons]

		return exons, exon_no

	def get_utr(self, region):
		viewport_utr = (region["start"] <= self.utr.left) & (self.utr.right <= region["end"])
		utr = self.utr[viewport_utr]

		return utr
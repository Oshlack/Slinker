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
import Canvas as cv

''' External '''
import sys
from plotly.subplots import make_subplots

import time

''' --------------------------------------------------------------------------------------------------------------------
Classes
---------------------------------------------------------------------------------------------------------------------'''

class Plot():

	def __init__(self, layout, region, highlights=False,
				 height=1500, width=1000, bgcolor="#efefef",
				 title=False, padding=50,
				 coord_map=False):

		# User supplied order of tracks and subsequent details
		self.layout = layout
		self.highlights = highlights
		self._track_info()

		# Physical properties of the canvas
		self.height = height
		self.width = width
		self.bgcolor = bgcolor
		self.padding = padding

		# Define the canvas
		self.title = title

		# Define properties of the visualisation
		self.region = region
		self.cov_range = []

		# Create the visualisation
		self._create()


	def _create(self):

		self.blank()
		self.paint()

		if self.highlights:
			self.highlight()

		self.touchup()
		self.sign()
		self.canvas.update_layout(plot_bgcolor=self.bgcolor, margin=dict(l=50, r=50), autosize=False)
		self.canvas.show()


	def _track_info(self):

		# Get information about the layout information
		tracks = sorted(self.layout.keys())
		self.no_tracks = len(self.layout.keys())
		self.track_sizes = [self.layout[i]["size"] for i in tracks]
		self.track_types = [self.layout[i]["type"] for i in tracks]

	def blank(self):

		'''Generate a blank canvas'''

		self.canvas = make_subplots(rows=self.no_tracks, cols=1,
									vertical_spacing=0.0, shared_xaxes=True,
									horizontal_spacing=0.2,
									row_heights=self.track_sizes,
									start_cell='top-left')

		self.canvas.update_layout(height=self.height, width=self.width)


	def paint(self):

		'''Iterate through every track in the layout and plot it according to the type.'''

		for track_no, track in self.layout.items():

			''' ## NEED A TRACK VALIDATION BIT ## '''

			if track["type"] == "axis":
				axis = cv.Axis(self.canvas, self.region, row=track_no)

			elif track["type"] == "coverage":
				track_coverage = cv.Coverage(sample=track["data"])
				track_coverage.print(self.canvas, self.region,
									 log=track["log"], cpm=track["cpm"],
									 coord_map=track["coord_map"], hover_template=track["hover_template"], row=track_no,
									 line=track["line"], fill=track["fill"])
				self.cov_range.append(track_coverage.cov_range)

			elif track["type"] == "gene":
				form = track["form"] if track["form"] else None
				colors = track["colors"] if "colors" in track.keys() else None

				'''If not a path, then a gtf object'''
				if "path" in track.keys():
					track_annotation = cv.Annotation(path=track["path"], form=form, colors=colors)
				else:
					track_annotation = cv.Annotation(gtf=track["gtf"], form=form, colors=colors)

				track_annotation.print(self.canvas, self.region,
									   row=track_no, bgcolor=track["bgcolor"], title=track["title"])

			elif track["type"] == "junctions":
				track_junctions = cv.Junctions(sample=track["data"])
				track_junctions.print(self.canvas, self.region, min_support=track['support'],
									  row=track_no, line=track["line"], bgcolor=track["bgcolor"])

			else:
				print("The track type:", track["type"], "does not exist. Terminating.")
				sys.exit()

			title_bgcolor = track["title_bgcolor"] if "title_bgcolor" in track else "#FFFFFF"
			self._track_title(track_no, bgcolor=title_bgcolor)


	def sign(self, bg_color="#000000", font_size=16, font_color="#FFFFFF"):

		self.canvas.add_shape(dict(type="rect", yref="paper", xref="paper", layer="above",
								   x0=-1.1, y0=1, x1=1, y1=1.07, fillcolor=bg_color),
							  line=dict(width=0))

		self.canvas.add_annotation(x=0.5, y=1.035, text=self.title, ax=0, ay=0,
								   showarrow=True, arrowcolor="#ffffff", yref="paper", xref="paper",
								   font={"color": font_color, "size": font_size})

	def touchup(self):

		# Share ylimits
		self._scale_ylim()
		self.canvas.update_xaxes(visible=False, range=[self.region["start"] - self.padding,
													   self.region["end"] + self.padding])

	def _track_title(self, row, bgcolor="rgba(87, 22, 162, 1)"):

		y1 = 1 - sum(self.track_sizes[0:row - 1]) / sum(self.track_sizes)
		y0 = 1 - sum(self.track_sizes[0:row]) / sum(self.track_sizes)

		self.canvas.add_shape(dict(type="rect", yref="paper", xref="paper", layer="below",
								   x0=0, y0=y0, x1=-120, y1=y1, fillcolor=bgcolor),
							  line=dict(width=0))

	def _scale_ylim(self):

		max_y = max(self.cov_range)
		for i in range(0, len(self.track_types)):
			if self.track_types[i] == "coverage":
				j = i + 1
				self.canvas['layout']['yaxis' + str(j)].update(range=[0, max_y + 0.5])

	def highlight(self):

		if not self.highlights:
			return False

		for coords in self.highlights:
			start = coords[0]
			end = coords[1]
			color = coords[2]

			self.canvas.add_shape(dict(x0=start, y0=0, x1=end, y1=1, type="rect", yref="paper", fillcolor=color),
								  line=dict(width=0))



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


''' --------------------------------------------------------------------------------------------------------------------
Classes
---------------------------------------------------------------------------------------------------------------------'''

class Axis():

	def __init__(self, canvas, region, row=1, marker_no=10, height=4):
		self.region = region
		self.region_size = self.region["end"] - self.region["start"]
		self.canvas = canvas
		self.track_no = row
		self.marker_no = marker_no

		'''Arrow points'''
		self.height = height
		self.middle = round(self.height/2, 2)
		self.offset_top = self.middle + self.height*0.025
		self.offset_bottom = self.middle - self.height*0.025
		self.arrow_top = self.offset_top + self.height*0.025
		self.arrow_bottom = self.offset_bottom - self.height*0.025
		self.annotation_top = self.arrow_top + self.height*0.075
		self.annotation_bottom = self.arrow_bottom - self.height*0.05
		self.arrow_instep_end = self.region["end"] - self.region_size*0.01
		self.arrow_instep_start = self.region["start"] + self.region_size*0.01

		self.print_axis()

	def _axis_arrows(self):

		""" Create the 5' and 3' axis arrows."""

		self.canvas.add_shape({'type': "rect", 'x0': 0, 'x1': 1, "yref": "paper", "xref": "paper",
							   'y0': 1, 'y1': 1 - self.height, 'fillcolor': '#ffffff'},
							  line={'width': 0}, layer="below")

		self.canvas.add_shape({'type': "line",
							   'x0': self.region["start"], 'x1': self.region["end"],
							   'y0': self.offset_top, 'y1': self.offset_top,
							   'line': {'color': '#555555', 'width': 1}}, row=self.track_no, col=1)

		self.canvas.add_shape({'type': "line",
							   'x0': self.region["start"], 'x1': self.region["end"],
							   'y0': self.offset_bottom, 'y1': self.offset_bottom,
							   'line': {'color': '#555555', 'width': 1}}, row=self.track_no, col=1)

		self.canvas.add_shape({'type': "line",
							   'x0': self.region["end"], 'x1': self.arrow_instep_end,
							   'y0': self.offset_top, 'y1': self.arrow_top,
							   'line': {'color': '#555555', 'width': 1}}, row=self.track_no, col=1)

		self.canvas.add_shape({'type': "line",
							   'x0': self.region["start"], 'x1': self.arrow_instep_start,
							   'y0': self.offset_bottom, 'y1': self.arrow_bottom,
							   'line': {'color': '#555555', 'width': 1}}, row=self.track_no, col=1)

	def _axis_ticks(self):

		""" Begin making some ticks, find a multiple for the markers """
		marker_increment = int(self.region_size/self.marker_no)
		marker_step = self.region["start"] + marker_increment
		while marker_step < self.region["end"]:
			tick = "<b>" + str(marker_step) + "</b>"
			self.canvas.add_annotation(x=marker_step, y=self.offset_top, ay=-15, text=tick, ax=0,
									   showarrow=True, row=self.track_no, col=1, font={"color": "#000000", "size": 8},
									   arrowcolor="black", arrowwidth=1)
			marker_step += marker_increment
			if marker_step + marker_increment > self.region["end"]:
				break


		self.canvas.add_annotation(x=self.region["end"], y=self.annotation_top, text="3'", ax=0, ay=0,
								   showarrow=True, arrowcolor="#ffffff", row=self.track_no, col=1,
								   font={"color": "#000000", "size": 8})

		self.canvas.add_annotation(x=self.region["end"], y=self.annotation_bottom, text="5'", ax=0, ay=0,
								   showarrow=True, arrowcolor="#ffffff", row=self.track_no, col=1,
								   font={"color": "#000000", "size": 8})

		self.canvas.add_annotation(x=self.region["start"], y=self.annotation_top, text="5'", ax=-10, ay=0,
								   showarrow=True, arrowcolor="#ffffff", row=self.track_no, col=1,
								   font={"color": "#000000", "size": 8})

		self.canvas.add_annotation(x=self.region["start"], y=self.annotation_bottom, text="3'", ax=-10, ay=0,
								   showarrow=True, arrowcolor="#ffffff", row=self.track_no, col=1,
								   font={"color": "#000000", "size": 8})

	def print_axis(self):

		self._axis_arrows()
		self._axis_ticks()
		self.canvas.update_xaxes(visible=False, row=self.track_no, col=1)
		self.canvas.update_yaxes(showticklabels=False, showgrid=False,
								 row=self.track_no, col=1, range=[0, self.height], autorange=False)


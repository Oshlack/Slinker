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

import plotly
import sys

if int(plotly.__version__.split(".")[0]) < 4:
	print("Plotly Version 4 and above required - exiting.")
	sys.exit()

''' Internal '''

# Managing Data
from Canvas.data.sample import Sample
from Canvas.data.sample import load_samples

# Organising Data
from Canvas.track.axis import Axis
from Canvas.track.coverage import Coverage
from Canvas.track.junction import Junctions
from Canvas.track.annotation import Annotation

# Visualising Data
from Canvas.plot import Plot





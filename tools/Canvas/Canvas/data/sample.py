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
import pandas as pd
import numpy as np
import pysam as py
import plotly.graph_objects as go
from glob import glob

''' --------------------------------------------------------------------------------------------------------------------
Classes
---------------------------------------------------------------------------------------------------------------------'''

class Sample():

	"""
	A class that represents a bam file and associated information.

	...

	Methods
	-------

	"""

	def __init__(self, bam_path):
		self.name = bam_path.split("/")[-1].split("_Aligned.sortedByCoord.out.bam")[0]
		self.path = "/".join(bam_path.split("/")[:-1])
		self.bam_path = bam_path
		self._get_idxstats()
		self.lib_size = self.bam_stats["mapped"].sum() + self.bam_stats["unmapped"].sum()
		self.sam = py.AlignmentFile(bam_path, "rb")
		self.junctions = self._get_junctions()

	def _get_idxstats(self):
		colnames = ["chrom", "length", "mapped", "unmapped"]
		idxstats = py.idxstats(self.bam_path).rstrip()
		chroms = []

		for l in idxstats.split("\n"):
			chroms.append(l.split("\t"))

		self.bam_stats = pd.DataFrame.from_records(chroms, columns=colnames)
		for col in colnames[1:]:
			self.bam_stats[col] = pd.to_numeric(self.bam_stats[col])

	def _get_junctions(self):
		return cv.Junctions(sample=self.sam)



''' --------------------------------------------------------------------------------------------------------------------
Functions
---------------------------------------------------------------------------------------------------------------------'''

def load_samples(bam):

	bam_paths = bam if isinstance(bam, list) else glob(bam + "/*.bam")
	samples = [Sample(bam_path) for bam_path in bam_paths]

	return samples
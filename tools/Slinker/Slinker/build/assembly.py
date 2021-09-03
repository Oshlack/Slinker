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
from typing import List, Any, Tuple

''' Internal '''
from Slinker.build.supertranscript import ST

''' External '''
import pandas as pd
from Canvas.track.annotation import GTF
from Canvas.track.annotation import Gene

''' --------------------------------------------------------------------------------------------------------------------
Classes
---------------------------------------------------------------------------------------------------------------------'''

class Assembly():

	# Arguably you could have a superTranscript class

	def __init__(self, gene_name, dir, colors=False):
		self.gene_name = gene_name

		'''Load gtf files'''
		self.reference = GTFE(dir+"/"+self.gene_name+".gtf")
		self.assembly = GTFE(dir+"/assembly.combined.gtf", gene=self.gene_name, skip_row=2, assembly=True, filter_gene=True)
		self.gene = GeneE(self.gene_name, assembly=self.assembly.table)

		'''superTranscript Transformation'''
		self.st = ST(assembly=self.assembly, reference=self.reference, gene=self.gene)

		''' Find novel bits '''
		self.colors = colors
		self.novel_regions = self.find_novel()


	def find_gaps(self, assembly):

		first = True
		gaps = []
		padding = 100000

		for i, e in assembly.iterrows():

			if first:
				start = e["end"]
				first = False
				#gaps.append((start - padding, start))
				continue

			if abs(self.st.st_map[start] - self.st.st_map[e["start"]]) > 1:
				gaps.append((start, e["start"]))
			start = e["end"]

		#gaps.append((start, start + padding))

		return list(set(gaps))

	def find_novel(self):

		# Novel, truncated, intron, extended, AS
		novel_exons, ref_exons = self.assembly.get_exons(split=True)

		ref_exons_merge = ref_exons.sort_values("start")
		ref_exons_merge["group"] = (ref_exons_merge["start"] > ref_exons_merge["end"].shift().cummax()).cumsum()
		ref_exons_merge = ref_exons_merge.groupby("group").agg({"start": "min", "end": "max"})

		# At this point we only care about differences between the reference and the build
		novel_intervals = []
		ref_intervals: List[Tuple[Any, Any]] = []

		for i, e in novel_exons.iterrows():
			novel_intervals.append((e["start"], e["end"]))

		''' Account for strandedness in ref '''
		for i, e in ref_exons.iterrows():
			ref_intervals.append((e["start"], e["end"]))

		novel_differences = set(novel_intervals).difference(set(ref_intervals))

		''' Let's consider deletion events before we think about specifics'''
		deletions = self.has_deletion(novel_exons, ref_exons)

		''' We now have all the novel things, we need to determine what they might be '''
		novel_results = {}
		for novelty in ["ne", "ee", "ri", "te", "se"]:
			novel_results[novelty] = {"pos": [], "c": self.colors[novelty]}

		for novel in novel_differences:

			if self.is_deletion(novel, deletions):
				continue

			novel_adj = (novel[1], novel[0]) if self.st.st_map[novel[1]] < self.st.st_map[novel[0]] else novel

			if self.is_novel_exon(novel_adj, ref_exons_merge):
				novel_results["ne"]["pos"].append(self.is_novel_exon(novel_adj, ref_exons_merge))
			elif self.is_retained_intron(novel_adj, ref_exons_merge):
				novel_results["ri"]["pos"].append(self.is_retained_intron(novel_adj, ref_exons_merge))
			elif self.is_extended_exon(novel_adj, ref_exons_merge):
				novel_results["ee"]["pos"].append(self.is_extended_exon(novel_adj, ref_exons_merge))
			elif self.is_truncated_exon(novel_adj, ref_exons_merge):
				novel_results["te"]["pos"].append(self.is_truncated_exon(novel_adj, ref_exons_merge))
			else:
				continue

		return novel_results

	def _strand_correct(self, start, end):
		if start > end:
			temp = end
			end = start
			start = temp

		return start, end

	def is_novel_exon(self, novel, reference):

		gaps = self.find_gaps(reference)

		for gap in gaps:
			n_start, n_end = self._strand_correct(novel[0],  novel[1])
			start, end = self._strand_correct(self.st.st_map[novel[0]],  self.st.st_map[novel[1]])

			if (n_start > gap[0]) & (n_end < gap[1]):
				return pd.Interval(start, end)

		return False

	def skipped_exon(self, samples, region, case_id, min_support=10):

		# Do some storage
		case_junctions = []
		control_junctions = []

		for sample in samples:
			sample_junctions = sample.junctions.get(region, min_support, sam=sample.sam)

			for index, junction in sample_junctions.iterrows():
				coords = index.split(":")[-1].split("-")
				start = int(coords[0])
				end = int(coords[1])

				if sample.name == case_id:
					case_junctions.append((start, end))
				else:
					control_junctions.append((start, end))

		for junction in control_junctions:
			if junction in case_junctions:
				case_junctions.remove(junction)

		''' Now a skipped exon must be... at known exon boundaries! Let's have a look.'''

		exons = self.assembly.table[self.assembly.table["feature"] == "exon"]
		ref_assembly = exons[exons["gene_name"] != "None"]

		boundaries = {}
		for i, values in ref_assembly.iterrows():
			boundaries[self.st.st_map[values["start"]]] = True
			boundaries[self.st.st_map[values["end"]]] = True

		skipped_exons = []
		for junction in case_junctions:
			try:
				if boundaries[junction[0]] and boundaries[junction[1]]:
					skipped_exons.append(pd.Interval(junction[0], junction[1]))

			except KeyError:
				continue

		return skipped_exons

	def is_extended_exon(self, novel, reference):

		fixed_downstream = reference[reference["end"] == novel[0]]
		fixed_upstream = reference[reference["start"] == novel[1]]

		ee_downstream = fixed_downstream[fixed_downstream["end"] >= novel[1]]
		ee_upstream = fixed_upstream[novel[0] <= fixed_upstream["end"]]

		if (ee_downstream.shape[0] > 0 or ee_upstream.shape[0] > 0):
			start, end = self._strand_correct(self.st.st_map[novel[0]], self.st.st_map[novel[1]])
			return pd.Interval(start, end)

	def is_truncated_exon(self, novel, reference, distance=100):

		fixed_downstream = reference[reference["end"] == novel[1]]
		fixed_upstream = reference[reference["start"] == novel[0]]
		trunc_downstream = fixed_downstream[fixed_downstream["start"] < novel[0]]
		trunc_upstream = fixed_upstream[novel[1] < fixed_upstream["end"]]

		if (trunc_downstream.shape[0] > 0 or trunc_upstream.shape[0] > 0):
			start, end = self._strand_correct(self.st.st_map[novel[0]], self.st.st_map[novel[1]])
			return pd.Interval(start, end)

	def is_deletion(self, novel, deletions):

		for deletion in deletions:
			if novel[0] == deletion[0] or novel[0] == deletion[1]:
				return True
			elif novel[1] == deletion[0] or novel[1] == deletion[1]:
				return True

		return False

	def has_deletion(self, n_reference, reference):


		gaps = self.find_gaps(n_reference)
		no_transcripts = reference["transcript_id"].unique().shape[0]
		deletions = []

		'''For a deletion event to be a deletion event, these "gaps" cannot be in any reference transcripts.'''

		for gap in gaps:

			''' Either the gap is contained within the '''
			n_start, n_end = self._strand_correct(gap[0], gap[1])
			potential_del = reference[(reference["start"] < n_start) & (n_end < reference["end"])]

			if potential_del["transcript_id"].unique().shape[0] == no_transcripts:
				deletions.append((n_start, n_end))

		return deletions

	def is_retained_intron(self, novel, reference):

		gaps = self.find_gaps(reference)

		for gap in gaps:

			''' Either the gap is contained within the '''
			start, end = self._strand_correct(novel[0], novel[1])
			n_start, n_end = self._strand_correct(gap[0], gap[1])

			if (start <= n_start) & (end >= n_end):
				start, end = self._strand_correct(self.st.st_map[n_start], self.st.st_map[n_end])
				return pd.Interval(start, end)

		return False

class GeneE(Gene):

	def __init__(self, name, assembly):
		self.assembly = assembly
		self.region = self.get_region()
		strand = self.get_strand()
		super().__init__(name, self.region["start"], self.region["end"], strand)

	def get_strand(self):
		return self.assembly.iloc[0]["strand"]  # might have multiple genes

	def get_region(self):
		transcripts = self.assembly[self.assembly["feature"] == "transcript"]
		return {"chr": transcripts["seqname"].iloc[0],
				"start": transcripts["start"].min(),
				"end": transcripts["end"].max()}


class GTFE(GTF):

	def __init__(self, path, gene=False, skip_row=0, assembly=False, filter_gene=False):
		self.skip_row = skip_row
		self.gene = gene
		self.assembly = assembly
		self.filter_gene = filter_gene
		super().__init__(path)

		if self.gene:
			self.assembly_name = self.get_gene_ref()

	def load(self):

		'''Necessary for the blocks'''
		gtf = pd.read_csv(self.path, sep="\t", skiprows=self.skip_row, header=None)
		gtf.rename(columns={gtf.columns[-1]: "attribute"}, inplace=True)

		''' Format columns '''
		attributes = ["gene_id", "transcript_id", "reference_id"]
		if self.assembly:
			attributes += ["gene_name", "cov"]
		else:
			attributes += ["gene_name"]

		for attribute in attributes:
			if not self.assembly:
				if attribute == "reference_id":
					gtf["reference_id"] = " "
					continue

			attribute_formatted = gtf["attribute"].apply(self._get_attribute, attribute=attribute)
			gtf[attribute] = attribute_formatted

		gtf.drop("attribute", axis=1, inplace=True)
		columns = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame",
				   "string_ref", "transcript_id", "reference_id", "gene_name"]


		if self.assembly:
			columns.append("coverage")


		gtf.columns = columns
		if self.filter_gene:
			m_tag = gtf.loc[gtf["gene_name"] == self.gene, "string_ref"].values[0]
			gtf = gtf[(gtf["gene_name"] == "None") | (gtf["gene_name"] == self.gene)]
			gtf = gtf[gtf["string_ref"] == m_tag]

		if self.assembly:
			gtf = gtf.drop("coverage", axis=1)

		return gtf


	def get_gene_ref(self):
		assembly_name = list(self.table["string_ref"].unique())
		return assembly_name

	def get_exons(self, split=False):

		exons = self.table[self.table["feature"] == "exon"]

		if split:
			novel_exons = exons[exons["gene_name"] == "None"]
			ref_exons = exons[exons["gene_name"] != "None"]
			return novel_exons, ref_exons
		else:
			return exons

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

'''Internal'''

from Slinker.build.flatten import flatten_exons

'''External'''
import numpy as np
import pandas as pd

''' --------------------------------------------------------------------------------------------------------------------
Classes
---------------------------------------------------------------------------------------------------------------------'''

class ST():

	def __init__(self, assembly=False, reference=False, gene=False):

		'''Create the supertranscript from the assembly'''

		'''Setup'''
		self.gene = gene
		self.transcripts = assembly.table[assembly.table["feature"] == "transcript"]
		self.exons = self.get_adj_exons(assembly=assembly.table)
		self.exons_ref = self.get_adj_exons(assembly=reference.table)
		self.size = self.exons["end"].max() - self.exons["start"].min() + 1
		self.supertranscript, self.st_map = self.make_supertranscript()

		'''Now make ST coordinate transcripts'''
		self.st_gtf = GTF_ST(assembly.table, self.st_map, self.gene)
		self.ref_st_gtf = GTF_ST(reference.table, self.st_map, self.gene, block=True)
		self.st_region = self.get_st_region(assembly)

		'''And as every ST is forward stranded'''
		self.gene.strand = "+"
		self.st_gtf.table["strand"] = "+"
		self.ref_st_gtf.table["strand"] = "+"

	def make_supertranscript(self):
		flat_ref = self._flat_reference()
		if self.gene.strand == "+":
			st_map = pd.Series(range(1, flat_ref.shape[0] + 1), index=flat_ref.index, name="st_coord")
		else:
			st_map = pd.Series(range(flat_ref.shape[0], 0, -1), index=flat_ref.index, name="st_coord")

		st = self._flat_to_st(flat_ref)

		return st, st_map

	def _flat_reference(self):

		# 1 = In exon
		# 2 = Exon boundaries

		# Make a flat representation that represents the entire region
		exon_region = []

		all_exons = pd.concat([self.exons, self.exons_ref], join="inner")

		for index, exon in all_exons.iterrows():

			start = self.gene.region["start"] + exon["start"]
			end = self.gene.region["start"] + exon["end"] + 1
			flat = np.ones(end-start, dtype=int)
			flat[0] = 2
			flat[-1] = 2
			exon_region.append(pd.DataFrame(flat, index=range(start, end)))

		exon_regions = pd.concat(exon_region, axis=1)
		exon_regions.fillna(0)
		flat = exon_regions.max(axis=1).astype(int)

		return flat

	def get_adj_exons(self, assembly=False):

		exons = assembly[assembly["feature"] == "exon"].copy()
		exons["start"] = list(exons["start"] - self.gene.region["start"])
		exons["end"] = list(exons["end"] - self.gene.region["start"])

		return exons

	def _make_st_map(self, flat):

		st_map_build = self._flat_to_st(flat)
		st_map_build = pd.DataFrame(st_map_build)
		st_map_build.columns = ["type"]
		st_map_build["st_coord"] = list(range(0, st_map_build.shape[0]))
		st_map = st_map_build["st_coord"]

		return st_map

	def _flat_to_st(self, flat):

		st = flat.copy()
		st.index = range(1, flat.shape[0] + 1)
		prev = 0
		drop = []
		for pos, value in st.iteritems():
			if prev == 2:
				drop.append(pos)
			prev = value

		# Drop points inter-boundary
		st.loc[drop] = 1
		st = st[st != 1]

		return st

	def get_st_region(self, assembly):

		offset = [1] # We are saying that these coordinates are 1-based (important for Pysam).

		return {"chr": assembly.assembly_name,
				"start": self.supertranscript.index[0],
				"end": self.supertranscript.index[-1] + 1,
				"offset": offset}

class GTF_ST():

	def __init__(self, gtf, st_map, gene, block=False):
		self.block = block
		self.st_map = st_map
		self.gene = gene
		self.table = self.make_st_gtf(gtf)

	def make_st_gtf(self, gtf):

		if "coverage" in gtf.columns:
			total = gtf.loc[gtf["feature"] == "transcript", "coverage"].astype(float).sum()
			columns = gtf.columns.drop("coverage")
		else:
			columns =gtf.columns

		table = pd.DataFrame(columns=columns)  # perhaps should be numpy for speed

		for i, entry in gtf.iterrows():

			modify = entry.copy()
			st_start = self.st_map[int(entry["start"])]
			st_end = self.st_map[int(entry["end"])]
			gene_ref = entry["gene_name"] if entry["gene_name"] != "None" else entry["string_ref"]

			if "coverage" in gtf.columns:
				coverage = round(float(entry["coverage"]) / total * 100, 2)
				trans_id = entry["transcript_id"] + " (" + str(coverage) + "%)"
			else:
				trans_id = entry["transcript_id"]

			modify["transcript_id"] = trans_id
			modify["start"] = st_start
			modify["end"] = st_end
			modify["gene_name"] = gene_ref
			modify = pd.DataFrame(modify).transpose()
			table = pd.concat([table, modify], join="inner")

		if self.gene.strand == "-":
			starts = table["start"].copy()
			ends = table["end"].copy()
			table["start"] = ends
			table["end"] = starts


		'''Now add a gene entry if none exit'''

		if table[table["feature"] == "gene"].shape[0] == 0:
			genes = list(table["gene_name"].value_counts().index)
			first = True
			for gene in genes:
				contents = table[table["gene_name"] == gene]
				gene_line = pd.DataFrame(contents.iloc[0]).transpose()
				gene_line["feature"] = "gene"
				gene_line["start"] = int(contents["start"].min())
				gene_line["end"] = int(contents["end"].max())
				gene_line["reference_id"] = "None"

				if first:
					final_table = pd.concat([gene_line, contents], join="inner")
					first = False
				else:
					final_table = pd.concat([final_table, gene_line, contents], join="inner")

		else:
			final_table = table

		''' Now finally, sort the exons by start position '''
		if not self.block:

			exons = final_table[final_table["feature"] == "exon"]
			final_table = pd.DataFrame(final_table.groupby("transcript_id").apply(lambda x: x.sort_values("start")))
		else:
			final_table.reset_index(inplace=True)

		if self.block:

			exons = final_table[final_table["feature"] == "exon"]
			exons_flat = flatten_exons(exons)

			model_exon = pd.DataFrame(exons.iloc[0]).transpose()
			table_block = pd.DataFrame(final_table.iloc[:2], columns=final_table.columns)

			for i, exon in exons_flat.iterrows():
				new_block = model_exon.copy()
				new_block["start"] = int(exon["start"])
				new_block["end"] = int(exon["end"])
				table_block = pd.concat([table_block, new_block], join="inner")
				table_block["reference_id"] = table_block["gene_name"][0] + " block"

			final_table = table_block.drop_duplicates()
			final_table = final_table[final_table["start"] != final_table["end"]]

		else:

			ref_transcripts = final_table["reference_id"] != "None"
			old_name = final_table.loc[:, "transcript_id"].str.split(' ', expand=True)[0]
			#coverage = final_table.loc[:, "transcript_id"].str.split(' ',expand=True)[1]
			final_table.loc[:, "transcript_id"] = old_name
			final_table.loc[ref_transcripts, "transcript_id"] = final_table.loc[ref_transcripts, "reference_id"]
			new_name = final_table.loc[:, "transcript_id"] 
			final_table.loc[:, "reference_id"] = new_name


		return final_table



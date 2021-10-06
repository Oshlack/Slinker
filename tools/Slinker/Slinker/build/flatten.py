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
import sys, argparse
import pandas as pd
import csv

''' --------------------------------------------------------------------------------------------------------------------
Functions
---------------------------------------------------------------------------------------------------------------------'''

def user_input():

	''' Get arguments and options from CLI '''

	cli = argparse.ArgumentParser(description="Slinker CLI")
	cli.add_argument('-gtf', '-g',
					 required=True,
					 help=("""Path to gtf needing flattening.
					 		  e.g. /path/to/to_flatten.gtf"""))


	cli.add_argument('-output', '-o',
					 required=True,
					 help=("""Output file for the flattened gtf.
					 		  e.g. /path/to/flattened.gtf"""))

	return cli.parse_args()


def flatten_exons(exons):

	exons_merge = exons.sort_values("start")
	exons_merge = exons_merge.drop_duplicates(subset=['start', 'end'])
	exons_merge["group"] = (exons_merge["start"] > exons_merge["end"].shift().cummax()).cumsum()
	exons_flat = []

	for i, exon in exons_merge.groupby('group'):
		if exon.shape[0] > 1:

			current_end = exon["end"]
			exon.loc[(exon["start"].shift(-1) < exon["end"]) &
					 (exon["start"].shift(-1) != exon["start"]), "end"] = exon["start"].shift(-1)
			ends_less = (exon["start"].shift(-1) < current_end) & \
						(current_end > current_end.shift(1))
			# exon.loc[ends_less, "start"] = current_end.shift(1).loc[ends_less] + 1
			exons_flat.append(exon)
		else:
			exons_flat.append(exon)

	exons_flat = pd.concat(exons_flat)
	exons_flat.drop("group", axis=1, inplace=True)

	return exons_flat

def gtf_to_flat(gtf_path=False, skip_row=0):

	if not gtf_path:
		print("No gtf was supplied. Exiting.")
		sys.exit()

	gtf = pd.read_csv(gtf_path, sep="\t", skiprows=skip_row, header=None)
	gtf.columns = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"]
	gtf_flat = pd.DataFrame(columns=gtf.columns)

	exons = gtf[gtf["feature"] == "exon"]
	exons_flat = flatten_exons(exons)

	block_no = 1
	model_exon = pd.DataFrame(exons.iloc[0]).transpose()
	for i, exon in exons_flat.iterrows():
		new_block = model_exon.copy()
		gene_name = str(new_block["attributes"]).split('gene_id "')[-1].split('"; transcript_id')[0]
		transcript_id = "block"
		attributes = 'gene_id "' + gene_name + \
					 '"; transcript_id "' + transcript_id + \
					 '"; block_number "' + str(block_no) + '";'

		new_block["start"] = int(exon["start"])
		new_block["end"] = int(exon["end"])
		new_block["attributes"] = attributes
		gtf_flat = pd.concat([gtf_flat, new_block], join="inner")
		block_no += 1

	gtf_flat = gtf_flat.drop_duplicates()
	gtf_flat = gtf_flat[gtf_flat["start"] != gtf_flat["end"]]

	return gtf_flat

''' --------------------------------------------------------------------------------------------------------------------
Run
---------------------------------------------------------------------------------------------------------------------'''

def main():
	user_args = user_input()
	gtf_flat = gtf_to_flat(gtf_path=user_args.gtf, skip_row=2)
	gtf_flat.to_csv(user_args.output, sep="\t", header=False, index=False, quoting=csv.QUOTE_NONE)

if __name__ == "__main__":
	main()




#!/usr/bin/python

import os
import sys
import re
import argparse
import numpy as np
import pandas as pd 

"""
This script is utilized to filter a given MZID TSV file from MSGF+ to a user-defined FDR. Default 0.01 FDR. 


i.e.
cd /waterlooms/scripts

# Invoke the script
python filter_fdr_mzid.py -mzid /data/PXD005573/Fig1_MP-DIA-120min-120kMS1-10W30k-14dppp_MHRM_R01_Q1.tsv
"""

# Parse the Arguments
class ErrorHandlingArgparse(argparse.ArgumentParser):
	def error(self, message):
		sys.stderr.write('error: %s\n' % message)
		self.print_help(sys.stderr)
		sys.exit(1)


parser = ErrorHandlingArgparse()
parser.add_argument("-mzid",
	help="MZID TSV from MSGF+",
	required=True
	)
parser.add_argument("-fdr",
	help="FDR Level (Default: 0.01) for 1% to Filter Peptides to",
	default=0.01
	)
args = parser.parse_args()

# Sanity check for file existance.
if not os.path.exists(args.mzid):
	parser.error("The MZID file %s does not exist!" % args.mzid)



# Read in the MSGF+ MZID TSV
mzid_df = pd.read_csv(args.mzid, sep="\t")

# Compute the negative log of the EValue
# Evalue is already sorted from lowest to highest with (lower) being better
# tmp = -1 * np.log(mzid_df.loc[:,"EValue"])

# Grep if the string contains "DeBruijn" 
mzid_df["Decoy"] = mzid_df.loc[:,"Protein"].str.contains('DeBruijn') 
mzid_df["Real"] = ~mzid_df["Decoy"]

mzid_df["Decoy"] = mzid_df["Decoy"].astype(int)
mzid_df["Real"] = mzid_df["Real"].astype(int)


# Compute Cumulative sum for Decoys (FALSE PSM)
mzid_df["Decoy"] = mzid_df["Decoy"].cumsum()

# FDR = # False PSM / # PSM Above Threshold
mzid_df["FDR"] = mzid_df["Decoy"] / pd.Series(range(1,len(mzid_df)))
trimmed_df = mzid_df[mzid_df['FDR'].lt(0.01)]


base_name = os.path.basename(args.mzid)
path_name = os.path.dirname(args.mzid)

# Export Peptides above 1% FDR
trimmed_df.to_csv(path_name + "/" + "0.01_fdr_" + base_name, sep="\t")








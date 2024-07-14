from pyteomics import mgf, pepxml, mass
import os
from urllib.request import urlopen, Request

import re

from pyteomics import parser

parser.parse('PEPTIDE')



# TODO:
# Convert the identification into modX notation
# 



# pyteomics takes the modX notation
# 

var reg = "/d\+1/"

if re.search('\+', "SDFSDF+") is not None:
	print(True)


# https://msgfplus.github.io/msgfplus/examples/Mods.txt

# C + C2H3N1O1

# C: +57.021
# M: +15.995
# N-Term: +42.011




def seq_to_modX(peptide_sequence)
	modX = ""
	return modX



aa_comp = dict(mass.std_aa_comp)
aa_comp['camC'] = mass.composition({'C': 2, 'H': 3, 'N': 1, 'O': 1})


pseq = "PEPTcamCIDE"
from pyteomics import mass
mass.calculate_mass(sequence=pseq, ion_type='b', charge=1)





# Residues with N-Term
A
P
M 





from pyteomics import mass
aa_comp = dict(mass.std_aa_comp)
aa_comp['camC'] = mass.Composition('C2H3N1O1')	# Carbamidomethyl Cap
aa_comp['oM'] = mass.Composition('O1') 			# Oxidation ofMethionine
aa_comp['onM'] = mass.Composition('C2H2O') 		# N-Term Acetylation
mass.calculate_mass('camCPEPTID', aa_comp=aa_comp)




mass.fast_mass2('H-PEPTIDE-OH')



# Extract out set of numbers +[0-9]
# Replace with _
# For each AA, calculate mass, 
# If encounter _, add the modification to rest of value(s)
# Repeat for y ions except subtract 

import re 
import pandas as pd
from pyteomics import mass

mod_pattern = "\+\d+\.\d+" # Grep for format +75.93 etc...
mod_peptide = "+42.011PE+15.995PTI+57.021DE"



# Get the modifications and convert them to float values
mods = re.findall(mod_pattern, mod_peptide)
#mods = [float(re.sub("\+", "", mod)) for mod in mods]
mods = [float(mod) for mod in mods]


substituted_peptide = re.sub(mod_pattern, "_", mod_peptide)

# Split string into characters
# For each character in the set
# Create a data frame

y_ions = list()
b_ions = list()
b_ions_mod = list()
y_ions_mod = list()
b_pre = list()
y_pos = list()
mod_index = 0
# For each residue in the set
for residue_index in range(len(substituted_peptide)):
	residue = substituted_peptide[residue_index]
	if residue != "_":
		# Compute fragment ions for the N->C [i] substring
		# Compute fragment ions for [i] N->C
		pref = re.sub("_", "", substituted_peptide[:residue_index])
		posf = re.sub("_", "", substituted_peptide[residue_index:])
		if (len(pref) == 0 or len(posf) == 0):
			print("continue:")
			continue
		#print ("b" + str(len(pref)) + ": " + pref)
		#print ("y" + str(len(posf)) + ": " + posf)
		# Add masses to b-ion
		b_pre.append(pref)
		y_pos.append(posf)
		b_mass = mass.fast_mass2(pref, ion_type="b", charge=1)
		b_ions.append(b_mass)
		b_ions_mod.append(b_mass + sum(mods[:mod_index]))
		y_mass = mass.fast_mass2(posf, ion_type="y", charge=1)
		y_ions.append(y_mass)
		y_ions_mod.append(y_mass + sum(mods[mod_index:]))
	else:
		print (residue)
		mod_index += 1

print(y_ions)
print(b_ions)

[x1 - x2 for (x1, x2) in zip(b_ions_mod, b_ions)]

d = {'y':y_ions_mod, 'b':b_ions_mod, 'b_fragment':b_pre,'y_fragment':y_pos, 'unmod_y':y_ions, 'unmod_b':b_ions}
peptide_df = pd.DataFrame(d)
peptide_df["y_diff"] = peptide_df["y"]-peptide_df["unmod_y"]
peptide_df["b_diff"] = peptide_df["b"]-peptide_df["unmod_b"]







================================================================================


# Extract DIA Scan


# Extract DIA Spectra By Window
"""
inspired by 0_build_spectra_by_window.py

input: MGF file
output: pkl file representing scans by window
"""

# Need lxml
# pip install lxml

from pyteomics import mzxml

mzxml_file = "/data/PXD005573/Fig1_MP-DIA-120min120kMS1-22W30k-8dppp_MHRM_R01.mzXML"

all_rts = list()
ms1_rts = list()
ms2_rts = list()
mzxml_in = mzxml.MzXML(mzxml_file)

# Can get index this way
# mzxml_in.get_by_index(0)

for scan in mzxml_in:
	all_rts.append(scan['retentionTime'])
	#print(scan)


ms1_indicies = list()
ms1_rts = list()
mzxml_in = mzxml.MzXML(mzxml_file)
for scan_index, scan in enumerate(mzxml_in):
	if scan['msLevel'] == 1:
		ms1_rts.append(scan['retentionTime'])
		ms1_indicies.append(scan_index)



ms1_indicies = list()
ms1_rts = list()
scan_index = 0
for scan in mzxml_in:
	if scan['msLevel'] == 1:
		ms1_rts.append(scan['retentionTime'])
		ms1_indicies.append(scan_index)
	scan_index += 1


# Average Frame Length in minute(s) Last Retention Time / num_retention_frames
# 10 seconds per peptide elution
# 6 elutions per minute
# (10 / 60)


# Minutes / Num_Frames
# ~ 0.001691 minutes per frame
# Need to get the number of FRAMES in 10 seconds


avg_frame_len = all_rts[-1] / len(all_rts)
avg_frame = avg_frame_len / (10 / 60 )
avg_frame
# Minutes / Frame
# Minutes


# 100th RT converted to seconds
# 1st RT conveted to seconds
# 100th - 1st = second difference
all_rts[99] * 60 - all_rts[0] * 60


















from pyteomics import mgf
mgf_file = "/data/PXD005573/Fig1_MP-DIA-120min120kMS1-22W30k-8dppp_MHRM_R01_Q1.mgf"
mgf_in = mgf.IndexedMGF(mgf_file)
#mgf_in.get_by_index(202)
mgf_in.get_by_index(202)['params']
rt_in_min = float(mgf_in.get_by_index(202)['params']['rtinseconds']) / 60


binary_search_list(rts, rt_in_min)




# {'params': {'pepmass': (1102.5203, None), 'charge': [2], 'rtinseconds': 2985.7922, 'title': 'Fig1_MP-DIA-120min120kMS1-22W30k-8dppp_MHRM_R01_Q1.203.203.2'}, 'm/z array': array([ 203.06662,  206.09227,  208.10786,  213.09735,  225.1237 ,




# select the ms1 scan
# select the m/z range that the precursor would've been acquired in
import pandas as pd
unscored_psm_file = "data/PXD005573/0.01_fdr_Fig1_MP-DIA-120min-120kMS1-10W30k-14dppp_MHRM_R01_Q1.tsv"
psm_in = pd.read_csv(unscored_psm_file, sep="\t")
peptide_lengths = list()
for index, row in psm_in.iterrows():
	print(index)
	protein = row.get("Protein")
	scan_num = row.get("ScanNum")
	title = row.get("Title")
	peptide = row.get("Peptide")
	peptide = re.sub("\+\d+\.\d+", "", peptide)
	# Sub out the mods
	peptide_lengths.append(len(peptide))



mzxml_in.get_by_index(203)



width = 0
for scan_index in range(1, 23):
	print(str(scan_index) + "/" + "23")
	width += int(mzxml_in.get_by_index(scan_index)['precursorMz'][0]['windowWideness'])
	print(mzxml_in.get_by_index(scan_index)['precursorMz'])







# Build mz lookup array
# assume scan 0 is ms1
scan_index = 1
mz_lower_bound = 0
mz_windows_lowerbound = list()
mz_windows_lowerbound.append(mz_lower_bound)
while (1):
	mzxml_scan = mzxml_in.get_by_index(scan_index)
	# Break when encountering 2nd 
	if (mzxml_scan['msLevel'] == 1):
		break
	mz_lower_bound = mzxml_scan['precursorMz'][0]['precursorMz']
	mz_windows_lowerbound.append(mz_lower_bound)
	scan_index += 1
	print(mzxml_scan['precursorMz'][0])



# Grab the lower offset for the window
[tmp, ms2_offset] = binary_search_list(mz_windows_lowerbound, 400)
if ms2_offset != 0:
	ms2_offset -= 1







#####=====#####=====#####=====#####=====#####=====#####=====#####=====#####=====

import pandas as pd 

mzxml_in = mzxml.MzXML(mzxml_file)

def build_ms1_rt_to_index_df(mzxml_in):
	""" Builds a DataFrame consisting of MS1_INDEX and MS1_RT
	:return a DataFrame with columns [MS1_INDEX, MS1_RT]
	"""
	ms1_indicies = list()
	ms1_rts = list()
	for scan_index, scan in enumerate(mzxml_in):
		if scan['msLevel'] == 1:
			ms1_rts.append(scan['retentionTime'])
			ms1_indicies.append(scan_index)
	d = {'MS1_INDEX':ms1_indicies, 'MS1_RT':ms1_rts}
	ms1_df = pd.DataFrame(d)
	return ms1_df


ms1_df = build_ms1_rt_to_index_df(mzxml_in)


[tmp_rt, index] = binary_search_list(ms1_df['MS1_RT'], 40.5)

ms1_offset = ms1_df['MS1_INDEX'][index]




#####=====#####=====#####=====#####=====#####=====#####=====#####=====#####=====

residue = "A"
def one_hot_encode_residue(residue):
	residue_one_hot_encoding_dict = {'A':0, 'R':1, 'N':2, 'D':3, 'C':4, 'Q':5, 'E':6, 'G':7, 'H':8, 'I':9, 'L':10, 'K':11, 'M':12, 'F':13, 'P':14, 'S':15, 'T':16, 'W':17, 'Y':18, 'V':19}
	# Empty array of 20, zeroes
	encoding = [0]*20
	# Set the encoding
	encoding[residue_one_hot_encoding_dict[residue]] = 1
	return encoding





one_hot_encode_residue("L")


charge = 2
pdf = generate_fragment_ions("PEPTIDE", charge)

pdf['b_h2o_loss'] = pdf['b'] - (float(CommonMass.H2O.value) / charge)
pdf.loc[len(pdf)-1, 'b_h2o_loss'] = 0
pdf['b_nh3_loss'] = pdf['b'] - (float(CommonMass.NH3.value) / charge)
pdf.loc[len(pdf)-1, 'b_nh3_loss'] = 0

pdf['y_h2o_loss'] = pdf['y'] - (float(CommonMass.H2O.value) / charge)
pdf.loc[0, 'y_h2o_loss'] = 0
pdf['y_nh3_loss'] = pdf['y'] - (float(CommonMass.NH3.value) / charge)
pdf.loc[0, 'y_nh3_loss'] = 0
pdf





# Building a custom tensor for each eluted peptide
import tensorflow as tf
import numpy as np 


# one_hot_encode_residue("L")
row = one_hot_encode_residue("L")
row.append(0.1) #b_mono			intensity
row.append(0.2) #y_mono			intensity
row.append(0.3) #b+1			intensity
row.append(0.4) #y+1			intensity
row.append(0.5) #b-H2O			intensity
row.append(0.6) #b-NH3			intensity
row.append(0.7) #y-H2O			intensity
row.append(0.8) #y-NH3 			intensity
row.append(0.9) #precursor 		intensity
row.append(0.95)#precursor+1 	intensity
# Done row operations
# intensity is normalized to the maximum intensity within the scan



tmp_tensor = tf.concat([row, row], 0)
tf.reshape(tmp_tensor, shape=(2, 30))
tf.reshape(tmp_tensor, shape=(2, 2, 15))
tf.reshape(tmp_tensor, shape=(1, 1, 2, 2, 15))





# Given a set of fragment ions
# Need to extract these from raw data
def extract_from_raw_data(mzxml_in, fragment_ion_df, ppm_tolerance, ms1_rt_index_df, ms1_window_offset, ms2_window_offset, max_peptide_len, n_features, n_ms2_windows):
	max_peptide_len = 40
	n_features = 30
	n_scans = 100
	ppm_tolerance = 10
	max_scan_index = len(mzxml_in)
	# ms1_index = ms1_window_offset + i * n_ms2_windows
	# ms2_index = ms1_window_offset + ms2_window_offset + i * n_ms2_windows
	extracted_features = list()
	for i in range(n_scans):
		print(i)
		# Compute MS1 Index
		# Compute MS2 Index
		# If they are out of range, append dummy data, else:
		ms1_scan_index = ms1_window_offset  + i * n_ms2_windows
		ms2_scan_index = ms1_scan_index + ms2_window_offset
		if (ms1_scan_index >= max_scan_index or ms2_scan_index >= max_scan_index):
			# Append dummy data
			# Still expect to have num features * max_peptide_len for this scan
			padding = [-1] * n_features * max_peptide_len
			continue
		else:
			# Load MS1 and MS2 scan from RAW
			ms1_raw = mzxml_in.get_by_index(ms1_scan_index)
			ms2_raw = mzxml_in.get_by_index(ms2_scan_index)
			max_ms1_intensity = np.max(ms1_raw['m/z array'])
			max_ms2_intensity = np.max(ms2_raw['m/z array'])
		# For each residue, extract the value by m/z
		for residue_index in range(max_peptide_len):
			if residue_index < len(fragment_ion_df.index):
				# Extract Encoding of that residue
				aa_encoding = fragment_ion_df.loc[residue_index]['residue']
				
				# Compute the closest match using binary search
				b_mz = fragment_ion_df.loc[residue_index]['b']
				[value, b_index] = binary_search_list(ms2_raw['m/z array'], b_mz)
				# Check if the closest match is within the ppm tolerance
				# If Yes, lookup intensity & add to the feature(s)
				# If No, append -1 as dummy data
				b_int = ms2_raw['intensity array'] / max_ms2_intensity
			else:
				padding = [-1] * n_features
				extracted_features.extend(padding)
	# For charges 1 and 2:
	# Iterate from 0,99
	# If within the valid length of the mzxml_in
	#	For each m/z extract the intensity
	# Else
	# Append dummy -1 data
	return extracted_features








#row.append(pdf.loc[0, 'b'])
#row.append(pdf.loc[0, 'y'])
#row.append(pdf.loc[0, 'b_h2o_loss'])
#row.append(pdf.loc[0, 'b_nh3_loss'])
#row.append(pdf.loc[0, 'y_h2o_loss'])
#row.append(pdf.loc[0, 'y_nh3_loss'])
#row.append(0.5)
#row.append(0.8)




"""
cd /home/jia/Documents/Code/waterlooms/lin_work/Data/pepfind/dia_dp_data_formating

python
"""
import pickle
import numpy as np

samplefile = "cas1_DIAUmpire_timeseries.pkl"
with open(samplefile, "rb") as load_file:
	realX=pickle.load(load_file)
	realY=pickle.load(load_file)
	decoyX=pickle.load(load_file)
	decoyY=pickle.load(load_file)

height = 66  # pep len
width = 42  # pep features not used
time_window = 3  # number of previous or following scans to be considered
depth = time_window*2+1

#realX.shape (15295, 19404)

X = np.vstack((realX, decoyX))
# X.shape (36557, 19404)


Y = np.vstack((realY, decoyY))
# combine X,Y and shuffle and then detach X,Y
XY = np.hstack((X,Y))
np.random.shuffle(XY)
np.random.shuffle(XY)

X_data = XY[:, 0:(XY.shape[1]-1)]
# X_data.shape (36557, 19404)

"""
36557 training samples of this shape: 
Data is just arranged contiguously 

>>> 19404 / 66
294.0
>>> 294/7
42.0
>>> 42/42
1.0
"""


training_num = round(0.9*X_data.shape[0]) #scalar

X_test = X_data[training_num:training_num*2, :]
# X_test.shape (3656, 19404)

X_train = X_data[:training_num, :]
# X_train.shape (32901, 19404)

X_train = X_train.reshape(-1, depth, height, width)
# X_train.shape (32901, 7, 66, 42)

X_test = X_test.reshape(-1, depth, height, width)
# (3656, 7, 66, 42)



X_train = np.swapaxes(X_train, 1, 2)
X_train = X_train.reshape(-1, height, depth, width, 1)





#####=====#####=====#####=====#####=====#####=====#####=====#####=====#####=====
# Tensorflow reshaping test
import tensorflow as tf


encoding = [0.1] * 20
encoding.extend([0.2] * 10)
# single peptide representation
encoding = encoding * 40

# peptide w/ charges 1,2 representation
encoding.extend([item * 2 for item in encoding])

# peptide matched to 100 scans representation
encoding = encoding * 100

# shaped final tensor product
encoding_tensor = tf.reshape(encoding, shape=(100, 2, 40,30))





encoding.extend([item * 3 for item in encoding])
encoding_tensor = tf.reshape(encoding, shape=(2, 100, 2, 40, 30))
# num_samples, num_scans, num_charges, peptide_len, feature_len




tmp = tf.concat([encoding_tensor, encoding_tensor], axis=0)






"""
>>> len(encoding_tensor)
100
>>> len(encoding_tensor[0])
40
>>> len(encoding_tensor[0][0])
30
"""









row.append(0.1) #b_mono			intensity
row.append(0.2) #y_mono			intensity
row.append(0.3) #b+1			intensity
row.append(0.4) #y+1			intensity
row.append(0.5) #b-H2O			intensity
row.append(0.6) #b-NH3			intensity
row.append(0.7) #y-H2O			intensity
row.append(0.8) #y-NH3 			intensity
row.append(0.9) #precursor 		intensity
row.append(0.95)#precursor+1 	intensity
# Done row operations
# intensity is normalized to the maximum intensity within the scan



tmp_tensor = tf.concat([row, row], 0)
tf.reshape(tmp_tensor, shape=(2, 30))
tf.reshape(tmp_tensor, shape=(2, 2, 15))
tf.reshape(tmp_tensor, shape=(1, 1, 2, 2, 15))





import argparse
from bisect import bisect_left
from enum import Enum
import numpy as np
import pandas as pd
import pickle
from pyteomics import mass, mgf, mzxml
import re
import sys

mgf_input = "/data/PXD005573/Fig1_MP-DIA-120min120kMS1-22W30k-8dppp_MHRM_R01_Q1.mgf"
mgf_in = mgf.IndexedMGF(mgf_input)

spec_id = 88505

pseudospectra = mgf_in.get_by_index(spec_id)






python /waterlooms/src/main/python/scoring/ExtractTrainingFeatures.py -start_index 0 -end_index 1000

python /waterlooms/src/main/python/scoring/ExtractTrainingFeatures.py -start_index 1001 -end_index 2000

python /waterlooms/src/main/python/scoring/ExtractTrainingFeatures.py -start_index 2001 -end_index 3000

python /waterlooms/src/main/python/scoring/ExtractTrainingFeatures.py -start_index 3001 -end_index 4000
python /waterlooms/src/main/python/scoring/ExtractTrainingFeatures.py -start_index 4001 -end_index 5000
python /waterlooms/src/main/python/scoring/ExtractTrainingFeatures.py -start_index 5001 -end_index 6000
python /waterlooms/src/main/python/scoring/ExtractTrainingFeatures.py -start_index 6001 -end_index 7000
python /waterlooms/src/main/python/scoring/ExtractTrainingFeatures.py -start_index 7001 -end_index 8000
python /waterlooms/src/main/python/scoring/ExtractTrainingFeatures.py -start_index 8001 -end_index 9000
python /waterlooms/src/main/python/scoring/ExtractTrainingFeatures.py -start_index 9001 -end_index 10000
python /waterlooms/src/main/python/scoring/ExtractTrainingFeatures.py -start_index 10001 -end_index 11000
python /waterlooms/src/main/python/scoring/ExtractTrainingFeatures.py -start_index 11001 -end_index 12000

python /waterlooms/src/main/python/scoring/ExtractTrainingFeatures.py -start_index 12001 -end_index 13000
python /waterlooms/src/main/python/scoring/ExtractTrainingFeatures.py -start_index 13001 -end_index 14000

python /waterlooms/src/main/python/scoring/ExtractTrainingFeatures.py -start_index 14001 -end_index 15000

python /waterlooms/src/main/python/scoring/ExtractTrainingFeatures.py -start_index 15001 -end_index 16000

python /waterlooms/src/main/python/scoring/ExtractTrainingFeatures.py -start_index 16001 -end_index 16655






import numpy as np
import pickle
training_pickle = "/data/PXD005573/training_features_16001_16655.pkl"
with open(training_pickle, "rb") as load_file:
	realX = pickle.load(load_file)
	realY = pickle.load(load_file)
	decoyX = pickle.load(load_file)
	decoyY = pickle.load(load_file)



realX_arr = np.array(realX)


flat_array = np.concatenate(realX_arr).ravel()


import tensorflow as tf

encoding_tensor = tf.reshape(flat_array, shape=(len(realX_arr), 100, 2, 40,30))


# For a given tensor, slice the first 20 elements
encoding_tensor[0][0][:][0][20:]

encoding_tensor[0][0][:][0][20:][20:]


encoding_tensor[0][0:10][:][:20][20:]



import sys
np.set_printoptions(threshold=sys.maxsize)

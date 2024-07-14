import sys
import os
from matplotlib import pyplot as plt
from math import ceil
"""
Caroline's Evaluation
# this script evaluates the quality of my match by plotting 
# 1) all peptide score distribution
# 2) high scoring peptide score distribution
# 3) up to 15% fdr distribution

"""



# post: postfix added after plot names
post = "sample-2"
# folder: output folder
# folder = "../data/Logs/windowsLogs/"
folder = "../../data/Logs/"

# if os.path.exists("../data/Logs/target.png"):
#   os.remove("../data/Logs/target.png")

# if os.path.exists("../data/Logs/decoy.png"):
#   os.remove("../data/Logs/decoy.png")


logFile = "../../data/Logs/log.txt"

# logFile = "results/log-mid-2.txt"
fh = open(logFile, "r")
source = fh.read()
fh.close()

lst = source.split('>')

# print # of unique peps and # of matched peps
pepsCount = lst.pop(0)
print(pepsCount)

# fdr
stop = False 
targetCount = 0
decoyCount = 0 
target = []
decoy = []
# together = []
def get_score_from_pep(pep):
    return float(pep.split('\n')[2])

def plot_data(target, decoy, start, file):
    plt.clf()
    score_bins = range(start, ceil(max_score) + 5, 1)
    
    (n_d, bins, patches) =plt.hist([target, decoy],
             bins = score_bins, label=['Target', 'Decoy'])

    plt.legend();
    plt.savefig(file)

max_score = get_score_from_pep(lst[0])

# separate into target vs decoy
for pep in lst:
    # together.append(get_score_from_pep(pep))
    if 'DeBruijn' in pep:
        decoy.append(get_score_from_pep(pep))
        # if not stop: 
        decoyCount+=1
        
    else:
        target.append(get_score_from_pep(pep))
        # if not stop: 
        targetCount+=1

# target vs decoy
print("all scores:")
print("target: " + str(len(target)))
print("decoy: " + str(len(decoy)))
# full score distribution plot
plot_data(target, decoy, 0, folder+"AllScores"+ post+".png")

# FDR
fdr_decoy = []
fdr_target = []
res = []
fdrs = []
fdr_threshold = .15 # plot up to 15% FDR
lowest = 0

for pep in lst:
    if 'DeBruijn' in pep:
        fdr_decoy.append(get_score_from_pep(pep))
        
    else:
        fdr_target.append(get_score_from_pep(pep))

    fdr = len(fdr_decoy)/(len(fdr_target))
    fdrs.append(fdr)
    res.append(len(fdr_target))

    if fdr > fdr_threshold:
        lowest = ceil(get_score_from_pep(pep) - 5) 
        break    

plt.clf()
plt.plot(res, fdrs)

plt.savefig(folder+"Fdr"+ post+".png")

# FDR 15% score distribution plot
plot_data(fdr_target, fdr_decoy, lowest, folder+"HighScores"+ post+".png")

import pandas as pd
import os

# this script shows high-scoring umpire outputs that my program fail to identify.
# the script filters umpire output to find 1% fdr results of purey tryptic peptides
# within selected mass window

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', None)
pd.options.mode.chained_assignment = None

# fdr cutoff of umpire peptides to consider
umpire_fdr = 0.01
# fdr cutoff of my peptides
my_fdr = .95
print("threshold: " + str(my_fdr))

# fn: given df, return df with rolling fdr column
def add_fdr_df(df, mode, fdr, trunk=True):
    colname = 'name' if mode=='mine' else 'Protein' 
    df['false_match'] = df[colname].str.startswith('DeBruijn').astype(int).cumsum()
    # print(false_match)
    df.reset_index(inplace=True,drop=True)
    # print(df.head(25))
    df[mode+'_fdr'] = df['false_match']/(df.index -df['false_match'])
    if trunk:
        last = df[df[mode+'_fdr']<=fdr].index[-1]+1
        return df.truncate(0,last)

    return df


# read Umpire results
charge = 2


# read my results, can use window high/low to filter umpire results
# to get comparable set of peptides

# windowLow = 477.0
# windowHigh = 491.0
# logFile = "results/log480-3.txt"

# windowLow = 733.0
# windowHigh = 755.0
# logFile = "results/log750-2.txt"

# windowLow = 1117.0 
# windowHigh = 1650.0
# logFile = "../../data/Logs/log.txt"
# logFile = "results/log1300-2.txt"

# windowLow = 634.0 
# windowHigh = 825.0
# logFile = "../../data/Logs/log.txt"
logFile = "results/log-sample-2.txt"


df = pd.read_csv("../../data/umpire_output/HELA-DIA-45C-2.tsv", sep="\t")
df = df[['ScanTime(Min)', 'Precursor','Charge', 'Peptide',
       'Protein', 'EValue']]
# print(df.dtypes)


# filter dataframe to get right mass window
# df = df[(df['Charge']==charge) & (df['Precursor']>=windowLow) & (df['Precursor']<=windowHigh)]
df = df[(df['Charge']==charge) & \
((df['Precursor']>=477.0) & (df['Precursor']<=491.0)) | \
((df['Precursor']>=553.0) & (df['Precursor']<=567.0)) |
((df['Precursor']>=733.0) & (df['Precursor']<=755.0)) |
((df['Precursor']>=923.0) & (df['Precursor']<=971.0)) |
((df['Precursor']>=1117.0) & (df['Precursor']<=1650.0))]

# filter for pure tryptic peptides
df_trp = df[((df['Peptide'].str.startswith('K')) | (df['Peptide'].str.startswith('R'))) & \
    ((df['Peptide'].str.endswith('K')) | (df['Peptide'].str.endswith('R')))]

# remove peptides with M adjustment since my code does not adjust for M
df_trp = df_trp[~df_trp['Peptide'].str.contains('M\+', regex=True)]

# remove head and tail and the last character (since I want no R/K in middle)
df_trp['PeptideProcessed'] = df_trp['Peptide'].str[2:-2]
# remove numbers and special characters
df_trp['PeptideProcessed'] = df_trp['PeptideProcessed'].replace('\+[0-9]*\.[0-9]*','',regex=True)

# make sure peptide does not start with P (not tryptic)
df_trp = df_trp[~df_trp['PeptideProcessed'].str.startswith('P')]
# make sure peptides has no R/K in middle (not tryptic)
df_trp = df_trp[~df_trp['PeptideProcessed'].str.contains('(K|R).*(K|R)', regex=True, case=False)]

# dedup
df_trp.drop_duplicates(subset=['PeptideProcessed'],inplace=True)

df_trp = df_trp.sort_values(['EValue'])

df_trp.reset_index(inplace=True,drop=True)
df_top = df_trp.truncate(0,5000)
# df_top.to_csv('umpire_top.csv')

# now read my results into df
fh = open(logFile, "r")
source = fh.read()
fh.close()

lst = source.split('>')

pepsCount = lst.pop(0)

# where things are stored
dicts = list()

for pep in lst:
    pep_lst = pep.split('\n')
    cur_dict = {
        'name' : pep_lst[0],
        'pep' : pep_lst[1],
        'score' : float(pep_lst[2]),
        'times' : pep_lst[3]
    }
    dicts.append(cur_dict)

my_match = pd.DataFrame(dicts)

print(my_match.shape)
# print(my_match.head(40))
my_match.reset_index(inplace=True,drop=True)
my_match_top = my_match.truncate(0,5000)
my_match_top.to_csv('my_match.csv')
my_match = add_fdr_df(my_match, 'mine', my_fdr)
# print(my_match.head(40))
print(my_match.shape)
print(my_match['score'].iloc[-1])

# join by peptide
umpire_res = pd.merge(my_match, df_trp, how='right',left_on='pep', right_on='PeptideProcessed')
umpire_res = umpire_res.sort_values(['EValue'])

umpire_res = add_fdr_df(umpire_res, 'umpire', umpire_fdr)

umpire_res = umpire_res[['name','pep','score', 'times', 'EValue', 'Peptide', 'Protein', 'ScanTime(Min)', 'umpire_fdr']]

print(umpire_res)

print('scored by umpire:')
# check ones that are in umpire
in_umpire = umpire_res["pep"].notna().sum()
print(in_umpire)
# check ones that are not in umpire
not_umpire = umpire_res["pep"].isna().sum()
print(not_umpire)

umpire_decoy = umpire_res[umpire_res['Protein'].str.startswith('DeBruijn')]['pep'].isna().sum()
print(umpire_decoy)

# check ones that are in my results
my_res = pd.merge(df_trp, my_match, how='right',right_on='pep', left_on='PeptideProcessed')
my_res = my_res.sort_values(['score'], ascending=False)

my_res = add_fdr_df(my_res, 'mine', umpire_fdr)

my_res = my_res[['name','pep','score', 'EValue', 'Peptide', 'Protein','mine_fdr']]

# print(my_res)

print('scored by me:')
in_mine = my_res["Peptide"].notna().sum()
print(in_mine)
# check ones that are not in umpire
not_mine = my_res["Peptide"].isna().sum()
print(not_mine)


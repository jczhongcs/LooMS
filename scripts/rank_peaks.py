import pandas as pd
import numpy as np

# this script calculates the half rank of a csv file of peaks,
# the order of peaks in given csv is preserved
# this code is very slow but I haven't found why

file = "../../data/spectrums/sample/spectrum.csv"
out = "../../data/spectrums/sample/spectrum_half_rank.csv"

df = pd.read_csv(file, header=None)
print(df)
df.columns = ['window','mz','rt','int']
df_sorted = df.sort_values(by=['window', 'rt', 'int'], ascending=False)

# first, calculate rank of each peak
df_sorted['rank'] = df_sorted.groupby(['window', 'rt']).cumcount()
print(df_sorted)

# this merge is used to get back original ordering
res = df.merge(df_sorted, how='left', on=['window','mz','rt', 'int'])
print(res)

res.to_csv(out, header=None)

# sometimes I stop here and run the second part later
# by reading in the csv generated above

# after sorting, get half rank by getting rank of mz/2
# this part can be optimized as we cna simply find max; no need to sort

# df = pd.read_csv(file, header=None, index_col=0)
# df.columns = ['window','mz','rt','int', 'rank']
df = res
df['half_rank'] = np.nan
groups = df.groupby(['window', 'rt'])

for index,row in df.iterrows():
    window = row['window']
    rt = row['rt']
    half_int = row['int']/2
    print(index)
    cur_group = groups.get_group((window, rt))
    ranked = cur_group.iloc[(cur_group['int'] - half_int).abs().argsort()].reset_index(inplace=False, drop=True)
    df.loc[index, 'half_rank'] = ranked['rank'][0]
    # print(ranked['rank'][0])
    # break

print(df)
df.to_csv(out)
# df['half_rank'] = 



# print(df)
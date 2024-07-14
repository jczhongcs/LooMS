import numpy as np
import pandas as pd
import math

import os
class RT:
    def __init__(self):

        self.iRT_max = 0
        self.RT_max = 0
        self.iRT_min = 0
        self.RT_min = 0

    def read_tsv(self,input_file_path, cols):
        results=pd.read_csv(input_file_path, sep='\t',encoding='ISO-8859-1', usecols=cols)
        return results

    def write_tsv(self,results, output_file_path):
        results.to_csv(output_file_path, sep='\t', index=False)

    def read_csv(self,input_file_path, cols):
        results=pd.read_csv(input_file_path, sep=',',encoding='ISO-8859-1', usecols=cols)
        return results



    def get_RTiRT_csv(self,peptide_input_file_path):
        #peptide_cols = [ 'RT', 'iRT' , 'decoyTarget']
        peptide_cols = [ 'RT', 'iRT' ]
        # peptide_cols = ['peptide', 'RT', 'iRT']
        peptide_input = self.read_tsv(peptide_input_file_path, peptide_cols)
        peptide_input = peptide_input.dropna()#remove nan rows

        # peptide_input = peptide_input.sort_values(by=['msoneQualityScore'])
        # peptide_input = peptide_input.tail(30000)
        global iRT_max, RT_max, iRT_min, RT_min
        peptide_input = peptide_input.sort_values(by=['iRT'])
        iRT_max = peptide_input.iloc[[len(peptide_input)-1]]['iRT'].values[0]
        iRT_min = peptide_input.iloc[[0]]['iRT'].values[0]
        peptide_input = peptide_input.sort_values(by=['RT'])
        RT_max = peptide_input.iloc[[len(peptide_input) - 1]]['RT'].values[0]
        RT_min = peptide_input.iloc[[0]]['RT'].values[0]
        return peptide_input


    def MSE_DP(self,X, Y, N):
        dx = (RT_max - RT_min) / N
        dy = (iRT_max - iRT_min) / N
        step = 4
        s_ij = [[0 for x in range(N)] for y in range(N)]
        a = [0 for x in range(N)]
        b = [0 for y in range(N)]

        # initialization
        for i in range(N):
            a[i] = dx * i + RT_min
        for j in range(N):
            b[j] = dy * j + iRT_min

        # divide into grid
        grid = [[[] for x in range(N)] for y in range(N)]
        for w in range(len(X)):
            # check grid
            x = X[w]
            y = Y[w]
            i = int((x - RT_min) / dx)
            j = int((y - iRT_min) / dy)
            i_step = step
            j_step = step
            if i + i_step - 1 > N:
                i_step = N - i + 1
            if j + j_step - 1 > N:
                j_step = N - j + 1
            for p in range(i_step):
                if i+p-1 < 0:
                    continue
                for q in range(j_step):
                    if j+q-1 < 0:
                        continue
                    grid[i+p-1][j+q-1].append(w)

        # DP
        for i in range(N):
            grid_i = grid[i][0]
            s = 0
            for p in range(len(grid_i)):
                x = X[grid[i][0][p]]
                y = Y[grid[i][0][p]]
                s += math.exp(-(((x - a[i]) / dx) ** 2 + ((y - b[0]) / dy) ** 2))
            s_ij[i][0] = s
        for j in range(N):
            grid_j = grid[0][j]
            s = 0
            for p in range(len(grid_j)):
                x = X[grid[0][j][p]]
                y = Y[grid[0][j][p]]
                s += math.exp(-(((x - a[0]) / dx) ** 2 + ((y - b[j]) / dy) ** 2))
            s_ij[0][j] = s
        for i in range(1, N):
            for j in range(1, N):
                grid_ij = grid[i][j]
                s0 = s_ij[i-1][j] # left
                s1 = s_ij[i-1][j-1] # left-bottom
                s2 = s_ij[i][j-1] # bottom
                if s0 > max(s1, s2):
                    s = s0
                elif s1 > max(s0, s2):
                    s = s1
                else:
                    s = s2
                for p in range(len(grid_ij)):
                    x = X[grid[i][j][p]]
                    y = Y[grid[i][j][p]]
                    s += math.exp(-(((x - a[i]) / dx)**2 + ((y - b[j]) / dy)**2))
                s_ij[i][j] = s
        # backtracking
        i = N - 1
        j = N - 1
        k = [[i, j]]
        while i >= 1 and j >= 1:
            s0 = s_ij[i - 1][j]  # left
            s1 = s_ij[i - 1][j - 1]  # left-bottom
            s2 = s_ij[i][j - 1]  # bottom
            if s0 > max(s1, s2):
                k.append([i - 1, j])
                i = i - 1
            elif s1 > max(s0, s2):
                k.append([i - 1, j - 1])
                i = i - 1
                j = j - 1
            else:
                k.append([i, j - 1])
                j = j - 1
        # k.pop()
        pos_x = []
        pos_y = []
        for e in k:
            # pos_x.append(a[e[0]])
            # pos_y.append(b[e[1]])
            pos_x.append(round(a[e[0]],1))
            pos_y.append(round(b[e[1]],1))
        return  pos_x, pos_y
    def run(self, feature_file):
        filepath = os.path.splitext(feature_file)[0]

        #R1_input_file_path = "/Users/jianzhong/Documents/uwaterloo_dia/DIA Data/humanR01.csv"
        # R1_input_file_path = filepath+"forRTminfit_maxfit.csv"
        R1_input_file_path = filepath+".csv"
        print(R1_input_file_path)
        result = self.get_RTiRT_csv(R1_input_file_path)
        # result = get_iRT_csv(R1_input_file_path)
        RT = result['RT'].values
        iRT = result['iRT'].values
        print("# of peptides: ", len(result))
        print("RT MAX:", RT_max)
        print("iRT MAX:", iRT_max)
        print("RT MIN:", RT_min)
        print("iRT MIN:", iRT_min)

        N = 3000
        predctedRT, predctediRT = self.MSE_DP(RT, iRT, N)

        percentile_list = pd.DataFrame(
            {'predctedRT': predctedRT,
             'predctediRT': predctediRT
             }).drop_duplicates()


        df = percentile_list.groupby(by=['predctediRT']).agg(autoRt=('predctediRT','min'),
                                                             fit_min=('predctedRT','min'),
                                                             fit_max=('predctedRT','max'))

        df.to_csv(filepath+"_percentile_list.csv", sep='\t', index=False)

import argparse # add to beginning of file
parser = argparse.ArgumentParser(description='Train on features extracted from mzXML.')
parser.add_argument("-features",
                    help="Features computed from FeatureDetect."
                    )
parser.add_argument("-parameters",
                    default="/data/featuredetect.params",
                    help="FeatureDetect parameter file."
                    )
args = parser.parse_args()
# Invoke the training
rtC = RT()
# TO-DO: change to training.run(args.features, parsed)
rtC.run(args.features)
#fig.savefig('/Users/mac/Desktop/dia_webapp/resources/output/RT_alignment_dp_n' + str(N))
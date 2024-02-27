'''
Purpose: Produce scatter plot of hap1 vs hap2 QV's, with line connecting "before" and "after" polishing
Author: Mira Mastoras, mmastora@ucsc.edu
Usage: python3 QV_scatter_plot.py -i input_data.csv -o outfile

python3 QV_scatter_plot.py -i /Users/miramastoras/Desktop/HPRC_polishing_QC.k21.csv -o /Users/miramastoras/Desktop/HPRC_polishing_QC.k21.QV_scatter.png
'''
import matplotlib.pyplot as plt
import numpy as np
import argparse
import pandas as pd
import matplotlib.patches as mplpatches

# parse inputs
parser = argparse.ArgumentParser()
parser.add_argument('--input_csv', '-i',
                    type=str, action='store', help='input csv with one row per sample and columns with data inputs.')
parser.add_argument('--outfile', '-o',
                    type=str, action='store', help='output filename')
args = parser.parse_args()

# read in data
data=pd.read_csv(args.input_csv, sep=',', header=0)

# subset by
# gather QVs for each sample before and after polishing
# QV_dict_polished={}
# QV_dict_raw={}
#
# for data_index, data_row in data.iterrows():
#     sample_id = data_row[0]
#
#     if

# Plot the figure
figureHeight=3
figureWidth=8
plt.figure(figsize=(figureWidth,figureHeight))
plt.rcParams.update({'font.size': 7.5})
#
x = np.array([5,7,8,7,2,17,2,9,4,11,12,9,6])
y = np.array([99,86,87,88,111,86,103,87,94,78,77,85,86])

plt.scatter(x, y)

# function to add arrow on a graph
plt.arrow(5,99,17,86,width=0.2)

plt.show()
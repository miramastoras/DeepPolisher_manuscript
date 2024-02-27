'''
Purpose: Produce scatter plot of hap1 vs hap2 QV's, with line connecting "before" and "after" polishing
Author: Mira Mastoras, mmastora@ucsc.edu
Usage: python3 QV_scatter_plot.py -i input_data.csv -o outfile
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

# gather QVs for each sample before and after polishing 
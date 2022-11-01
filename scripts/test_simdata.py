#! /usr/bin/env python

from treeswift import *
from math import *
from random import random
from scipy.optimize import fsolve
import sys
sys.path.append("/Users/gillianchu/raphael/repos/problin")
from problin_libs.sequence_lib import read_sequences
import numpy as np

outdir = "/Users/gillianchu/raphael/repos/problin/MP_inconsistent/"
m = 10


for k in [20,30,40,50,100,200,300,400,500,1000,5000]:
    count = 0
    outfile = outdir + "seqs_m" + str(m) + "_k" + str(k) + ".txt"

    count_db = 0
    count_bd = 0
    count_ac = 0
    count_ca = 0
    # calculate average number of changes for each sequence
    S = read_sequences(outdir + "seqs_m10_k" + str(k) + ".txt")
    d = dict()
    # D = S[0]
    for D in S:
    # print(D) 
        for x in D:
            seq = D[x]
            # if x not in d.keys():
            #     d[x] = []
            d[x] = sum([int(y != 0) for y in seq]) / len(seq)
        if d['b'] < d['d']:
            count_db += 1
        elif d['b'] > d['d']:
            count_bd += 1
        if d['c'] < d['a']:
            count_ac += 1
        elif d['c'] > d['a']:
            count_ca += 1
    
    #for x in D:
    #    print(k, "characters; sequence", x, np.mean(np.array(d[x])), "changes")

    print(k, "characters", "diff bd", (count_db - count_bd)/1000, "diff ca", (count_ac - count_ca)/1000 )

import numpy as np
from scipy import optimize
import math
from math import log, exp

import sys
sys.path.append("/Users/gillianchu/raphael/repos/problin/problin_libs")
# from distance_based_lib import ML_pairwise_estimate
from ml import wrapper_felsenstein
from sequence_lib import read_sequences

true_tree = ''
T = '[&R] ((0,1),(2,3));'

msa = np.array([[1,0,2], # species 1
                [1,0,0], # species 2
                [0,1,2], # species 3
                [0,1,0]] # species 4
            )

# prob of transitioning 0 -> 0, 0 -> 1, 0 -> 2
# Q = [{0: 0.2, 1:0.3, 2:0.5}, # site 1
#      {0: 0.2, 1:0.3, 2:0.5}, # site 2
#      {0: 0.2, 1:0.3, 2:0.5}, # site 3
#     ]

# for each k
#   for each of the 1k replicates:
#       generate all fifteen topologies
#       optimize the branch lengths for each of those topologies
#       pick the best one
#   compute the percent correct for each k
# q = 0.1

# x_star, f_star = wrapper_felsenstein(T, Q, msa, 3, root_edge_len=0.2)
# print(x_star, f_star)

# enumerate all fifteen topologies for 4 leaves
topologies = ["((a,b),(c,d));","((a,c),(b,d));","((a,d),(b,c));",
            "(a,(b,(c,d)));","(a,(c,(b,d)));","(a,(d,(b,c)));",
            "(b,(a,(c,d)));","(b,(c,(a,d)));","(b,(d,(a,c)));",
            "(c,(a,(b,d)));","(c,(b,(a,d)));","(c,(d,(a,b)));",
            "(d,(a,(b,c)));","(d,(b,(a,c)));","(d,(c,(a,b)));"]

true_topology = '((a,b),(c,d));',

m = 10
with open("ML_felsenstein_results.txt",'w') as fout:
    for k in [5000]: #,30,40,50,100,200,300,400,500,1000,5000]:
        correct = 0
        total = 0

        S = read_sequences("/Users/gillianchu/raphael/repos/problin/MP_inconsistent/seqs_m10_k" + str(k) + ".txt")
        n_total = 0
        n_correct = 0
        for D in S[:20]:

            Q = []
            for i in range(k):
                q = {j+1:1/m for j in range(m)}
                q[0] = 0
                Q.append(q)

            top_x_star =  -float("inf")
            top_f_star =  -float("inf")
            top_topology = None
            for topology in topologies: 
                x_star, f_star = wrapper_felsenstein(topology, Q, D, use_log=True)
                # print(topology, f_star)
                if f_star > top_f_star:
                    top_f_star = f_star
                    top_x_star = x_star
                    top_topology = topology
            if top_topology == true_topology:
                correct += 1
            total += 1
            print(k, "characters", correct, total, top_topology)
        print(k, "characters", correct/total, "true topologies selected.")
        # fout.write(str(k) + " " + str(n_correct/n_total) + "\n")

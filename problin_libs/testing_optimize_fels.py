
### Testing
import numpy as np
import sys
sys.path.append("/Users/gillianchu/raphael/repos/problin")
from problin_libs.sequence_lib import read_sequences
from problin_libs.ml import wrapper_felsenstein, optimize_len
from treeswift import *
# k=20
m=10 # 10
Q = []
datadir = "/Users/gillianchu/raphael/repos/problin/MP_inconsistent/"

# k = 3

# for k in [3,10,20,30,40,50,100,200,300,400,500,1000,5000]:
for k in [5000]:

    S = read_sequences(datadir + "seqs_m10_k" + str(k) + ".txt")
    
    corr_d = dict()
    wrong_d = dict()

    num_true, total, num_incorrect, num_same = 0, 0, 0, 0
    for i in range(k):
        q = {j+1:1/m for j in range(m)}
        # q = {j+1:1 for j in range(m)}
        q[0] = 0
        Q.append(q)

    for idx, D in enumerate(S):
        print("k", k, "rep", idx)
        for x in D:
            print(x, len([y for y in D[x] if D[x][y] == 0] ) )
        # # true_tree = "((a:0.0360971597765934,b:3.339535381892265)e:0.0360971597765934,(c:0.0360971597765934,d:3.339535381892265)f:0.0360971597765934)r;"
        # true_tree = "(a:0.0360971597765934,b:3.339535381892265)e:0.0360971597765934;"
        # opt_likelihood, opt_tree, opt_branches = wrapper_felsenstein(true_tree, Q, D, use_log=True, optimize_branchlengths=True, initials=5)
        # print(opt_likelihood, opt_tree, opt_branches)
        # est_tree = "((a:0.00045143472400357813,b:7.683405830494674)e:0.00010001000133352235,(c:0.00010001000133352235,d:7.1648549920570614)f:0.00010001000133352235)r:1.132411178316622;"
        
        true_tree = "((a:0.0360971597765934,b:3.339535381892265)e:0.0360971597765934,(c:0.0360971597765934,d:3.339535381892265)f:0.0360971597765934)r:0.0;"
        # noopt_likelihood, noopt_tree, noopt_branches = wrapper_felsenstein(true_tree, Q, D, use_log=True, optimize_branchlengths=True, init_tree=true_tree)
        # print(noopt_likelihood, noopt_tree, noopt_branches)
        
        # noopt_likelihood, noopt_tree, noopt_branches = wrapper_felsenstein(true_tree, Q, D, use_log=True, optimize_branchlengths=False)
        # print(noopt_likelihood, noopt_tree, noopt_branches)

        # x0 = (0.0360971597765934, 3.339535381892265, 0.0360971597765934)
        # fout = optimize_len(k, D['a'], D['b'], x0)
        # print(sorted(fout, key=lambda x: x[0])[-1])
        # break
        # noopt_likelihood_1, noopt_tree_1, noopt_branches_1 = wrapper_felsenstein(true_tree, Q, D, use_log=True, optimize_branchlengths=False)
        # true_tree = "(a:0.0360971597765934,b:3.339535381892265)e:1.0;"
        # noopt_likelihood_2, noopt_tree_2, noopt_branches_2 = wrapper_felsenstein(true_tree, Q, D, use_log=True, optimize_branchlengths=False)
        # true_tree = "(a:0.0360971597765934,b:3.339535381892265)e:10.0;"
        # noopt_likelihood_3, noopt_tree_3, noopt_branches_3 = wrapper_felsenstein(true_tree, Q, D, use_log=True, optimize_branchlengths=False)
        # print(noopt_likelihood_1, noopt_likelihood_2, noopt_likelihood_3)

        # print(noopt_branches_1, noopt_branches_2)
        # print(fels_true_n)
        # false_tree = "((a:0.0360971597765934,d:3.339535381892265)e:0.0360971597765934,(c:0.0360971597765934,b:3.339535381892265)f:0.0360971597765934)r;"
        # fels_false_likelihood = wrapper_felsenstein(false_tree, Q, D, use_log=True, optimize_branchlengths=False)
        
        # print(opt_likelihood > noopt_likelihood)
        # for x, y in zip(opt_branches, noopt_branches):
        #     print(x, y)

        # if idx > 2:
        #     break
        # if opt_likelihood > noopt_likelihood:
        #     num_true += 1
        # elif opt_likelihood < noopt_likelihood: 
        #     num_incorrect += 1
        # else:
        #     num_same +=1

        # total += 1
        # break  
    # print(k, "characters", num_true/total, "preferred true topology.", num_same/total, "num same")
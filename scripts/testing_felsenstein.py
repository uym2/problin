import numpy as np
import dendropy
import math
from random import random,seed
from math import log,exp
from scipy import optimize


def prob_same(node_likelihood, char, site, curr_node, use_log):
    # print("prob same")
    c = curr_node.child_nodes()[0]
    tp = math.exp(-get_branchlen(c))
    return tp * node_likelihood[c][site][char]

def prob_change(msa, q_dict, node_likelihood, site, curr_node, child_states, use_log):
    all_child_prob = 1.0
    for c in curr_node.child_nodes():
        # print("[prob_change]: call to get_char()", c)
        if c.is_leaf():
            char = get_char(msa, c, site)
            if char != 0: 
                # print(q_dict[site], site, char)
                q_ialpha = q_dict[site][char] 
                tp = q_ialpha * (1 - math.exp(-get_branchlen(c)))
                all_child_prob *= tp * node_likelihood[c][site][char]
            else:
                q_ialpha = q_dict[site][0]
                tp = math.exp(-get_branchlen(c))
                all_child_prob *= tp * node_likelihood[c][site][0]
        else:
            for char in node_likelihood[c][site].keys():
                # print("[prob_change]: internal child char", char)
                if char != 0: 
                    # print(q_dict[site], site, char)
                    q_ialpha = q_dict[site][char] 
                    if node_likelihood[c][site][char] > 0:
                        tp = q_ialpha * (1 - math.exp(-get_branchlen(c)))
                        # print("[prob_change]: tp", tp)
                        all_child_prob += tp * node_likelihood[c][site][char] # GC: Edited this
                else:
                    q_ialpha = q_dict[site][0]
                    if node_likelihood[c][site][0] > 0:
                        tp = math.exp(-get_branchlen(c))
                        all_child_prob += tp * node_likelihood[c][site][0] # GC: Edited this
                # print("[prob_change]:", all_child_prob)

    return all_child_prob



def likelihood_under_n(node_likelihood, n, site, msa, q_dict, is_root, use_log):
    # print("[likelihood_under_n]", n, "site", site)
    child_states = set()
    # print(n)
    if n not in node_likelihood:
        node_likelihood[n] = dict()
        node_likelihood[n][site] = dict()
        
    child_states = []
    for child in n.child_nodes():
        if child.is_leaf():
            child_states.append(get_char(msa, child, site))
        else:
            for x in node_likelihood[child][site]:
                state_prob = node_likelihood[child][site][x]
                if state_prob > 0.0:
                    child_states.append(x)

    # print("[likelihood_under_n]: child_states", child_states)
    parent_poss_states = dict()
    if 0 in set(child_states): # probability 0 -> 0
        if len(set(child_states)) == 1: # both children are state 0 
            tmp = prob_same(node_likelihood, 0, site, n, use_log)
            if is_root:
                parent_poss_states[0] = tmp
            else:
                parent_poss_states[0] = tmp **2
        else: 
            for c in child_states: # probability c -> c != 0
                parent_poss_states[c] = 0.0
            # probability 0 -> c (alpha)
            parent_poss_states[0] = prob_change(msa, q_dict, node_likelihood, site, n, child_states, use_log)  
    else:
        if len(set(child_states)) == 1: # both children are same nonzero state
            c = child_states[0]
            parent_poss_states[c] = 1.0 
            parent_poss_states[0] = prob_change(msa, q_dict, node_likelihood, site, n, child_states, use_log)
        else:
            parent_poss_states[0] = 1.0

    for x in parent_poss_states.keys():
        node_likelihood[n][site][x] = parent_poss_states[x]

    return node_likelihood

def get_branchlen(child_node):
    if child_node.edge_length is None:
        print(child_node.child_nodes())
    return child_node.edge_length

def get_char(msa, leaf_node, site):
    # print(msa)
    return msa[leaf_node.taxon.label][site]

def felsenstein(T, Q, msa, use_log=False): #, root_edge_len=0.0):
    # print("MSA", msa)
    numsites = len(msa[next(iter(msa.keys()))])
    # numsites = len(msa.key[0])

    alphabet = dict()
    for site in range(numsites):
        alphabet[site] = Q[site].keys()

    # print("q_dict", Q)
    nwkt = dendropy.Tree.get(data=T, schema="newick", rooting="force-rooted")
    # print(nwkt)

    # for n in nwkt.leaf_node_iter():
    #     print(n.taxon, ''.join([str(get_char(msa, n, s)) for s in range(numsites)]))

    node_likelihood = dict()

    ## CALCULATE THE LIKELIHOOD
    for n in nwkt.postorder_node_iter():
        # print("node:", n)
        if n.is_leaf(): # n.taxon is not None: # must be a leaf node, set up 
            node_likelihood[n] = dict()
            for site in range(numsites):
                node_likelihood[n][site] = dict()
                for char in alphabet[site]:
                    node_likelihood[n][site][char] = 0.0
                # print("[felsenstein]: in site", site)
                char_state = get_char(msa, n, site)
                node_likelihood[n][site][char_state] = 1.0
            
        elif n.is_internal(): # n.taxon is None: # must be an internal node
            for site in range(numsites):
                if n not in node_likelihood.keys():
                    node_likelihood[n] = dict()
                node_likelihood[n][site] = dict()
                for char in alphabet[site]:
                    node_likelihood[n][site][char] = 0.0

                node_likelihood = likelihood_under_n(node_likelihood, n, site, msa, Q, nwkt.seed_node is n, use_log)
    # print(node_likelihood)
    if use_log:
        # print("Using log.")
        tree_likelihood = 0.0
    else:
        # print("NOT using log.")
        tree_likelihood = 1.0

    if nwkt.is_rooted:
        # print("Tree provided was rooted.")
        # print("Tree is rooted, here is likelihood at root.", node_likelihood[n])
        for site in range(numsites):
            if use_log:
                tree_likelihood += np.log(node_likelihood[n][site][0])
            else:
                tree_likelihood *= node_likelihood[n][site][0]
        
    elif not nwkt.is_rooted:
        # print("Tree provided was NOT rooted.")
        for site in range(numsites):
            for rootchar in node_likelihood[n][site].keys():
                prob_rootchar = node_likelihood[n][site][rootchar]
                # print(rootchar, prob_rootchar)
                if prob_rootchar > 0.0: 
                    q_ialpha = Q[site][rootchar]
                    if rootchar == 0:
                        if use_log:
                            tree_likelihood += (-root_edge_len) + np.log(prob_rootchar) # + np.log(q_ialpha) 
                        else:
                            tree_likelihood *= (math.exp(-root_edge_len)) * prob_rootchar # * q_ialpha 
                        
                    else:
                        if use_log:
                            tree_likelihood += np.log((1 - math.exp(-root_edge_len))) + np.log(q_ialpha) + np.log(prob_rootchar)
                        else:
                            tree_likelihood *= ((1 - math.exp(-root_edge_len)) * q_ialpha * prob_rootchar)
    
    return tree_likelihood


def wrapper_felsenstein(T, Q, msa, root_edge_len=0.2, use_log=False, initials=20):
    numsites = len(msa[next(iter(msa.keys()))])
    alphabet = dict()
    for site in range(numsites):
        alphabet[site] = Q[site].keys()

    # print("q_dict", Q)
    nwkt = dendropy.Tree.get(data=T, schema="newick")
    num_edges = len(list(nwkt.postorder_edge_iter()))
    # print(nwkt)

    # for n in nwkt.leaf_node_iter():
    #     print(n.taxon, ''.join([str(get_char(msa, n, s)) for s in range(numsites)]))

    def felsenstein(x): 

        for i, e in enumerate(nwkt.postorder_edge_iter()):
            e.length = x[i]

        alphabet = dict()
        for site in range(numsites):
            alphabet[site] = Q[site].keys()

        node_likelihood = dict()

        for n in nwkt.postorder_node_iter():
            if n.is_leaf(): 
                node_likelihood[n] = dict()
                for site in range(numsites):
                    node_likelihood[n][site] = dict()
                    for char in alphabet[site]:
                        node_likelihood[n][site][char] = 0.0
                    char_state = get_char(msa, n, site)
                    node_likelihood[n][site][char_state] = 1.0
                
            elif n.is_internal(): 
                for site in range(numsites):
                    if n not in node_likelihood.keys():
                        node_likelihood[n] = dict()
                    node_likelihood[n][site] = dict()
                    for char in alphabet[site]:
                        node_likelihood[n][site][char] = 0.0

                    node_likelihood = likelihood_under_n(node_likelihood, n, site, msa, Q, nwkt.seed_node is n, use_log)
        if use_log:
            tree_likelihood = 0.0
        else:
            tree_likelihood = 1.0

        if nwkt.is_rooted:
            for site in range(numsites):
                if use_log:
                    tree_likelihood += np.log(node_likelihood[n][site][0])
                else:
                    tree_likelihood *= node_likelihood[n][site][0]
            
        elif not nwkt.is_rooted:
            for site in range(numsites):
                for rootchar in node_likelihood[n][site].keys():
                    prob_rootchar = node_likelihood[n][site][rootchar]
                    if prob_rootchar > 0.0: 
                        q_ialpha = Q[site][rootchar]
                        if rootchar == 0:
                            if use_log:
                                tree_likelihood += (-root_edge_len) + np.log(prob_rootchar) # + np.log(q_ialpha) 
                            else:
                                tree_likelihood *= (math.exp(-root_edge_len)) * prob_rootchar # * q_ialpha 
                        else:
                            if use_log:
                                tree_likelihood += np.log((1 - math.exp(-root_edge_len))) + np.log(q_ialpha) + np.log(prob_rootchar)
                            else:
                                tree_likelihood *= ((1 - math.exp(-root_edge_len)) * q_ialpha * prob_rootchar)
        
        return tree_likelihood
        
    # def felsenstein(x):
    #     # x is a vector containing all branch lengths
    #     # map branch lengths to tree edges
    #     for i, e in enumerate(nwkt.postorder_edge_iter()):
    #         e.length = x[i]

    #     node_likelihood = dict()

    #     for n in nwkt.postorder_node_iter():
    #         if n.taxon is not None: # must be a leaf node, set up 
    #             node_likelihood[n] = dict()
    #             for site in range(numsites):
    #                 node_likelihood[n][site] = dict()
    #                 for char in alphabet[site]:
    #                     node_likelihood[n][site][char] = 0.0
    #                 char_state = get_char(msa, n, site)
    #                 node_likelihood[n][site][char_state] = 1.0
                
    #         elif n.taxon is None: # must be an internal node
    #             for site in range(numsites):
    #                 if n not in node_likelihood:
    #                     node_likelihood[n] = dict()
    #                 node_likelihood[n][site] = dict()
    #                 for char in alphabet[site]:
    #                     node_likelihood[n][site][char] = 0.0
                    
    #                 node_likelihood = likelihood_under_n(node_likelihood, n, site, msa, Q, use_log)

    #     tree_likelihood = 1.0
    #     for site in range(numsites):
    #         for rootchar in node_likelihood[n][site].keys():
    #             prob_rootchar = node_likelihood[n][site][rootchar]
    #             # print(rootchar, prob_rootchar)
    #             if prob_rootchar > 0.0: 
    #                 q_ialpha = Q[site][rootchar]
    #                 if rootchar == 0:
    #                     if use_log:
    #                         tree_likelihood += (-root_edge_len) + np.log(prob_rootchar) # + np.log(q_ialpha) 
    #                     else:
    #                         tree_likelihood *= (math.exp(-root_edge_len)) * prob_rootchar # * q_ialpha 
                        
    #                 else:
    #                     if use_log:
    #                         tree_likelihood += np.log((1 - math.exp(-root_edge_len))) + np.log(q_ialpha) + np.log(prob_rootchar)
    #                     else:
    #                         tree_likelihood *= ((1 - math.exp(-root_edge_len)) * q_ialpha * prob_rootchar)
    #     return tree_likelihood

#     x_star = []
#     dmax = -log(1/numsites)*2
#     dmin = -log(1-1/numsites)/2
#     bound = (dmin, dmax)
#     x_star = None
#     f_star = float("inf")

#     for i in range(initials):
#         x0 = [random()] * num_edges
#         # print("x0", x0)
#         out = optimize.minimize(felsenstein, x0, method="SLSQP", options={'disp':False,'maxiter':1000}, bounds=[bound]*num_edges)
#         if out.success and out.fun < f_star:
#             x_star = out.x
#             f_star = out.fun

#     return x_star, f_star

# def __sets__(seq_a, seq_b):
#     # get the msa
#     k = len(seq_a)
    
#     ## calculate the sets  
#     s_0, s_1a, s_1b, s_2, s_3 = set(), set(), set(), set(), set()
#     for idx in range(len(seq_a)):
#         c_a, c_b = seq_a[idx], seq_b[idx]
#         if c_a == c_b:
#             if c_a == 0:
#                 s_0.add(idx)
#             else:
#                 s_2.add(idx)
#         elif c_a == 0: # then c_b != 0
#             s_1b.add(idx)
#         elif c_b == 0: # then c_a != 0
#             s_1a.add(idx)
#         else:
#             s_3.add(idx)
    
#     assert len(s_0) + len(s_1a) + len(s_1b) + len(s_2) + len(s_3) == k
#     return [s_0, s_1a, s_1b, s_2, s_3,k]

# def l(x, seqa, seqb, q=0.1): # negative log-likelihood
#     s_0, s_1a, s_1b, s_2, s_3, k = __sets__(seqa, seqb)
#     s_0_len,s_1a_len,s_1b_len,s_2_len,s_3_len = len(s_0), len(s_1a), len(s_1b), len(s_2), len(s_3)

#     d_a, d_b, d_r = x
#     # print(d_a, d_b, d_r)
#     # print("p1:", log(1 - exp(-d_a)))

#     p1 = -(s_1b_len + s_0_len) * d_a + (s_1a_len + s_3_len) * log(1 - exp(-d_a))
#     p2 = -(s_1a_len + s_0_len) * d_b + (s_1b_len + s_3_len) * log(1 - exp(-d_b)) - (k - s_2_len) * d_r
#     p3 = 0.0
    
#     for i in range(s_2_len): # assuming that prob to any alpha is the same
#         # q_ialpha is the transition rate at site i from 0 -> character at node a at site i
#         # iterate over sites in s_2, get q for each alpha
#         # p3 += log(q**2 * (1 - exp(-d_a)) * (1 - exp(-d_b)) * exp(-d_r) + q*(1 - exp(-d_r)))
#         p3 += log(q * (1 - exp(-d_a)) * (1 - exp(-d_b)) * exp(-d_r) + (1 - exp(-d_r)))

#     return (p1 + p2 + p3)
    
# import sys
# sys.path.append("/Users/gillianchu/raphael/repos/problin")
# from problin_libs.sequence_lib import read_sequences
# from treeswift import *
# # k=20
# m=10 # 10
# Q = []
# datadir = "/Users/gillianchu/raphael/repos/problin/MP_inconsistent/"

# # k = 3

# # for k in [3,10,20,30,40,50,100,200,300,400,500,1000,5000]:
# for k in [5000]:

#     S = read_sequences(datadir + "seqs_m10_k" + str(k) + ".txt")
    
#     corr_d = dict()
#     wrong_d = dict()

#     num_true, total, num_incorrect, num_same = 0, 0, 0, 0
#     for i in range(k):
#         q = {j+1:1/m for j in range(m)}
#         # q = {j+1:1 for j in range(m)}
#         q[0] = 0
#         Q.append(q)

#     # D = S[1]
#     # print(D)
#     for idx, D in enumerate(S):
#     # D = dict()
#     # D['a'] = [0] * k
#     # D['b'] = [0] * k
#     # D['d'] = [0] * k
#     # print(D)
#         # true_tree = "((a:0.0360971597765934,b:3.339535381892265)e:0.0360971597765934)r;" # read_tree_newick(datadir + "m10.tre") 
#         true_tree = "((a:0.0360971597765934,b:3.339535381892265)e:0.0360971597765934,(c:0.0360971597765934,d:3.339535381892265)f:0.0360971597765934)r;"
#         # true_tree = true_tree.newick()
#         # print(D)
#         fels_true_likelihood = felsenstein(true_tree, Q, D, use_log=True) # 0.0360971597765934))
#         # x = (0.0360971597765934, 3.339535381892265, 0.0360971597765934)
#         # true_likelihood = l(x, D['a'], D['b'], q=0.1)
#         # true_norm = true_likelihood # / k 

#         false_tree = "((a:0.0360971597765934,d:3.339535381892265)e:0.0360971597765934,(c:0.0360971597765934,b:3.339535381892265)f:0.0360971597765934)r;"
#         # false_tree = "((a:0.0360971597765934,d:3.339535381892265)e:0.0360971597765934)r;"
#         # false_tree = "(a:0.07219431955,d:3.37563254167)r;" # 0.086
#         # false_tree = "((a:0.0360971597765934,d:3.339535381892265)e:0.0360971597765934,(c:0.0360971597765934,b:3.339535381892265)f:0.0360971597765934)r;" # 0.086
#         # false_tree = "((a:0.0360971597765934,b:0360971597765934)e:0360971597765934,(c:3.339535381892265,d:3.339535381892265)f:0.0360971597765934)r;"
#         # false_tree = "((a:0.0360971597765934,b:3.339535381892265)e:0.0360971597765934,(d:0.0360971597765934,c:3.339535381892265)f:0.0360971597765934)r;"
#         # T = "((a:2.5,c:3.5)e:2.5,(b:2.5,d:3.5)f:2.5)r;" # 1.0
#         # false_likelihood = felsenstein(false_tree, Q, D, use_log=True) # 0.0360971597765934))
        
#         fels_false_likelihood = felsenstein(false_tree, Q, D, use_log=True)
#         print(fels_true_likelihood > fels_false_likelihood)
#         # x = (0.07219431955, 3.37563254167, 0.0)
#         # false_likelihood = l(x, D['a'], D['d'], q=0.1)
#         # false_topo = false_likelihood # / k
#         # if fels_true_likelihood < fels_false_likelihood and true_likelihood > false_likelihood: 
#         #     print(fels_true_likelihood, fels_false_likelihood, true_likelihood, false_likelihood)

#             # print(fels_true_likelihood > fels_false_likelihood, true_likelihood > false_likelihood)
#             # print((fels_true_likelihood > fels_false_likelihood), (true_likelihood > false_likelihood))
#     #     d = dict()
#     #     for x in D:
#     #         seq = D[x]
#     #         if x not in corr_d.keys():
#     #             corr_d[x] = 0 # []
#     #             wrong_d[x] = 0 #[]
#     #         d[x] = sum([y != 0 for y in seq]) / len(seq)


#     #     print(true_norm, false_topo)
#         if fels_true_likelihood > fels_false_likelihood:
#             num_true += 1
#     #         #for x in D:
#     #             # corr_d[x].append(d[x])
#     #         # if d['b'] > d['d']:
#     #         # if d['a'] > d['c']:
#     #         #     corr_d[x] += 1
#         elif fels_true_likelihood < fels_false_likelihood: 
#             num_incorrect += 1
#     #         # for x in D:
#     #         #     wrong_d[x].append(d[x])
#     #         # if d['b'] > d['d']:
#     #         # if d['c'] < d['a']:
#     #         #     wrong_d[x] += 1
#         else:
#             num_same +=1

#         total += 1

#     print(k, "characters", num_true/total, "preferred true topology.", num_same/total, "num same")
    # print("average wrong", [(x, np.mean(np.array(wrong_d[x]))) for x in wrong_d], "over", len(wrong_d['a']))
    # print("average correct", [(x, np.mean(np.array(corr_d[x]))) for x in corr_d], "over", len(corr_d['a']))
    
    # print("average wrong", [(x, wrong_d[x]) for x in wrong_d], "over")
    # print("average correct", [(x, corr_d[x]) for x in corr_d], "over")

### Small example

# T = '[&R] (((1:0.96,3:0.96):0.14,2:1.10,4:1.10):0.0360971597765934);' # non-star tree
# msa = np.array([[1,0,2], # species 1
#                 [1,0,0], # species 2
#                 [0,1,2], # species 3
#                 [0,1,0]] # species 4
#             )

# # prob of transitioning 0 -> 0, 0 -> 1, 0 -> 2
# Q = [{0: 0, 1:0.3, 2:0.5}, # site 1
#         {0: 0, 1:0.3, 2:0.5}, # site 2
#         {0: 0, 1:0.3, 2:0.5}, # site 3
#     ]
# print("likelihood",felsenstein(T, Q, msa, use_log=True))
# T = '[&R] (((1:0.96,2:0.96):0.14,3:1.10,4:1.10):0.0360971597765934);'
# print("likelihood",felsenstein(T, Q, msa, use_log=True))
# T = '[&R] (((1:3,3:3):1,2:3,4:3):1);' # non-star tree
# print("likelihood",felsenstein(T, Q, msa, use_log=True))

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


def wrapper_felsenstein(T, Q, msa, use_log=True, initials=20, optimize_branchlengths=False):
    numsites = len(msa[next(iter(msa.keys()))])
    alphabet = dict()
    for site in range(numsites):
        alphabet[site] = Q[site].keys()

    nwkt = dendropy.Tree.get(data=T, schema="newick", rooting="force-rooted")
    num_edges = len(list(nwkt.postorder_edge_iter()))
    
    # check what the root edge length is. if none, set to 0.0
    # print("branch length of seed node", get_branchlen(nwkt.seed_node))
    root_edge_len = get_branchlen(nwkt.seed_node)
    if root_edge_len == None:
        root_edge_len = 0.0

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
            
            # print(node_likelihood.keys())
        if use_log:
            tree_likelihood = 0.0
        else:
            tree_likelihood = 1.0

        # if nwkt.is_rooted:
        #     print("tree is rooted.", nwkt)
        #     for site in range(numsites):
        #         if use_log:
        #             tree_likelihood += np.log(node_likelihood[n][site][0])
        #         else:
        #             tree_likelihood *= node_likelihood[n][site][0]
            
        # elif not nwkt.is_rooted:
            # print("tree is not rooted.", nwkt)

        # assuming that n is the seed_node (forced or otherwise)
        for site in range(numsites):
            for rootchar in node_likelihood[n][site].keys():
                prob_rootchar = node_likelihood[n][site][rootchar]
                if prob_rootchar > 0.0: 
                    q_ialpha = Q[site][rootchar]
                    if rootchar == 0:
                        if use_log:  # standard log is base e
                            tree_likelihood += np.log(math.exp(-root_edge_len)) + np.log(prob_rootchar) # + np.log(q_ialpha) 
                        else:
                            tree_likelihood *= (math.exp(-root_edge_len)) * prob_rootchar # * q_ialpha 
                    else:
                        if use_log:
                            tree_likelihood += np.log((1 - math.exp(-root_edge_len))) + np.log(q_ialpha) + np.log(prob_rootchar)
                        else:
                            tree_likelihood *= ((1 - math.exp(-root_edge_len)) * q_ialpha * prob_rootchar)
        
        return -tree_likelihood

    if optimize_branchlengths: 
        x_star = []
        dmax = -log(1/numsites)*2
        dmin = -log(1-1/numsites)/2
        bound = (dmin, dmax)
        # print("bound", bound)

        x_star = None
        f_star = float("inf")

        # x0 = []
        # for i, e in enumerate(nwkt.postorder_edge_iter()):
        #     x0.append(e.length)
        # print("x0", x0)
        # GC: add this back in later
        for i in range(initials):
            x0 = [random() * (dmax - dmin) + dmin] * num_edges
            out = optimize.minimize(felsenstein, x0, method="SLSQP", options={'disp':False,'maxiter':1000}, bounds=[bound]*num_edges)
            if out.success and out.fun < f_star:
                x_star = out.x
                f_star = out.fun
            # print(out.fun, out.x)
        for i, e in enumerate(nwkt.postorder_edge_iter()):
            e.length = x_star[i]
        return -f_star, nwkt.as_string("newick"), x_star

    else:
        x0 = []
        # same way that we put it into the tree
        for i, e in enumerate(nwkt.postorder_edge_iter()):
            x0.append(e.length)
        return -felsenstein(x0), nwkt.as_string("newick"), x0

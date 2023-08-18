#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This code tests the performance of Algorithms 3, 4 for Scenario 2 in my paper
https://arxiv.org/abs/2301.03827

From NCBI sequencing database, we obtain gene sequences of three individuals 
of Escherichia coli strain ST540 (GenBank CP007265.1, GenBank CP007390.1, 
GenBank CP007391.1) and three individuals of Escherichia coli strain ST2747 
(GenBank CP007392.1, GenBank CP007393.1, GenBank CP007394.1). Genes that do not
appear in all sequences of the same strain have been removed. Genes that have 
multiple copies are also removed.

The gene names are first transformed to different numbers. Then we run the 
algorithm to find non-transposons, quasi-transposons, and proper-transposons.
"""

"""
Consider genes 1, 2,..., n - 2.
The input is m integer sequences, representing gene sequences. Each integer
only appears once in the same sequence. We need to find the longest common 
subsequence. 
Notice that all sequences and subsequences are cyclic, meaning that they can
be rotated to have different equivalent forms. 
For example, [3, 1, 2], [1, 2, 3], [2, 3, 1] are all equivalent forms of the
same cyclic sequence. 
Algorithm 3 finds one longest common subsequence. 

The longest common subsequence might not be unique. Algorithm 4 determines 
whether each gene appears in all/some/none of the longest common subsequences.
They are called non-transposons, quasi-transposons, proper-transposons 
accordingly.

@author: yuewang
"""
def Algorithm1(sequences): # duplicate of Algorithm 1 in the paper
# see the file NewScenario1.py for more explanations
    m = len(sequences) # number of sequences
    n = len(sequences[0]) # number of genes

    standard_subseq = [[[] for i in range(n)] for j in range(n)]
    # standard_subseq[i][j] (i < j) is the subsequence of the first gene sequence
    # that only consists of gene i and gene j.
    # it is the standard that we will use to compare to other sequences.
    # standard_subseq[i][j] for i >= j is not used

    # construct standard_subseq
    for curr_gene in sequences[0]: # for each gene copy in the first sequence
        for i in range(curr_gene): # for i smaller than curr_gene
            standard_subseq[i][curr_gene].append(curr_gene)
        # add the current gene copy to the subsequence of i and curr_gene
        for i in range(curr_gene+1, n): # for i larger than curr_gene
            standard_subseq[curr_gene][i].append(curr_gene)
        # add the current gene copy to the subsequence of i and curr_gene

    graph_edge = [[False for i in range(n)] for j in range(n)]
    # graph_edge[i][j] means whether i appears before j for all sequences.
    # if yes, there is an edge from i to j in the auxiliary graph.
    # here graph_edge is initialized. we will update it later.
    for i in range(n):
        for j in range(i+1,n): # for each pair i, j
            graph_edge[standard_subseq[i][j][0]][standard_subseq[i][j][1]] = True
    # check whether i appears before or after j in the first sequence, and add it 
    # to graph_edge

    # for each gene sequence (except the first one), construct the subsequences
    # and compare them with standard_subseq to see if the subsequence is the same
    # in the current gene sequence and the first gene sequence
    for counter in range(1, m): # for each gene sequence except the first one
        curr_subseq = [[[] for i in range(n)] for j in range(n)]
        # curr_subseq[i][j] (i < j) is the subsequence of the current gene sequence
        # that only consists of gene i and gene j.
        # curr_subseq[i][j] for i >= j is not used
        
        # construct curr_subseq
        for curr_gene in sequences[counter]: # for each gene copy in the current sequence
            if curr_gene >= n: 
                # if the current gene does not appear in the first sequence,
                # skip this gene, since it cannot form a subsequence that is 
                # the same for all sequences
                continue
            for i in range(curr_gene): # for i smaller than curr_gene
                curr_subseq[i][curr_gene].append(curr_gene)
            # add the current gene copy to the subsequence of i and curr_gene
            for i in range(curr_gene+1, n): # for i larger than curr_gene
                curr_subseq[curr_gene][i].append(curr_gene)
            # add the current gene copy to the subsequence of i and curr_gene
        
        # compare curr_subseq and standard_subseq, and update the result in
        # subseq_compare
        for i in range(n):
            for j in range(i+1, n): # for each gene pair j < k
                if curr_subseq[i][j] != standard_subseq[i][j]:
                    graph_edge[i][j] = False 
                    graph_edge[j][i] = False
                # if the order of i and j from the current sequence is not
                # same as the order from the first sequence, update
                # subseq_compare[i][j]

    # now we have a directed acyclic graph (DAG) with edges in graph_edge.
    # next, we perform a topological sorting for vertices 0,1,...,n-1.
    # after this sorting, there is no directed edge from a back vertex to a
    # front vertex

    degree = [0] * n 
    # the degree (number of edges starting from this vertex) of each vertex    
    for i in range(n):
        for j in range(n):
            if graph_edge[i][j] == True: # if there is an edge from i to j
                degree[i] += 1 # update degree         
    # now perform topological sort. the idea is to put one vertex with degree 0
    # to the end, and delete edges that point to this vertex, and update degree.
    # repeat this procedure to pop out vertices one by one
    topological_sort = [] # store the result of topological sorting
    vertices_0_degree = [] # the list of vertices with 0 degree
    for i in range(n):
        if degree[i] == 0:
            vertices_0_degree.append(i) # find all vertices with 0 degree
    # initially, since n-1 is at the end of each sequence, every other vertex has
    # an edge pointing to n-1. thus initially vertices_0_degree has only n-1

    while len(vertices_0_degree) >0: # while there is a vertex with 0 degree
        curr_vertex = vertices_0_degree.pop() # choose a vertex with 0 degree
        topological_sort.append(curr_vertex) # add it to the sorting result
        for i in range(n):
            if graph_edge[i][curr_vertex] == True:
                degree[i] -= 1 
                # delete edges pointing to curr_vertex and update degree
                if degree[i] == 0:
                    vertices_0_degree.append(i)
                    # update vertices_0_degree
    if len(topological_sort) < n: # if topological sorting terminates before all
    # n vertices are popped out, then the graph has a directed cycle.
        print('Topological sorting error!')
    topological_sort = [topological_sort[i] for i in range(n-1, -1, -1)]
    # final result of topological sorting

    # with the result of topological sorting, we can calculate some quantities:
    # longest_path_to_tail: the length of the longest path from vertex i to n - 1
    # in the graph
    # next_vertex_to_tail: the next vertex on the longest path from vertex i to 
    # n - 1 in the graph. for n - 1, this quantity is -1.
    longest_path_to_tail = [0] * n
    next_vertex_to_tail = [-1] * n
    for i in range(n-2, -1, -1):
        # by definition, longest_path_to_tail[n-1] == 0, 
        # next_vertex_to_tail[n-1] == -1
        # notice that topological_sort[n-1] == n - 1
        curr_vertex = topological_sort[i] # we go back in the result of 
        # topological sorting
        for j in range(i+1,n): # for each vertex that might be the next vertex 
        # of curr_vertex
            next_vertex = topological_sort[j]
            if graph_edge[curr_vertex][next_vertex] == True: # there is an edge
                if longest_path_to_tail[next_vertex] + 1 > \
                    longest_path_to_tail[curr_vertex]:
                        # if the length of path from curr_vertex through 
                        # next_vertex to n - 1 is longer
                    longest_path_to_tail[curr_vertex] = \
                        longest_path_to_tail[next_vertex] + 1
                    next_vertex_to_tail[curr_vertex] = next_vertex
                    # update longest_path_to_tail and next_vertex_to_tail
                    
    # with the result of topological sorting, we can calculate some quantities:
    # longest_path_from_head: the length of the longest path from 0 to vertex i
    # in the graph
    # previous_vertex_from_head: the previous vertex on the longest path from 0 to
    # vertex i in the graph. for 0, this quantity is -1.
    longest_path_from_head = [0] * n
    previous_vertex_from_head = [-1] * n
    for i in range(1, n):
        # by definition, longest_path_from_head[0] == 0, 
        # previous_vertex_from_head[0] == -1
        # notice that topological_sort[0] == 0
        curr_vertex = topological_sort[i] # we go forth in the result of 
        # topological sorting
        for j in range(i): # for each vertex that might be the previous vertex 
        # of curr_vertex
            previous_vertex = topological_sort[j]
            if graph_edge[previous_vertex][curr_vertex] == True: # there is an edge
                if longest_path_from_head[previous_vertex] + 1 > \
                    longest_path_from_head[curr_vertex]:
                        # if the length of path from 0 through 
                        # previous_vertex to curr_vertex is longer
                    longest_path_from_head[curr_vertex] = \
                        longest_path_from_head[previous_vertex] + 1
                    previous_vertex_from_head[curr_vertex] = previous_vertex
                    # update longest_path_from_head and previous_vertex_from_head

    # now construct one longest path from 0 to n - 1
    longest_path_0 = [False] * n # whetehr each vertex is in the longest path. 
    # this corresponds to a longest common subsequence for the input.
    longest_path_0[0] = True
    curr_vertex = 0 # most recent vertex found in generating the longest path
    while next_vertex_to_tail[curr_vertex] != -1: # not reaching n - 1
        curr_vertex = next_vertex_to_tail[curr_vertex]
        longest_path_0[curr_vertex] = True
    return longest_path_to_tail, next_vertex_to_tail, longest_path_from_head, \
        previous_vertex_from_head, longest_path_0
        

def Algorithm2(sequences): # duplicate of Algorithm 2 in the paper
# see the file NewScenario1.py for more explanations
    longest_path_to_tail, next_vertex_to_tail, longest_path_from_head, \
        previous_vertex_from_head, longest_path_0 = Algorithm1(sequences)
    transposon = [1] * n # whether each gene is a proper/quasi/non-transposon
    # 1 is non-transposon; 0 is quasi-transposon; -1 is proper-transposon
    # non-transposon means this gene appears in all longest common subsequences.
    # quasi-transposon means this gene appears in some but not all longest common
    # subsequences.
    # proper-transposon means this gene does not appear in any longest common
    # subsequences.

    existed_path_length = [-1] * n 
    # if the current gene is not a propere-transposon, record if the 
    # longest_path_to_tail of the current gene has been seen. If yes, both genes 
    # are quasi-transposons. If longest_path_to_tail of one gene only appears once,
    # this gene is a non-transposon.
    # existed_path_length[i] == -1 if we have not seen a gene with 
    # longest_path_to_tail == i. Otherwise, existed_path_length[i] = gene name.
    for i in range(1, n-1):
        if longest_path_to_tail[i] + longest_path_from_head[i] \
            < longest_path_to_tail[0]: 
                # the longest path from 0 through i to n - 1 is shorter than the 
                # longest path from 0 to n - 1. thus gene i is a proper-transposon
            transposon[i] = -1
        else:
            if existed_path_length[longest_path_to_tail[i]] == -1:
                existed_path_length[longest_path_to_tail[i]] = i
                # we have not seen existed_path_length[longest_path_to_tail[i]]
                # with this longest_path_to_tail. Record the current gene name
            else:
                transposon[i] = 0
                transposon[existed_path_length[longest_path_to_tail[i]]] = 0
                # we have seen gene existed_path_length[longest_path_to_tail[i]]
                # with this longest_path_to_tail, so both genes are 
                # quasi-transposons           

    # a gene that is not determined to be proper/quasi-transposon here
    # is a non-transposon, so that transposon[i] = 1.
    # in Lemma 1 of the paper, we show that for each quasi-transposon or 
    # non-transposon, its position in any LCS is fixed, determined by 
    # longest_path_to_tail. Thus after excluding proper-transposons, if two genes 
    # have the same longest_path_to_tail, they are both quasi-transposable. If a 
    # gene has a unique longest_path_to_tail, it is a non-transposon.

    return transposon


def cut_circular_seq(i, sequences): # cut the sequences at i and add auxiliary 
    #0 and n-1
    m = len(sequences) # number of sequences
    n = len(sequences[0]) + 2 # number of genes, including auxiliary 0 and 
    # n - 1
    cut_sequence_i = [[-1] * n for j in range(m)] # initialize the cut sequences
    for j in range(m): # for each sequence
        cut_sequence_i[j][0] = 0
        cut_sequence_i[j][n-1] = n - 1 # add auxiliary genes
        location_i = sequences[j].index(i) # where i appears in sequence j
        for k in range(location_i, n-2):
            cut_sequence_i[j][k-location_i+1] = sequences[j][k]
        for k in range(0, location_i):
            cut_sequence_i[j][n-location_i+k-1] = sequences[j][k]
            # fill in genes to proper locations
    return cut_sequence_i

#########################

dataset = 2
# dataset=1 means all three variants of ST540;
# 2 means last two variants of ST540
# 3 means all three variants of ST2747

# read gene sequences
if dataset <= 2:
    with open('no_repeats_CP007265.1.txt') as f:
        lines = f.readlines()
    gene_name_sequence_1 = lines
    with open('no_repeats_CP007390.1.txt') as f:
        lines = f.readlines()
    gene_name_sequence_2 = lines
    with open('no_repeats_CP007391.1.txt') as f:
        lines = f.readlines()
    gene_name_sequence_3 = lines
else:
    with open('no_repeats_CP007392.1.txt') as f:
        lines = f.readlines()
    gene_name_sequence_1 = lines
    with open('no_repeats_CP007393.1.txt') as f:
        lines = f.readlines()
    gene_name_sequence_2 = lines
    with open('no_repeats_CP007394.1.txt') as f:
        lines = f.readlines()
    gene_name_sequence_3 = lines
  
name_dict = {} # dictionary of mapping gene names to numbers (nicknames)
inverse_name_dict = {}
# dictionary of mapping gene numbers (nicknames) to names
counter = 1
gene_number_sequence_1 = [] 
for gene_name in gene_name_sequence_1:
    if gene_name not in name_dict: 
        # for a new gene, give it a number as its nickname   
        name_dict[gene_name] = counter
        inverse_name_dict[counter] = gene_name
        gene_number_sequence_1.append(counter)
        counter += 1  
    else:
        # for a encountered gene, directly add its nickname into the sequence
        gene_number_sequence_1.append(name_dict[gene_name])
gene_number_sequence_2 = []
for gene_name in gene_name_sequence_2:
    gene_number_sequence_2.append(name_dict[gene_name])
gene_number_sequence_3 = []
for gene_name in gene_name_sequence_3:
    gene_number_sequence_3.append(name_dict[gene_name])
    

if dataset != 2:
    sequences = [gene_number_sequence_1, gene_number_sequence_2, 
                 gene_number_sequence_3]
else:
    sequences = [gene_number_sequence_2, gene_number_sequence_3]
# sequences is the input of m circular sequences, where each sequence
# has n-2 genes 1 to n-2 (no replicated genes)
# sequences can rotate

m = len(sequences) # number of sequences
n = len(sequences[0]) + 2 # length of sequence (adding auxiliary 0 and n-1)

chosen_cut = [False] * n # chosen_cut[i] == True if gene i has been chosen and cut

length_lcs = [-1] * n # length_lcs[i] is the length of the longest common subsequence
# when the circular sequences are cut at i

cut_sequence_i = cut_circular_seq(1, sequences) # cut at gene 1

_, _, _, _, longest_path = Algorithm1(cut_sequence_i) # find a longest common 
# subsequence in the linear case

max_path_length = sum(longest_path) # the length of the current longest common
# subsequence
chosen_cut[1] = True # we have cut at gene 1
length_lcs[1] = max_path_length

while sum(chosen_cut) + sum(longest_path) < n + 1: # total number of cut is k + 1
# where n - k is the length of the longest common subsequence
    for i in range(n):
        if chosen_cut[i] == False and longest_path[i] == False: 
        # find a gene not in longest_path which has not been cut
            break
    cut_sequence_i = cut_circular_seq(i, sequences) # cut at i
    _, _, _, _, longest_path_i = Algorithm1(cut_sequence_i)
    # find the longest common subsequence for cut sequences
    path_length_i = sum(longest_path_i)
    chosen_cut[i] = True
    length_lcs[i] = path_length_i
    if path_length_i > max_path_length: # update the longest common subsequence
        max_path_length = path_length_i
        longest_path = longest_path_i
# Algorithm 3 terminates here. longest_path is a longest common subsequence
# for the circular sequences

# Algorithm 4 starts here. try to determine whether each gene is in 
# all/some/none of longest common subsequences.

transposon = [1] * n # whether each gene is a proper/quasi/non-transposon
# 1 is non-transposon; 0 is quasi-transposon; -1 is proper-transposon
# non-transposon means this gene appears in all longest common subsequences.
# quasi-transposon means this gene appears in some but not all longest common
# subsequences.
# proper-transposon means this gene does not appear in any longest common
# subsequences.

for i in range(n):
    if longest_path[i] == False: # not in longest_path
        if length_lcs[i] < max_path_length: # the longest common subsequence
        # for sequences cut at i is shorter than the actual lcs
            transposon[i] = -1 # determine proper-transposons
        else:
            transposon[i] = 0 # genes not in longest_path cannot be 
            # non-transposons
            cut_sequence_i = cut_circular_seq(i, sequences)
            cut_transposon = Algorithm2(cut_sequence_i)
            # cut at i and check if some genes in longest_path are determined
            # as quasi-transposons.
            for j in range(n): # determine quasi-transposons
                if longest_path[j] == True and cut_transposon[j] <= 0:
                    transposon[j] = 0
# by Lemma 2 in the paper, every quasi-transposon can be found in this way.
# thus the remaining genes are non-transposons.

# output the final results of Algorithm 4. Notice that 0 and n - 1 are
# auxiliary, not actual genes.
non_transposons = []
quasi_transposons = []
proper_transposons=[]
for i in range(1, n-1): # output
    if transposon[i] == 1:
        non_transposons.append(inverse_name_dict[i])
    elif transposon[i] == 0:
        quasi_transposons.append(inverse_name_dict[i])
    else:
        proper_transposons.append(inverse_name_dict[i])    
            
print('non-transposons: ', non_transposons)        
print('quasi-transposons: ', quasi_transposons)        
print('proper-transposons: ', proper_transposons)  

print('number of non-transposons: ', len(non_transposons))        
print('number of quasi-transposons: ', len(quasi_transposons))  
print('number of proper-transposons: ', len(proper_transposons))  

"""
For dataset 1:
    number of non-transposons:  389
    number of quasi-transposons:  50
    number of proper-transposons:  129
For dataset 2:
    quasi-transposons:  ['hpaC\n', 'iraD\n', 'fbpC\n', 'psiB\n']
    number of non-transposons:  564
    number of quasi-transposons:  4
    number of proper-transposons:  0
For dataset 3:
    number of non-transposons:  573
    number of quasi-transposons:  0
    number of proper-transposons:  0
"""   

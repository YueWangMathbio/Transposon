# -*- coding: utf-8 -*-
"""
This code implements Algorithm 5 for Scenario 3 in my paper
https://arxiv.org/abs/2301.03827

Consider genes 0,1,2,...,n-1.
The input is m integer sequences, representing gene sequences. Each integer
can appear multiple times in the same sequence. This means each sequence 
contains different numbers of copies of genes 0,1,2,...,n-1. We need to find 
the longest common subsequence. Here we only consider common subsequences that 
consist of all or none copies of the same gene, and the subsequence length is 
calculated by genes, not gene copies.
After finding the longest common subsequence, we output that genes in this 
longest common subsequence are transposons, and other genes are 
non-transposons.

The basic idea is to construct a graph and transform the problem of finding the
longest common subsequence into finding a certain structure in the graph. See 
below for more details.

This algorithm does not always produce the correct result. 
See below for a counterexample.
@author: yuewang
"""

sequences = [[0, 1, 2, 3, 3, 2, 1], 
             [1, 2, 3, 3, 2, 1, 0]]
# sequences is the input, consisting of of m integer sequences.
# for each sequence, each integer i represents a copy of gene i.
# the same gene i might have multiple copies in each sequence.

"""
The following example fails this algorithm.
The correct outcome should be 7, 8, 9, 0. 
The outcome of this algorithm is 2, 4, 6.

sequences = [[7, 8, 9, 0, 1, 1, 2, 3, 3, 4, 5, 5, 6],
   [1, 2, 1, 3, 4, 3, 5, 6, 5, 7, 8, 9, 0]]
"""

"""
The first half is to construct an auxiliary graph, where each vertex is a gene, 
and there is an undirected edge between gene i and gene j if and only if the 
subsequence consisting of i and j is the same for all gene sequences.
For example, if we have two sequences [0, 1, 2] and [0, 2, 1], the subsequences
of 0 and 1 are [0, 1] and [0, 1], and there is an edge between 0 and 1. The 
subsequences of 1 and 2 are [1, 2] and [2, 1], and there is no edge between 0
and 1.
See below for how this auxiliary graph is used to find the longest common 
subsequence.
"""

n = max(sequences[0]) + 1 # number of different genes in the first sequence
m = len(sequences) # number of sequences in the input

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

subseq_compare = [[True for i in range(n)] for j in range(n)]
# subseq_compare[i][j] (i < j) means whether the subsequence consisting of 
# gene i and gene j is the same for all sequences.
# subseq_compare[i][j] for i >= j is not used

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
                subseq_compare[i][j] = False 
            # if the subsequence of i and j from the current sequence is not
            # same as the subsequence from the first sequence, update
            # subseq_compare[i][j]

# with the result in subseq_compare, we can construct an undirected graph,
# where each vertex is a gene, and there is an edge between gene i and gene j
# if and only if the subsequence of i and j is the same for all gene sequences,
# that is, subseq_compare[i][j] == True
edges = [[] for i in range(n)] # if there is an edge between i and j,
# add i to edges[j], and add j to edges[i]

# construct edges
for i in range(n):
    for j in range(i+1, n): # for each gene pair i < j
        if subseq_compare[i][j] == True: 
            # if the subsequence of i and j is the same for all gene sequences
            edges[i].append(j)
            edges[j].append(i)
            # add i to edges[j], and add j to edges[i]
            
"""
We have constructed an auxiliary graph. According to the paper, finding the 
longest common subsequence is equivalent to finding the maximum clique of 
this graph. Here a clique is a subgraph that every two vertices have an edge
that links them. The maximum clique is the clique that has the largest number
of vertices.

Since the problem of finding the maximum clique is NP-hard, we adopt a 
heuristic algorithm that repeatedly delete the vertex with the minimal 
degree, until there is a clique.
"""

degree = [len(edges[i]) for i in range(n)]
# degree[i] is the degree of vertex i in this graph, namely the number of 
# edges linking i

gene_number = n # number of genes that have not been deleted
deleted = [False] * n # whether gene i has been deleted
while gene_number > 0: # when there is at least one gene left
    min_degree = float('inf') # the minimum degree have seen
    min_degree_index = -1 # the index of the minimum degree
    for i in range(n):
        if deleted[i] == False:
            if degree[i] < min_degree:
                min_degree = degree[i]
                min_degree_index = i
    if min_degree == gene_number - 1: # we find a clique
        break
    else:
        # min_degree_index is the gene that will be deleted    
        for vertex in edges[min_degree_index]:
            if deleted[vertex] == False:
                degree[vertex] -= 1
        deleted[min_degree_index] = True
        gene_number -= 1
# a single vertex is a clique. we can always break the while loop with a clique
# we have found a clique (genes that are not deleted), and we hope that it is 
# the maximum clique

# output the results: genes in that clique should be non-transposons,
# other genes should be transposons
transposons = [] # list of possible transposons
non_transposons = [] # list of possible non-transposons
for i in range(n):
    if deleted[i] == True:
        transposons.append(i)
    else:
        non_transposons.append(i)

print('transposons: ', transposons)
print('non-transposons: ', non_transposons)
                
        
    
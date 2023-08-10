# -*- coding: utf-8 -*-
"""
This code implements Algorithm 6 for Scenario 4 in my paper
https://arxiv.org/abs/2301.03827

Consider genes 0,1,2,...,n-1.
The input is m integer sequences, representing gene sequences. Each integer
can appear multiple times in the same sequence. This means each sequence 
contains different numbers of copies of genes 0,1,2,...,n-1. We need to find 
the longest common subsequence. Here we only consider common subsequences that 
consist of all or none copies of the same gene, and the subsequence length is 
calculated by genes, not gene copies.
Notice that all sequences and subsequences are cyclic, meaning that they can
be rotated to have different equivalent forms. 
For example, [0, 1, 2], [1, 2, 0], [2, 0, 1] are all equivalent forms of the
same cyclic sequence.
After finding the longest common subsequence, we output that genes in this 
longest common subsequence are transposons, and other genes are 
non-transposons.

The basic idea is to construct a 3-uniform hypergraph and transform the problem of 
finding thelongest common subsequence into finding a certain structure in the 
3-uniform hypergraph. See below for more details.

For a 3-uniform hypergraph, there are some vertices, and certain vertex triples are
linked by 3-hyperedges.

This algorithm does not always produce the correct result. 
See below for a counterexample.
@author: yuewang
"""

sequences = [[0, 1, 2, 3, 4, 1],
   [1, 1, 2, 0, 3, 4]]
# sequences is the input, consisting of of m integer cyclic sequences.
# for each sequence, each integer i represents a copy of gene i.
# the same gene i might have multiple copies in each sequence.


"""
The following example fails this algorithm.
The correct outcome should be 7, 8, 9, 0.
The outcome of this algorithm is 2, 4, 6.
sequences = [[1, 2, 7, 3, 4, 8, 9, 5, 6, 0],
   [2, 1, 0, 4, 3, 7, 8, 6, 5, 9],
   [1, 2, 9, 3, 4, 0, 7, 5, 6, 8],
   [2, 1, 8, 4, 3, 9, 0, 6, 5, 7]]
"""

"""
The first half is to construct an auxiliary 3-uniform hypergraph, where each vertex is 
a gene, and there is an undirected 3-hyperedge among genes i, j, k if and only 
if the subsequence consisting of i, j, k is the same for all gene sequences in 
cyclic sense.
For example, if we have two sequences [0, 1, 2, 3] and [3, 0, 2, 1], the 
subsequences of 0, 1, 3 are [0, 1, 3] and [3, 0, 1], being the same in cyclic 
sense. Thus there is a 3-hyperedge among 0, 1, 3. The subsequences of 0, 1, 2 
are [0, 1, 2] and [0, 2, 1], being different in cyclic sense. Thus there is no 
edge among 0, 1, 2.
See below for how this auxiliary 3-uniform hypergraph is used to find the longest 
common subsequence.
"""

def shift_sort(arr):
    # for a cyclic integer sequence arr, from all its equivalent forms,
    # find the smallest alphabetically
    min_form = arr[:]
    l = len(arr)
    for i in range(1, l):
        temp_arr = arr[i:] + arr[:i] # for each equivalent form
        if temp_arr < min_form:
            min_form = temp_arr
    return min_form

n = max(sequences[0]) + 1 # number of different genes in the first sequence
m = len(sequences) # number of sequences in the input

standard_subseq = [[[[] for i in range(n)] for j in range(n)] for k in range(n)]
# standard_subseq[i][j][k] (i < j < k) is the subsequence of the first gene 
# sequence that only consists of genes i, j, k.
# it is the standard that we will use to compare to other sequences.
# if i < j < k does not hold, standard_subseq[i][j][k] is not used


# construct standard_subseq
for curr_gene in sequences[0]: # for each gene copy in the first sequence
    for i in range(curr_gene): # for i < j < curr_gene
        for j in range(i+1, curr_gene):
            standard_subseq[i][j][curr_gene].append(curr_gene)
    # add the current gene copy to the subsequence of i, j, curr_gene
    for i in range(curr_gene): # for i < curr_gene < j
        for j in range(curr_gene+1, n):
            standard_subseq[i][curr_gene][j].append(curr_gene)
    # add the current gene copy to the subsequence of i, j, curr_gene
    for i in range(curr_gene+1, n): # for curr_gene < i < j
        for j in range(i+1, n):
            standard_subseq[curr_gene][i][j].append(curr_gene)
    # add the current gene copy to the subsequence of i, j, curr_gene

# since the subsequence is cyclic, we use the function shift_sort to transform
# it into the smallest form alphabetically.
for i in range(n):
    for j in range(i+1, n):
        for k in range(j+1, n): # for each gene triple i < j < k
            standard_subseq[i][j][k] = shift_sort(standard_subseq[i][j][k])

subseq_compare = [[[True for i in range(n)] for j in range(n)] for k in range(n)]
# subseq_compare[i][j][k (i < j < k) means whether the subsequence consisting  
# of genes i, j, k is the same for all sequences.
# if i < j < k does not hold, subseq_compare[i][j][k] is not used

# for each gene sequence (except the first one), construct the subsequences
# and compare them with standard_subseq to see if the subsequence is the same
# in the current gene sequence and the first gene sequence
for counter in range(1, m): # for each gene sequence except the first one
    curr_subseq = [[[[] for i in range(n)] for j in range(n)] for k in range(n)]
    # curr_subseq[i][j][k] (i < j < k) is the subsequence of the current gene sequence
    # that only consists of genes i, j, k.
    # if i < j < k does not hold, curr_subseq[i][j][k] is not used
    
    # construct curr_subseq
    for curr_gene in sequences[counter]: # for each gene copy in the current sequence
        if curr_gene >= n: 
            # if the current gene does not appear in the first sequence,
            # skip this gene, since it cannot form a subsequence that is 
            # the same for all sequences
            continue
        for i in range(curr_gene): # for i < j < curr_gene
            for j in range(i+1, curr_gene):
                curr_subseq[i][j][curr_gene].append(curr_gene)
        # add the current gene copy to the subsequence of i, j, curr_gene
        for i in range(curr_gene): # for i < curr_gene < j
            for j in range(curr_gene+1, n):
                curr_subseq[i][curr_gene][j].append(curr_gene)
        # add the current gene copy to the subsequence of i, j, curr_gene
        for i in range(curr_gene+1, n): # for curr_gene < i < j
            for j in range(i+1, n):
                curr_subseq[curr_gene][i][j].append(curr_gene)
        # add the current gene copy to the subsequence of i, j, curr_gene
        
    # since the subsequence is cyclic, we use the function shift_sort to transform
    # it into the smallest form alphabetically.
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n): # for each gene triple i < j < k
                curr_subseq[i][j][k] = shift_sort(curr_subseq[i][j][k])
            
    # compare curr_subseq and standard_subseq, and update the result in
    # subseq_compare
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n): # for each gene triple i < j < k
                if curr_subseq[i][j][k] != standard_subseq[i][j][k]:
                    subseq_compare[i][j][k] = False 
                # if the subsequence of j and k from the current sequence is not
                # same as the subsequence from the first sequence, update
                # subseq_compare[j][k]

# with the result in subseq_compare, we can construct an undirected 3-uniform hypergraph,
# where each vertex is a gene, and there is a 3-hyperedge among genes i, j, k
# if and only if the subsequence of i, j, k is the same for all gene sequences,
# that is, subseq_compare[i][j][k] == True
            
"""
We have constructed an auxiliary 3-uniform hypergraph. According to the paper, finding 
the longest common subsequence can be transformed into finding the maximum 
clique of this graph. Here a clique is a subgraph that every three vertices 
have a 3-hyperedge that links them. The maximum clique is the clique that has 
the largest number of vertices.

Since the problem of finding the maximum clique in 3-uniform hypergraph is NP-hard, we 
adopt a heuristic algorithm that repeatedly delete the vertex with the minimal 
degree, until there is a clique.
"""

degree = [0] * n
# degree[i] is the degree of vertex i in this graph, namely the number of 
# 3-hyperedges linking i
for i in range(n):
    for j in range(i+1, n):
        for k in range(j+1, n): # for each gene triple i, j, k
            if subseq_compare[i][j][k] == True: # if they form a 3-hyperedge
                degree[i] += 1 # update degree
                degree[j] += 1
                degree[k] += 1
                
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
    if min_degree == (gene_number - 1) * (gene_number - 2) //2: 
        # we find a clique
        break
    else:
        # min_degree_index is the gene that will be deleted   
        for i in range(min_degree_index): # for i < j < min_degree_index
            for j in range(i+1, min_degree_index):
                if subseq_compare[i][j][min_degree_index] == True \
                    and deleted[i] == False and deleted[j]==False:
                    # if there is a 3-hyperedge that has not been deleted
                    degree[i] -= 1
                    degree[j] -= 1
                    degree[min_degree_index] -= 1
        # remove the 3-hyperedge of i, j, min_degree_index
        for i in range(min_degree_index): # for i < j < min_degree_index
            for j in range(min_degree_index+1, n):
                if subseq_compare[i][min_degree_index][j] == True \
                    and deleted[i] == False and deleted[j]==False:
                    # if there is a 3-hyperedge that has not been deleted
                    degree[i] -= 1
                    degree[j] -= 1
                    degree[min_degree_index] -= 1
        # remove the 3-hyperedge of i, j, min_degree_index
        for i in range(min_degree_index+1, n): # for i < j < min_degree_index
            for j in range(i+1, n):
                if subseq_compare[min_degree_index][i][j] == True \
                    and deleted[i] == False and deleted[j]==False:
                    # if there is a 3-hyperedge that has not been deleted
                    degree[i] -= 1
                    degree[j] -= 1
                    degree[min_degree_index] -= 1
        # remove the 3-hyperedge of i, j, min_degree_index        
        deleted[min_degree_index] = True # mark this gene as deleted
        gene_number -= 1 # update the number of remaining genes
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
                
        
    
# Transposon

This repository contains Python code for all algorithms in my paper https://arxiv.org/abs/2301.03827.
Besides, it contains all the real world data that I use to apply the algorithms 

NewScenario1.py implements Algorithms 1,2 for Scenario 1.
Consider m gene sequences, where each sequence is a permutation of 0, 1, 2,..., n - 1. Algorithm 1 finds one longest common subsequence. 
The longest common subsequence might not be unique. Algorithm 2 determines whether each gene appears in all/some/none of the longest common subsequences. They are called non-transposons, quasi-transposons, proper-transposons accordingly.

NewScenario2.py implements Algorithms 3,4 for Scenario 2.
Consider m gene sequences, where each sequence is a permutation of 0, 1, 2,..., n - 1. Algorithm 3 finds one longest common subsequence. Here all sequences and subsequences are cyclic, meaning that they can be rotated to have different equivalent forms. 
The longest common subsequence might not be unique. Algorithm 4 determines whether each gene appears in all/some/none of the longest common subsequences. They are called non-transposons, quasi-transposons, proper-transposons accordingly.

NewScenario3.py implements Algorithm 5 for Scenario 3.
Consider genes 0,1,2,...,n-1. The input is m integer sequences, representing gene sequences. Each integer can appear multiple times in the same sequence. This means each sequence contains different numbers of copies of genes 0,1,2,...,n-1. Algorithm 5 finds the longest common subsequence. Here we only consider common subsequences that consist of all or none copies of the same gene, and the subsequence length is calculated by genes, not gene copies. After finding the longest common subsequence, we output that genes in this longest common subsequence are transposons, and other genes are non-transposons. 
The basic idea is to construct a graph and transform the problem of finding the longest common subsequence into finding a certain structure in the graph.
This algorithm does not always produce the correct result. 

NewScenario4.py implements Algorithm 6 for Scenario 4.
Consider genes 0,1,2,...,n-1. The input is m integer sequences, representing gene sequences. Each integer can appear multiple times in the same sequence. This means each sequence contains different numbers of copies of genes 0,1,2,...,n-1. Algorithm 6 finds the longest common subsequence. Here we only consider common subsequences that consist of all or none copies of the same gene, and the subsequence length is calculated by genes, not gene copies. Notice that all sequences and subsequences are cyclic, meaning that they can be rotated to have different equivalent forms. For example, [0, 1, 2], [1, 2, 0], [2, 0, 1] are all equivalent forms of the same cyclic sequence. After finding the longest common subsequence, we output that genes in this longest common subsequence are transposons, and other genes are non-transposons.
The basic idea is to construct a 3-uniform hypergraph and transform the problem of finding thelongest common subsequence into finding a certain structure in the 3-uniform hypergraph.
For a 3-uniform hypergraph, there are some vertices, and certain vertex triples are linked by 3-hyperedges.
This algorithm does not always produce the correct result. 

Scenario3 random test.py tests the performance of Algorithm 5 on random graphs. We randomly generate graphs, and check whether Algorithm 5 can produce the correct result (compared to a brute-force searching algorithm). When generating the random graph, we test different vertex numbers and probabilities to have an edge. The success rates are outputed.

Scenario4 random test.py tests the performance of Algorithm 6 on random graphs. We randomly generate graphs, and check whether Algorithm 6 can produce the correct result (compared to a brute-force searching algorithm). When generating the random graph, we test different vertex numbers and probabilities to have an edge. The success rates are outputed.

NewScenario1 test.py runs Algorithms 1,2 on the gene sequence data (no_repeats_CP007xxx.1.txt) to find transposons.

NewScenario2 test.py runs Algorithms 3,4 on the gene sequence data (no_repeats_CP007xxx.1.txt) to find transposons.

NewScenario3 test.py runs Algorithm 5 on the gene sequence data (with_repeats_CP007xxx.1.txt) to find transposons.

NewScenario4 test.py runs Algorithm 6 on the gene sequence data (with_repeats_CP007xxx.1.txt) to find transposons.


About gene sequence data:
From NCBI sequencing database, we obtain gene sequences of three individuals of Escherichia coli strain ST540 (GenBank CP007265.1, GenBank CP007390.1, GenBank CP007391.1) and three individuals of Escherichia coli strain ST2747 (GenBank CP007392.1, GenBank CP007393.1, GenBank CP007394.1). We remove genes that do not appear in all sequences of the same strain, and record the remaining gene sequence in files with_repeats_CP007xxx.1.txt for all six gene sequence datasets. These six files are used in NewScenario1 test.py and NewScenario2 test.py. Then we also remove genes that have multiple copies, and record the remaining gene sequence in files no_repeats_CP007xxx.1.txt for all six gene sequence datasets. These six files are used in NewScenario3 test.py and NewScenario4 test.py.
Each file contains a list of gene names in the gene sequence.



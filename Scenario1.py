#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Given some linear gene sequences, where each gene has only one copy 
(meaning that each sequence is a permutation of 0 to n-1),this code 
determines whether each gene is a proper/quasi/non-transposon. 
This code implements Algorithms 1,2 for Scenario 1 in my paper
https://arxiv.org/abs/1506.02424

@author: Yue Wang
"""


def fp(G,i): # calculate f+(i) for gene i, 
# where G is the auxiliary graph
    n=len(G)
    if i==n-1:
        return 0
    else:
        return 1+fp(G,hp(G,i))

def hp(G,i): # calculate h+(i) for gene i, 
# where G is the auxiliary graph
    n=len(G)
    if i==n-1:
        return -1
    else:
        ma=-1
        ar=-1
        for j in range(n):
            if G[i][j]==1:
                if fp(G,j)>ma:
                    ma=fp(G,j)
                    ar=j
        return ar

def fm(G,i): # calculate f-(i) for gene i, 
# where G is the auxiliary graph
    if i==0:
        return 0
    else:
        return 1+fm(G,hm(G,i))

def hm(G,i): # calculate h-(i) for gene i, 
# where G is the auxiliary graph
    n=len(G)
    if i==0:
        return -1
    else:
        ma=-1
        ar=-1
        for j in range(n):
            if G[j][i]==1:
                if fm(G,j)>ma:
                    ma=fm(G,j)
                    ar=j
        return ar

#########################

A=[[0,1,2,3,4,5,6,7],
   [0,1,4,2,3,5,6,7],
   [0,1,2,4,3,6,5,7]]
# A is the input of m sequences, where each sequence
# has n genes 0,...,n-1 (no replicated genes)
# here 0 (always at the head) and n-1 (always at the tail) 
# are auxiliary

m=len(A) # number of sequences
n=len(A[0]) # number of genes

G=[([0]*n) for i in range(n)]
# the auxiliary graph G. G[i][j]=1 means there is an edge from i to j
for k in range(m):
    for i in range(n):
        for j in range(i+1,n):
            c=A[k][i]
            d=A[k][j]
            G[c][d]=G[c][d]+1 
for i in range(n):
    for j in range(n):
        if G[i][j]==m:
            G[i][j]=1
        else:
            G[i][j]=0
    
f1=([0]*n) # f+(i) for all genes
h1=([0]*n) # h+(i) for all genes
f2=([0]*n) # f-(i) for all genes
h2=([0]*n) # h-(i) for all genes
for i in range(n):
    f1[i]=fp(G,i)
    h1[i]=hp(G,i)    
    f2[i]=fm(G,i) 
    h2[i]=hm(G,i)

l0=([0]*n) # l0[i]=1 means gene i is in the 
# longest common subsequence L_0
l0[0]=1;
cl=0; # most recent gene found in generating L_0
while h1[cl]!=-1:
    cl=h1[cl]
    l0[cl]=1
    
tr=([1]*n) # whether each gene is a proper/quasi/non-transposon
# 1 is non-transposon; 0 is quasi-transposon; -1 is proper-transposon
#C0=sum(l0)

for i in range(1,n-1):
    if l0[i]==1: # only consider genes not in L_0
        continue
    if f1[i]+f2[i]<f1[0]: # i is a proper-transposon
        tr[i]=-1
    else:
        tr[i]=0 # i is a quasi-transposon
        lt=([0]*n) # the longest common subsequence L_i
        lt[i]=1
        cl=i
        while h1[cl]!=-1: # find L_i forward
            cl=h1[cl]
            lt[cl]=1
        cl=i
        while h2[cl]!=-1: # find L_i backward
            cl=h2[cl]
            lt[cl]=1
        for j in range(n):
            if l0[j]==1 and lt[j]==0:
                tr[j]=0
        # a gene in L_0 but not L_i is also a quasi-transposon

# a gene that is not determined to be proper/quasi-transposon here
# is a non-transposon, so that tr[i]=1

for i in range(1,n-1): # output
    if tr[i]==1:
        print(i,'is a non-transposon')
    else:
        if tr[i]==0:
            print(i,'is a quasi-transposon')
        else:
            print(i,'is a proper-transposon')




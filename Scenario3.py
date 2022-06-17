#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Given some linear gene sequences, where each gene can have multiple 
copies, this code determines whether each gene is a transposon. 
This algorithm does not always produce the correct result.
See below for a counterexample.
This code implements Algorithm 5 for Scenario 3 in my paper
https://arxiv.org/abs/1506.02424

@author: yuewang
"""

def match(A,n): # determine if each gene pair i,j have their copies keep
# relative locations in all sequences in A
# A is m sequences of n genes, and the length of each sequence is l
    m=len(A)
    l=len(A[0])
    B=[([1]*n) for i in range(n)] # B[i][j]=1 if copies of genes i,j
    # keep their relative locations in all sequences
    for i in range(n):
        for j in range(n):
            ar=[] # extract copies of genes i,j in A[0]
            for k in range(l):
                if A[0][k]==i or A[0][k]==j:
                    ar.append(A[0][k])
            for c in range(1,m):
                ar2=[] # extract copies of genes i,j in A[c] and compare
                # with that of A[0]
                for k in range(l):
                    if A[c][k]==i or A[c][k]==j:
                        ar2.append(A[c][k])
                if ar!=ar2:
                    B[i][j]=0
                    break
    return B

#########################

A=[[0,1,2,3,3,2,1],
   [1,2,3,3,2,1,0]]

"""
# the following example fails this algorithm
# the correct outcome should be 7,8,9,0
A=[[7,8,9,0,1,1,2,3,3,4,5,5,6],
   [1,2,1,3,4,3,5,6,5,7,8,9,0]]
"""

# input of m sequences of n genes, numbered 0,...,n-1.
n=max(A[0])+1
B=match(A,n)
ind=([0]*n) # genes that have not been deleted
for i in range(n):
    ind[i]=i
while len(ind)>0:
    nt=len(ind) # number of genes not deleted
    d=([0]*nt)
    for i in range(nt):
        d[i]=sum(B[i]) # calculate the degree of each gene
    minv=min(d)
    if minv==nt: # the minimal degree+1=number of genes not deleted,
    # meaning that the genes not deleted form a clique
        break
    else:
        indv=d.index(minv) # delete the gene with the minimal degree
        # and delete corresponding elements in B
        del ind[indv]
        for j in range(nt):
            del B[j][indv]
        del B[indv]
    
tr=([0]*n) # tr[i]=1 if gene i is not deleted,
# i.e., is in the "maximum" clique
for i in range(len(ind)):
    tr[ind[i]]=1
for i in range(n):
    if tr[i]==1:
        print(i,'is not a transposon')
    else:
        print(i,'is a transposon')




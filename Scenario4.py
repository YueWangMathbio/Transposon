#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Given some circular gene sequences, where each gene can have multiple 
copies, this code determines whether each gene is a transposon by
finding the largest complete subgraph.
B is actually a 3-uniform hypergraph. 
This algorithm does not always produce the correct result.
See below for a counterexample.
This code implements Algorithm 6 for Scenario 4 in my paper
https://arxiv.org/abs/1506.02424

@author: yuewang
"""

import numpy as np

def cycle(a): # rotate a circular sequence
    l=len(a)
    b=([0]*l)
    for i in range(1,l):
        b[i]=a[i-1]
    b[0]=a[l-1]
    return b

def compare(a1,a2): # compare two circular sequences by rotating them
    l1=len(a1)
    l2=len(a2)
    cr=0 # 1 means a1 and a2 are equal
    if l1==l2:
        for i in range(l1):
            if a1==a2:
                cr=1
                break
            else:
                a1=cycle(a1)
    return cr   

def match(A,n): # determine if each gene triple i,j,k have their 
# copies keep relative locations in all sequences in A
# A is m sequences of n genes, and the length of each sequence is l
    m=len(A)
    l=len(A[0])
    B=np.ones((n,n,n)) # B[i,j,k]=1 if copies of genes i,j,k
    # keep their relative locations in all sequences
    for i in range(n):
        for j in range(n):
            for k in range(n):
                ar=[] # extract copies of genes i,j,k in A[0]
                for p in range(l):
                    if A[0][p]==i or A[0][p]==j or A[0][p]==k:
                        ar.append(A[0][p])
                for c in range(1,m):
                    ar2=[] # extract copies of genes i,j,k in A[c] 
                    # and compare with that of A[0]
                    l=len(A[c])
                    for q in range(l):
                        if A[c][q]==i or A[c][q]==j or A[c][q]==k:
                            ar2.append(A[c][q])
                    if compare(ar,ar2)==0:
                        B[i,j,k]=0
                        break
    return B

#########################

A=[[0,1,2,3,4,1],
   [1,1,2,0,3,4]]

"""
# the following example fails this algorithm
# the correct outcome should be 7,8,9,0
A=[[1,2,7,3,4,8,9,5,6,0],
   [2,1,0,4,3,7,8,6,5,9],
   [1,2,9,3,4,0,7,5,6,8],
   [2,1,8,4,3,9,0,6,5,7]]
"""

n=max(A[0])+1
B=match(A,n) # the 3-uniform hypergraph
ind=([0]*n) # genes that have not been deleted
for i in range(n):
    ind[i]=i
while len(ind)>0:
    nt=len(ind) # number of genes not deleted
    d=([0]*nt)
    for i in range(nt):
        d[i]=sum(sum(B[i,:,:])) # calculate the degree of each gene
    minv=min(d)
    if minv==nt**2: # the genes not deleted form a clique
        break
    else:
        indv=d.index(minv) # delete the gene with the minimal degree
        # and delete corresponding elements in B
        #print(d)
        del ind[indv]
        B=np.delete(B,indv,0)
        B=np.delete(B,indv,1)
        B=np.delete(B,indv,2)
        #print(B)
    
tr=([0]*n) # tr[i]=1 if gene i is not deleted,
# i.e., is in the "maximum" clique
for i in range(len(ind)):
    tr[ind[i]]=1
for i in range(n):
    if tr[i]==1:
        print(i,'is not a transposon')
    else:
        print(i,'is a transposon')




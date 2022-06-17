#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Given some circular gene sequences, where each gene has only one copy 
(meaning that each sequence is a permutation of 1 to n-2, but we can 
rotate them),this code determines whether each gene is a 
proper/quasi/non-transposon. 
This code implements Algorithms 3,4 for Scenario 2 in my paper
https://arxiv.org/abs/1506.02424

@author: yuewang
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
 
def ag(A): # construct the auxiliary graph G for Scenario 1
    m=len(A)
    n=len(A[0])
    G=[([0]*n) for i in range(n)]
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
    return G

def alg1(G): # apply Algorithm 1 to determine one longest common
# subsequence from 0 to n-1, given the auxiliary graph G
    n=len(G)
    f1=([0]*n) # f+(i) for all genes
    h1=([0]*n) # h+(i) for all genes
    for i in range(n):
        f1[i]=fp(G,i)
    for i in range(n):    
        h1[i]=hp(G,i)   
    l0=([0]*n)
    l0[0]=1;
    cl=0;
    while h1[cl]!=-1:
        cl=h1[cl]
        l0[cl]=1
    return l0

def alg2(G): # apply Algorithms 1,2 to determine whether each gene is
# proper/quasi/non-transposon
# in the returned tr, 1 is non-transposon; 0 is quasi-transposon; 
# -1 is proper-transposon
    n=len(G)
    f1=([0]*n)
    h1=([0]*n)
    f2=([0]*n)
    h2=([0]*n)
    for i in range(n):
        f1[i]=fp(G,i)
        h1[i]=hp(G,i)
        f2[i]=fm(G,i)
        h2[i]=hm(G,i)
    l0=([0]*n) 
    l0[0]=1;
    cl=0; 
    while h1[cl]!=-1:
        cl=h1[cl]
        l0[cl]=1        
    tr=([1]*n) 
    for i in range(1,n-1):
        if l0[i]==1: 
            continue
        if f1[i]+f2[i]<f1[0]: 
            tr[i]=-1
        else:
            tr[i]=0 
            lt=([0]*n) 
            lt[i]=1
            cl=i
            while h1[cl]!=-1:
                cl=h1[cl]
                lt[cl]=1
            cl=i
            while h2[cl]!=-1:
                cl=h2[cl]
                lt[cl]=1
            for j in range(n):
                if l0[j]==1 and lt[j]==0:
                    tr[j]=0
    return tr

def cut(i,A): # cut the sequences at i and add auxiliary 0 and n-1
    m=len(A)
    n=len(A[0])+2
    B=[([-1]*n) for j in range(m)]
    for j in range(m):
        B[j][0]=0
        B[j][n-1]=n-1
        loc=A[j].index(i) # where i appears in sequence j
        for k in range(loc,n-2):
            B[j][k-loc+1]=A[j][k]
        for k in range(0,loc):
            B[j][n-loc+k-1]=A[j][k]
    return B

#########################

A=[[6,7,1,2,3,4,5],
   [1,4,2,3,5,6,7],
   [6,5,7,1,2,4,3]]
# A is the input of m circular sequences, where each sequence
# has n-2 genes 1 to n-2 (no replicated genes)
# sequences can rotate

m=len(A) # number of sequences
n=len(A[0])+2 # length of sequence (adding auxiliary 0 and n-1)

lc=([0]*n) # lc[i]=1 if gene i has been chosen and cut

cc=([-1]*n) # cc[i] is the length of the longest common subsequence
# when the circular sequences are cut at i

B=cut(1,A) # cut at gene 1
G=ag(B) # generate the auxiliary graph
L=alg1(G) # find a longest common subseequence in the linear case
# this is the complement of S
C=sum(L)
lc[1]=1
cc[1]=C

while sum(lc)+sum(L)<n+1: # total number of cut is k+1
    for i in range(n):
        if lc[i]==0 and L[i]==0: 
        # find a gene in S which has not been cut
            break
    B=cut(i,A)
    G=ag(B)
    Lt=alg1(G)
    Ct=sum(Lt)
    lc[i]=1
    cc[i]=Ct
    if Ct>C: # update C and L (the complement of S)
        C=Ct
        L=Lt

tr=([1]*n) # for proper-transposon, tr[i]=-1; for quasi-transposon, 
# tr[i]=0; for non-transposon, tr[i]=-1
for i in range(n):
    if L[i]==0:
        if cc[i]<C:
            tr[i]=-1 # determine proper-transposons
        else:
            tr[i]=0
            B=cut(i,A)
            G=ag(B)
            trt=alg2(G)
            for j in range(n): # determine quasi-transposons
                if L[j]==1 and trt[j]<=0:
                    tr[j]=0

for i in range(1,n-1): # output
    if tr[i]==1:
        print(i,'is a non-transposon')
    else:
        if tr[i]==0:
            print(i,'is a quasi-transposon')
        else:
            print(i,'is a proper-transposon')





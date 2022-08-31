# -*- coding: utf-8 -*-
"""
Given a random undirected graph with 15 vertices, the first half calculates
the size of the maximum clique. The second half uses a heuristic algorithm
to approximately calculate the size of the maximum clique. Then compare them
to see if the heuristic algorithm produces the correct result.

This code tests the performance of Algorithm 5 for Scenario 3 in my paper
https://arxiv.org/abs/1506.02424

@author: yuewang
"""
import random


n=15 # number of vertices in this random graph
tt=10000 # repeat for 10000 times
c=0 # number of runs that two methods match
for x in range(tt):
    #if x%100==0:
    #    print(x)
    B=[([1]*n) for i in range(n)]
    for i in range(n):
        for j in range(i+1,n):
            B[i][j]=random.randint(0, 1)
            B[j][i]=B[i][j] # generate the connect matrix
    
    lc=0          
    for i in range(2**n):
        ind=([0]*n)
        temp=i
        for j in range(n): # generate each subset
            ind[j]=temp%2
            temp=temp//2
        te=1
        for j in range(n):
            for k in range(n):
                if ind[j]==1 and ind[k]==1 and B[j][k]==0: # check if a clique
                    te=0
                    break
            if te==0:
                break
        if te==0:
            continue
        else:
            tc=sum(ind)
            if tc>lc:
                lc=tc
    talg=lc # true size of the maximum clique


    ind=([0]*n)
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
    alg5=len(ind) # result of Alg. 5
    
    if talg==lc:
        c=c+1
    else:
        print(B)

print(c,tt)
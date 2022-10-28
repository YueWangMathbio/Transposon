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

def tsu(n,G,c,v,s): # auxiliary for topological sorting
    v[c]=True
    for i in range(n):
        if G[c][i]==1 and v[i]==False:
            tsu(n,G,i,v,s)
    s.append(c)
 
def ts(n,G): # topological sorting for genes according to the DAG
    v=[False]*n
    s=[]
    for i in range(n):
        if v[i]==False:
            tsu(n,G,i,v,s)
    r=[]
    while s:
        z=s.pop()
        r.append(z)
    return r

def fh(n,G): # calculate f+, h+, f-, h- 
    r=ts(n,G)
    f1=[-1]*n
    h1=[-1]*n
    f1[n-1]=0
    for i in range(n-2,-1,-1):
        x=r[i]
        for j in range(i+1,n):
            y=r[j]
            if G[x][y]==1:
                if f1[y]+1>f1[x]:
                    f1[x]=f1[y]+1
                    h1[x]=y
    f2=[-1]*n
    h2=[-1]*n
    f2[0]=0
    for i in range(1,n):
        x=r[i]
        for j in range(i):
            y=r[j]
            if G[y][x]==1:
                if f2[y]+1>f2[x]:
                    f2[x]=f2[y]+1
                    h2[x]=y
    return f1,h1,f2,h2
 

#########################

dataset=1 
# dataset=1 means all three variants of ST540;
# 2 means last two variants of ST540
# 3 means all three variants of ST2747

if dataset<=2:
    with open('CP007265.1.txt') as f:
        lines = f.readlines()
    r1=lines
    with open('CP007390.1.txt') as f:
        lines = f.readlines()
    r2=lines
    with open('CP007391.1.txt') as f:
        lines = f.readlines()
    r3=lines
else:
    with open('CP007392.1.txt') as f:
        lines = f.readlines()
    r1=lines
    with open('CP007393.1.txt') as f:
        lines = f.readlines()
    r2=lines
    with open('CP007394.1.txt') as f:
        lines = f.readlines()
    r3=lines
  
dic={}
p=1
s1=[]
for w in r1:
    dic[w]=p
    s1.append(p)
    p=p+1    
s2=[]
for w in r2:
    s2.append(dic[w])
s3=[]
for w in r3:
    s3.append(dic[w])
s1=[0]+s1+[len(r1)]
s2=[0]+s2+[len(r1)]
s3=[0]+s3+[len(r1)]
invdic={}
for w in dic:
    invdic[dic[w]]=w

if dataset!=2:
    A=[s1,s2,s3]
else:
    A=[s2,s3]



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
    

f1,h1,f2,h2=fh(n,G) # f+(i), h+(i), f-(i), h-(i) for all genes

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
nt=[]
qt=[]
pt=[]
for i in range(1,n-1): # output
    if tr[i]==1:
        nt.append(invdic[i])
    elif tr[i]==0:
        qt.append(invdic[i])
    else:
        pt.append(invdic[i])    
            
print('non-transposons: ',nt)        
print('quasi-transposons: ',qt)        
print('proper-transposons: ',pt)  

print('number of non-transposons: ',len(nt))        
print('number of quasi-transposons: ',len(qt))        
print('number of proper-transposons: ',len(pt))  
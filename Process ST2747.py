#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Process three sequencing files of Escherichia coli strain ST2747.
First take out all genes. Then delete genes that appear more than once
in one sequence, and genes that do not appear in all three sequences.

@author: yuewang
"""

# select genes
with open('Escherichia coli strain ST2747 GenBank CP007392.1.txt') as f:
    lines = f.readlines()
l=len(lines)
r1=[]
for i in range(l):
    tn=len(lines[i])
    for j in range(tn-6):
        if lines[i][j:j+6]=='[gene=':
            z=''
            p=j+6
            while lines[i][p]!=']':
                z=z+lines[i][p]
                p=p+1
            r1.append(z)

with open('Escherichia coli strain ST2747 GenBank CP007393.1.txt') as f:
    lines = f.readlines()
l=len(lines)
r2=[]
for i in range(l):
    tn=len(lines[i])
    for j in range(tn-6):
        if lines[i][j:j+6]=='[gene=':
            z=''
            p=j+6
            while lines[i][p]!=']':
                z=z+lines[i][p]
                p=p+1
            r2.append(z)

with open('Escherichia coli strain ST2747 GenBank CP007394.1.txt') as f:
    lines = f.readlines()
l=len(lines)
r3=[]
for i in range(l):
    tn=len(lines[i])
    for j in range(tn-6):
        if lines[i][j:j+6]=='[gene=':
            z=''
            p=j+6
            while lines[i][p]!=']':
                z=z+lines[i][p]
                p=p+1
            r3.append(z)

# delete genes that appear more than once
dic={}
for w in r1:
    if w not in dic:
        dic[w]=1
    else:
        dic[w]=2
        
for i in range(len(r1)-1,-1,-1):
    if dic[r1[i]]==2:
        r1.pop(i)

dic={}
for w in r2:
    if w not in dic:
        dic[w]=1
    else:
        dic[w]=2
for i in range(len(r2)-1,-1,-1):
    if dic[r2[i]]==2:
        r2.pop(i)
                    
dic={}
for w in r3:
    if w not in dic:
        dic[w]=1
    else:
        dic[w]=2
for i in range(len(r3)-1,-1,-1):
    if dic[r3[i]]==2:
        r3.pop(i)

# delete genes that do not appear in all sequences
for i in range(len(r1)-1,-1,-1):
    if r1[i] not in r2 or r1[i] not in r3:
        r1.pop(i)
        
for i in range(len(r2)-1,-1,-1):
    if r2[i] not in r1 or r2[i] not in r3:
        r2.pop(i)
        
for i in range(len(r3)-1,-1,-1):
    if r3[i] not in r1 or r3[i] not in r2:
        r3.pop(i)

# output
with open(r'CP007392.1.txt', 'w') as fp:
    for item in r1:
        fp.write("%s\n" % item)
        
with open(r'CP007393.1.txt', 'w') as fp:
    for item in r2:
        fp.write("%s\n" % item)
        
with open(r'CP007394.1.txt', 'w') as fp:
    for item in r3:
        fp.write("%s\n" % item)
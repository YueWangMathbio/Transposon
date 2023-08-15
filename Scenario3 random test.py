# -*- coding: utf-8 -*-
"""
Given a random undirected graph with n vertices, the first half calculates
the size of the maximum clique. Here a clique is a subset of vertices that 
each pair has an edge that links them.
The second half implements Algorithm 5 for Scenario 3 in my paper
https://arxiv.org/abs/2301.03827 which is a heuristic algorithm
to approximately calculate the size of the maximum clique. 
Then we compare them to see if the heuristic algorithm produces the correct 
result.
Since determining the size of the true maximum clique is very slow, we only 
test for small graphs.
@author: yuewang
"""
import random

edge_probability = 0.7 # probability that one edge exists in the random graph
success_rate = [] # record the successs rate (Algorithm 5 produces the correct
# result) for n=4 to n=12 (larger n needs much more time)
for n in range(4, 13):
    # n is the number of vertices in this random graph
    total_simulation_count = 10000 # repeat for 10000 times
    total_success_count = 0 # number of runs that two methods match
    for x in range(total_simulation_count):
        # run the simulation total_simulation_count times
        if x % 1000 == 0:
            print(n, x)
        graph_edge = [[1] * n for i in range(n)] # whether there is an edge between
        # two vertices. graph_edge[i][j] == 1 means there is an edge between i and
        # j, and graph_edge[i][j] == 0 otherwise.
        for i in range(n):
            for j in range(i+1, n): # for each pair i, j
                graph_edge[i][j] = random.choices([0, 1], \
                        weights=(1-edge_probability, edge_probability))[0]
                graph_edge[j][i] = graph_edge[i][j] # generate a random edge
        
        maximum_clique_size = 0 # the size of the largest clique found so far  
    # for each subset of 1,...,n, check whether every two vertices in this subset
    # has an edge linking them, so that the subset is a clique.
        for i in range(2**n): # for each i from 0 to 2 ^ n - 1
            binary_expression = [0] * n # generate the binary expression of i,
            # which represents a subset of 1,...,n.
            temp = i
            for j in range(n): # generate each subset
                binary_expression[j] = temp % 2
                temp = temp // 2 # use the mod 2 method to generate the binary
                # expression of i
            termination_marker = False # whether we have found that the current
            # subset is not a clique
            for j in range(n):
                for k in range(n):
                    if binary_expression[j] == 1 and binary_expression[k] == 1 \
                        and graph_edge[j][k] == 0: # check if it is not a clique
                        termination_marker = True
                        break
                if termination_marker == True:
                    break
            if termination_marker == True:
                continue
            else:
                vertex_count = sum(binary_expression)
                if vertex_count > maximum_clique_size: # check whether the 
                # new-found clique is the largest
                    maximum_clique_size = vertex_count
        true_result = maximum_clique_size # true size of the maximum clique
    
        
    # the following is Algorithm 5, which delete the vertex with the smallest 
    # degree until we have a clique
        degree = [sum(graph_edge[i]) - 1 for i in range(n)]
        # degree[i] is the degree of vertex i in this graph, namely the number of 
        # edges linking i
        vertex_number = n # number of vertices that have not been deleted
        deleted = [False] * n # whether gene i has been deleted
        while vertex_number > 0: # when there is at least one gene left
            min_degree = float('inf') # the minimum degree have seen
            min_degree_index = -1 # the index of the minimum degree
            for i in range(n):
                if deleted[i] == False:
                    if degree[i] < min_degree:
                        min_degree = degree[i]
                        min_degree_index = i
            if min_degree == vertex_number - 1: # we find a clique
                break
            else:
                # min_degree_index is the gene that will be deleted    
                for vertex in range(n):
                    if deleted[vertex] == False \
                        and graph_edge[min_degree_index][vertex] == 1:
                        degree[vertex] -= 1 # update degree
                deleted[min_degree_index] = True
                vertex_number -= 1
        alg5_result = vertex_number
        
        
    
        
        if true_result == alg5_result:
            total_success_count = total_success_count + 1 
            # the result of Alg. 5 matches the true result.
            # add the success count by 1.
        else:
            pass
            #print(graph_edge) # show the counterexample
            #print(true_result, alg5_result)
    success_rate.append(float(total_success_count) / float(total_simulation_count))
print('success rate for n=4 to n=12: ', success_rate)
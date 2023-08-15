# -*- coding: utf-8 -*-
"""
Given a random undirected 3-uniform hypergraph with n vertices, the first half 
calculates the size of the maximum clique. Here a clique is a subset of 
vertices that each triple has a hyperedge that links them.
The second half implements Algorithm 6 for Scenario 4 in my paper
https://arxiv.org/abs/2301.03827 which is a heuristic algorithm
to approximately calculate the size of the maximum clique. 
Then we compare them to see if the heuristic algorithm produces the correct 
result.
Since determining the size of the true maximum clique is very slow, we only 
test for small graphs.
@author: yuewang
"""
import random


edge_probability = 0.01 # probability that one hyperedge exists in the random 
# graph
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
        graph_edge = [[[1] * n for i in range(n)] for j in range(n)] 
        # whether there is a hyperedge among
        # three vertices. graph_edge[i][j][k] == 1 means there is a hyperedge 
        # among i, j, k, and graph_edge[i][j][k] == 0 otherwise.
        for i in range(n):
            for j in range(i+1, n): # for each triple i, j, k
                for k in range(j+1, n):
                    graph_edge[i][j][k] = random.choices([0, 1], \
                            weights=(1-edge_probability, edge_probability))[0]
        
        maximum_clique_size = 0 # the size of the largest clique found so far  
    # for each subset of 1,...,n, check whether every two vertices in this subset
    # has an edge linking them, so that the subset is a clique.
        for number in range(2**n): # for each number from 0 to 2 ^ n - 1
            binary_expression = [0] * n # generate the binary expression of i,
            # which represents a subset of 1,...,n.
            temp = number
            for j in range(n): # generate each subset
                binary_expression[j] = temp % 2
                temp = temp // 2 # use the mod 2 method to generate the binary
                # expression of i
            termination_marker = False # whether we have found that the current
            # subset is not a clique
            for i in range(n):
                for j in range(i+1, n):
                    for k in range(j+1, n):
                        if binary_expression[i] == 1 \
                            and binary_expression[j] == 1 \
                                and binary_expression[k] == 1 \
                            and graph_edge[i][j][k] == 0: # check if it is not a clique
                            termination_marker = True
                            break
                    if termination_marker == True:
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
    
        
    # the following is Algorithm 6, which delete the vertex with the smallest 
    # degree until we have a clique
    
        degree = [0] * n
        # degree[i] is the degree of vertex i in this graph, namely the number of 
        # 3-hyperedges linking i
        for i in range(n):
            for j in range(i+1, n):
                for k in range(j+1, n): # for each gene triple i, j, k
                    if graph_edge[i][j][k] == True: # if they form a 3-hyperedge
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
                        if graph_edge[i][j][min_degree_index] == True \
                            and deleted[i] == False and deleted[j]==False:
                            # if there is a 3-hyperedge that has not been deleted
                            degree[i] -= 1
                            degree[j] -= 1
                            degree[min_degree_index] -= 1
                # remove the 3-hyperedge of i, j, min_degree_index
                for i in range(min_degree_index): # for i < j < min_degree_index
                    for j in range(min_degree_index+1, n):
                        if graph_edge[i][min_degree_index][j] == True \
                            and deleted[i] == False and deleted[j]==False:
                            # if there is a 3-hyperedge that has not been deleted
                            degree[i] -= 1
                            degree[j] -= 1
                            degree[min_degree_index] -= 1
                # remove the 3-hyperedge of i, j, min_degree_index
                for i in range(min_degree_index+1, n): # for i < j < min_degree_index
                    for j in range(i+1, n):
                        if graph_edge[min_degree_index][i][j] == True \
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
        alg6_result = gene_number
        
    
        
        if true_result == alg6_result:
            total_success_count = total_success_count + 1 
            # the result of Alg. 5 matches the true result.
            # add the success count by 1.
        else:
            pass
            #print(graph_edge) # show the counterexample
            #print(true_result, alg5_result)
    success_rate.append(float(total_success_count) / float(total_simulation_count))
print('edge generating probability: ', edge_probability)
print('success rate for n=4 to n=12: ', success_rate)


"""
resutls

edge generating probability:  0.01
success rate for n=4 to n=12:  [1.0, 1.0, 1.0, 1.0, 1.0, 
                                1.0, 1.0, 1.0, 1.0]

edge generating probability:  0.1
success rate for n=4 to n=12:  [1.0, 1.0, 1.0, 1.0, 0.9996, 
                                0.9993, 0.9969, 0.9947, 0.9893]

edge generating probability:  0.3
success rate for n=4 to n=12:  [1.0, 1.0, 0.9971, 0.972, 0.9217, 
                                0.8381, 0.7637, 0.7138, 0.7111]

edge generating probability:  0.5
success rate for n=4 to n=12:  [1.0, 1.0, 0.9676, 0.9167, 0.9089, 
                                0.921, 0.919, 0.8737, 0.7896]

edge generating probability:  0.7
success rate for n=4 to n=12:  [1.0, 1.0, 0.9741, 0.9652, 0.9142, 
                                0.8537, 0.8439, 0.8406, 0.7919]

edge generating probability:  0.9
success rate for n=4 to n=12:  [1.0, 1.0, 0.9998, 0.9854, 0.9559, 
                                0.9255, 0.898, 0.8402, 0.8177]

edge generating probability:  0.99
success rate for n=4 to n=12:  [1.0, 1.0, 1.0, 1.0, 0.9999, 
                                0.9997, 0.9973, 0.993, 0.9864]

"""







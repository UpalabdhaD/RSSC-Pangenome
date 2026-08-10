import sys
import matplotlib.pyplot as plt
import networkx as nx
# import pandas as pd
# from functools import reduce

pw_ani = sys.argv[1]
threshold = int(sys.argv[2])
cluster_membership = sys.argv[3]
filtered_network = sys.argv[4]


g = nx.Graph()

cnt = 0

new_network = []

# read the df, parse and add to network in weight is > threshold
with open(pw_ani) as f:
    next(f)
    for line in f:
        # print(line)
        toks = line.strip().split(',')
        node1 = toks[0]
        node2 = toks[1]
        weight = float(toks[2])
        # print(node1, node2, weight)
        
        # add the node
        g.add_node(node1)
        g.add_node(node2)

        if weight >= threshold:
            g.add_edge(node1, node2, weight = weight)
            new_network.append(','.join(line.split(' ')[:3]))
        
        cnt += 1


# write the pairwise filtered network
with open(filtered_network, 'w') as nw:
    nw.write("node1,node2,weight\n")
    for line in new_network:
        nw.write(line)



cc = nx.connected_components(g)

# write cluster for each genome/node
cluster_num = 0
with open(cluster_membership, 'w') as cf:
    for c in cc:
        c = list(c)
        cluster_num += 1
        for node in c:
            cf.write(node + ',' + str(cluster_num) + '\n')
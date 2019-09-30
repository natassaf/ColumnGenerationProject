import pandas as pd
import pickle
from itertools import chain
import time
from sage.all import *
import sys
import numpy as np
import local_search_new as ls

def readcsv(name_of_file):
    with open(name_of_file + '.csv') as csv_file:
        data_pd =  pd.read_csv(csv_file, skiprows=1, 
        names=['day_string','month', 'day_num', "time", "time_zone",  "year", "seconds" ,"departure_node", "arrival_node"], sep=' |;', engine='python')
    return data_pd

def create_edges(vertices, P, D):
    edges = []
    for v in vertices:
        for u in vertices:
            if v != u:
                if v == "s":
                    if u in D:
                        continue
                    edges.append((v,u))
                elif v in P:
                    if u == "t" or u == 's':
                        continue
                    edges.append((v,u))
                elif v in D:
                    if u == 's':
                        continue
                    edges.append((v,u))
                elif v != "t":
                    edges.append((v,u)) 
    return edges

def delete_edges(g, P, D, n):
    g.delete_edge(2*n+1,0)
    for d in D:
        g.delete_edge(0 , d)
        g.delete_edge(2*n+1,d)
        g.delete_edge(d,0)

    for p in P:
        g.delete_edge(n+p,p) #new
        g.delete_edge(p, 2*n+1)
        g.delete_edge(2*n+1, p)
        g.delete_edge(p,0)
    return g
    
def travel_time(g, from_node, to_node):
        """
        Gets the travel times between two locations.
        In these data we get as costs in the file the travel time. If we had the distance in meters we would use the function in the comment
        """
        #travel_time = (((g.edge_label(from_node,to_node)/1000)*3600)/80)
        travel_time = g.edge_label(from_node,to_node)
        return travel_time 

def eliminate_arc_1(g, n, P, D, data_dict):
    #deleting edges of the form (j, n+i) if the path (j,i,n+j,n+i) is infeasible
    unable_edges = set()
    for iter1 in P:
        i= iter1
        for iter2 in P:
            j = iter2
            if i != j:
                U = []
                
                path_checking = [i,j,n+i,n+j]

                path = []
                L = ls.Label(i, data_dict["time_windows"][i][0], data_dict["load"][i], 0, [], [], P, D, -1, n, 0, 0, {})
                path.append(L)
                U.append(path)
                
                    
                while U:
                    path = U.pop() #U is a queue first in first out!
                    L = path[-1]
                    current_node = path[-1].id
                    if current_node == path_checking[0]:
                        max_travel_time_of_user = data_dict["max_ride_time_per_user"][i]
                        next_node = path_checking[1]
                    elif current_node == path_checking[1]:
                        max_travel_time_of_user = data_dict["max_ride_time_per_user"][j]
                        next_node = path_checking[2]
                    elif current_node == path_checking[2]:
                        max_travel_time_of_user = min(data_dict["max_ride_time_per_user"][i], data_dict["max_ride_time_per_user"][j])
                        next_node = path_checking[3]
                    else:
                        continue
                    if (L.t + travel_time(g, current_node, next_node) > data_dict["time_windows"][next_node][1] or L.l + data_dict["load"][current_node] > data_dict["vehicle_capacity"] 
                    or L.t + travel_time(g, current_node, next_node) > max_travel_time_of_user):
                    
                        g.delete_edge(j, n+i)
                        unable_edges.add((j,n+i))
                    else:
                        time_node = L.t + travel_time(g, current_node, next_node)
                        cummulative_load = L.l + data_dict["load"][next_node]
                        L_new = ls.Label(next_node, time_node, cummulative_load, 0, [], [], P, D, current_node, n, 0, 0, {})
                        path.append(L_new)
                        U.append(path)
    return g,unable_edges 
    
def eliminate_arc_2(g, n, P, D, data_dict, unable_edges):
    #deleting edges of the form (n+i, j) if the path (i,n+i, j ,n+j) is infeasible
    for iter1 in P:
        i= iter1
        for iter2 in P:
            j = iter2
            if i != j:
                U = []
                if (i,j) in unable_edges:
                    continue
                path_checking = [i,n+i, j ,n+j]
                path = []
                L = ls.Label(i, data_dict["time_windows"][i][0], data_dict["load"][i], 0, [], [], P, D, -1, n, 0, 0, {})
                path.append(L)
                U.append(path)
                  
                while U:
                    path = U.pop() #U is a queue first in first out!
                    L = path[-1]
                    current_node = path[-1].id
                    if current_node == path_checking[0]:
                        max_travel_time_of_user = data_dict["max_ride_time_per_user"][i]
                        next_node = path_checking[1]
                    elif current_node == path_checking[1]:
                        max_travel_time_of_user = data_dict["max_ride_time_per_user"][j]
                        next_node = path_checking[2]
                    elif current_node == path_checking[2]:
                        max_travel_time_of_user = min(data_dict["max_ride_time_per_user"][i], data_dict["max_ride_time_per_user"][j])
                        next_node = path_checking[3]
                    else:
                        continue
                    if (L.t + travel_time(g, current_node, next_node) > data_dict["time_windows"][next_node][1] or 
                    L.l + data_dict["load"][current_node] > data_dict["vehicle_capacity"]):
                        g.delete_edge(n+i, j)
                        unable_edges.add((n+i, j))
                    else:
                        time_node = L.t + travel_time(g, current_node, next_node)
                        cummulative_load = L.l + data_dict["load"][next_node]
                        L_new = ls.Label(next_node, time_node, cummulative_load, 0, [], [], P, D, current_node, n, 0, 0, {})
                        path.append(L_new)
                        U.append(path)
    return g,unable_edges       

def eliminate_arc_3(g, n, P, D, data_dict, unable_edges):
    #deleting edges of the form (i, j) if the paths (i,n+i, j ,n+j) and the (i, j, n+j, n+i) are infeasible 
    
    for iter1 in P:
        i= iter1
        for iter2 in P:
            j = iter2
            if i != j:
                flag = 0
                U = []
                if (i,j) in unable_edges:
                    continue
                elif (j,n+i) in unable_edges or (n+i,n+j) in unable_edges:
                    unable_edges.add((i, j))
                    g.delete_edge(i,j)
                    
                path_checking = [i, j, n+i ,n+j]
                if not g.has_edge(i,j):
                    continue
                

                path = []
                L = ls.Label(i, data_dict["time_windows"][i][0], data_dict["load"][i], 0, [], [], P, D, -1, n, 0, 0, {})
                path.append(L)
                U.append(path)
                  
                while U:
                    path = U.pop() #U is a queue first in first out!
                    L = path[-1]
                    current_node = path[-1].id
                    if current_node == path_checking[0]:
                        max_travel_time_of_user = data_dict["max_ride_time_per_user"][i]
                        next_node = path_checking[1]
                    elif current_node == path_checking[1]:
                        max_travel_time_of_user = data_dict["max_ride_time_per_user"][j]
                        next_node = path_checking[2]
                    elif current_node == path_checking[2]:
                        max_travel_time_of_user = min(data_dict["max_ride_time_per_user"][i], data_dict["max_ride_time_per_user"][j])
                        next_node = path_checking[3]
                    else:
                        continue
                    if (L.t + travel_time(g, current_node, next_node) > data_dict["time_windows"][next_node][1] or
                    L.t + travel_time(g, current_node, next_node) > max_travel_time_of_user or 
                    L.l + data_dict["load"][current_node] > data_dict["vehicle_capacity"]):
                        flag +=1
                    else:
                        time_node = L.t + travel_time(g, current_node, next_node)
                        cummulative_load = L.l + data_dict["load"][next_node]
                        L_new = ls.Label(next_node, time_node, cummulative_load, 0, [], [], P, D, current_node, n, 0, 0, {})
                        path.append(L_new)
                        U.append(path)
                        
                path_checking = [i, j, n+j ,n+i]
               
                U = []
                path = []
                L = ls.Label(i, data_dict["time_windows"][i][0], data_dict["load"][i], 0, [], [], P, D, -1, n, 0, 0, {})
                path.append(L)
                U.append(path)
                  
                while U:
                    path = U.pop() #U is a queue first in first out!
                    L = path[-1]
                    current_node = path[-1].id
                    if current_node == path_checking[0]:
                        max_travel_time_of_user = data_dict["max_ride_time_per_user"][i]
                        next_node = path_checking[1]
                    elif current_node == path_checking[1]:
                        max_travel_time_of_user = data_dict["max_ride_time_per_user"][j]
                        next_node = path_checking[2]
                    elif current_node == path_checking[2]:
                        max_travel_time_of_user = min(data_dict["max_ride_time_per_user"][i], data_dict["max_ride_time_per_user"][j])
                        next_node = path_checking[3]
                    else:
                        continue
                        
                    if (L.t + travel_time(g, current_node, next_node) > data_dict["time_windows"][next_node][1] or
                    L.t + travel_time(g, current_node, next_node) > max_travel_time_of_user or 
                    L.l + data_dict["load"][current_node] > data_dict["vehicle_capacity"]):
                        flag +=1
                    else:
                        time_node = L.t + travel_time(g, current_node, next_node)
                        cummulative_load = L.l + data_dict["load"][next_node]
                        L_new = ls.Label(next_node, time_node, cummulative_load, 0, [], [], P, D, current_node, n, 0, 0, {})
                        path.append(L_new)
                        U.append(path)
                if flag == 2:
                    unable_edges.add((i, j))
                    g.delete_edge(i,j)
    return g,unable_edges      

def eliminate_arc_4(g, n, P, D, data_dict, unable_edges):
    #deleting edges of the form (n+i, n+j) if the paths (i, j, n+i ,n+j) and the (j, i, n+i ,n+j) are infeasible 

    for iter1 in P:
        i= iter1
        for iter2 in P:
            j = iter2
            if (i,j) in unable_edges:
                continue
            flag = 0
            if i != j:
                U = []

                path_checking = [i, j, n+i ,n+j]
                if not g.has_edge(i,j):
                    continue
                
                path = []
                L = ls.Label(i, data_dict["time_windows"][i][0], data_dict["load"][i], 0, [], [], P, D, -1, n, 0, 0, {})
                path.append(L)
                U.append(path)
                  
                while U:
                    path = U.pop() #U is a queue first in first out!
                    L = path[-1]
                    current_node = path[-1].id
                    if current_node == path_checking[0]:
                        max_travel_time_of_user = data_dict["max_ride_time_per_user"][i]
                        next_node = path_checking[1]
                    elif current_node == path_checking[1]:
                        max_travel_time_of_user = data_dict["max_ride_time_per_user"][j]
                        next_node = path_checking[2]
                    elif current_node == path_checking[2]:
                        max_travel_time_of_user = min(data_dict["max_ride_time_per_user"][i], data_dict["max_ride_time_per_user"][j])
                        next_node = path_checking[3]
                    else:
                        continue
                        
                    if (L.t + travel_time(g, current_node, next_node) > data_dict["time_windows"][next_node][1] or
                    L.t + travel_time(g, current_node, next_node) > max_travel_time_of_user or 
                    L.l + data_dict["load"][current_node] > data_dict["vehicle_capacity"]):
                        flag +=1
                        #g.delete_edge(n+i, j)
                    else:
                        time_node = L.t + travel_time(g, current_node, next_node)
                        cummulative_load = L.l + data_dict["load"][next_node]
                        L_new = ls.Label(next_node, time_node, cummulative_load, 0, [], [], P, D, current_node, n, 0, 0, {})
                        path.append(L_new)
                        U.append(path)
                        
                path_checking = [j, i, n+i ,n+j]

                U = []
                path = []
                L = ls.Label(i, data_dict["time_windows"][i][0], data_dict["load"][i], 0, [], [], P, D, -1, n, 0, 0, {})
                path.append(L)
                U.append(path)
                  
                while U:
                    path = U.pop() #U is a queue first in first out!
                    L = path[-1]
                    current_node = path[-1].id
                    if current_node == path_checking[0]:
                        max_travel_time_of_user = data_dict["max_ride_time_per_user"][i]
                        next_node = path_checking[1]
                    elif current_node == path_checking[1]:
                        max_travel_time_of_user = data_dict["max_ride_time_per_user"][j]
                        next_node = path_checking[2]
                    elif current_node == path_checking[2]:
                        max_travel_time_of_user = max(data_dict["max_ride_time_per_user"][i], data_dict["max_ride_time_per_user"][j])
                        next_node = path_checking[3]
                    else:
                        continue
                        
                    if (L.t + travel_time(g, current_node, next_node) > data_dict["time_windows"][next_node][1] or
                    L.l + data_dict["load"][next_node] > data_dict["vehicle_capacity"]):
                        flag +=1
                    else:
                        time_node = L.t + travel_time(g, current_node, next_node)
                        cummulative_load = L.l + data_dict["load"][next_node]
                        L_new = ls.Label(next_node, time_node, cummulative_load, 0, [], [], P, D, current_node, n, 0, 0, {})
                        path.append(L_new)
                        U.append(path)
                if flag == 2:
                    unable_edges.add((n+i, n+j))
                    g.delete_edge(n+i,n+j)
    return g,unable_edges       
    
def graph_preprocessing(g, n, P, D, data_dict):
    """
    eliminates edges applying criteria mentioned in the paper
    """
    g,unable_edges = eliminate_arc_1(g, n, P, D, data_dict)
    g,unable_edges = eliminate_arc_2(g, n, P, D, data_dict, unable_edges)
    g,unable_edges = eliminate_arc_3(g, n, P, D, data_dict, unable_edges)    
    g,unable_edges = eliminate_arc_4(g, n, P, D, data_dict, unable_edges)
    return g
    
def create_time_windows_per_node(g, time_win, P, D, n,  W = 1800):
    """
    :g: sagemath's default graph currently c_graph
    :time_win: dictionary containing the time demanded for the request to be picked up
    """
    time_windows = {}
    min_time_window = time_win[1]
    for i in P:
        if time_win[i] < min_time_window:
            min_time_window = time_win[i]
        time_windows[i] = [max(time_win[i]-W,0), time_win[i] + W]
        time_windows[n+i] = [max(time_win[i]- W, 0), time_win[i] + 5 * g.edge_label(i,n+i) + W]
    
    
    time_windows[0] = [0, 86400.0]
    time_windows[2*n+1] = [0, 86400.0]

    return time_windows, min_time_window

def create_max_ridetime_per_user(g, P, n):
    """
    :g: sagemath's default graph currently c_graph
    :P: set of pickup nodes
    :n: number of requests
    max travel time is the time to go from pickup point to delivery point + 45 minutes
    """
    max_travel_time = {}
    for request in P:
        max_travel_time[request] = travel_time(g, request, request+n) + 40*60 
    return max_travel_time
     
def create_data_dict(g, P, D, time_window, n, T, W, vehicle_capacity = 5, min_load = 1, max_load = 2):
    """
    :g: sagemath default graph currently c_graph
    :P: set of pickup nodes
    :D: set of delivery nodes
    :n: num of requests
    :T: planning horizon
    :W: time window interval
    :returns
    """
    data_dict = {}
    data_dict["P"] = P
    data_dict["D"] = D
    data_dict["time_windows"], data_dict["min_time_window"] = create_time_windows_per_node(g, time_window, P, D, n)
    data_dict["load"]  = create_load(P, D, min_load, max_load)
    data_dict["num_of_requests"] = (g.order()-2) /2
    data_dict["vehicle_speed"] = 80
    data_dict["vehicle_maximum_waiting_time"] = 45*60
    data_dict["client_maximum_waiting_time"] = 20*60
    data_dict["short_planning_horizon"] = T 
    data_dict["vehicle_capacity"] = vehicle_capacity
    data_dict["time_window_interval"] = W 
    data_dict["max_ride_time_per_user"] = create_max_ridetime_per_user(g, P, n)
    
    return data_dict
    
def set_edges_cost(vertices_enc, mapping_code_to_num, n, graph):
    '''
    :vertices_encoded: set of vertices names as written in the file not as numbers
    :mapping_code_to_num: dictionary with a number as key and the corresponding node code as value
    :n: num of requests
    :graph: sagemath's default graph currently c_graph
    :returns: graph with the cost mentioned in the file as edges_weights
    '''
    counter = 0
    rnd = current_randstate().python_random()

    num_of_nodes = 2*n+1
    min_dist = -1
    with open('C:/Users/nataz/Desktop/Github_repositories/PDPTW/distanceMatrix.csv') as distances_file:
        line = distances_file.readline()
        if counter == 0:
            line = distances_file.readline().split(";")
            counter += 1
        while line != ['']: 
            if line[0] in vertices_enc and line[3] in vertices_enc:
                nodes_list_from = mapping_code_to_num[line[0]] #list of nodes corresponding to line[0] code
                nodes_list_to = mapping_code_to_num[line[3]]
                for node1 in nodes_list_from:
                    for node2 in nodes_list_to:
                        graph.set_edge_label(node1, node2, float(line[6])) 
                        
            line = distances_file.readline().split(";")
            
    for nod in graph.vertices():
        graph.set_edge_label(0, nod, 10000 ) #10000 because we want to minimize the number of vehicles to use
        graph.set_edge_label(nod, 2*n+1, 1000)
        
    graph.set_edge_label(0, 2*n+1, np.inf) 
    
    #In this part I put in the edge between two nodes that are actually the same the min_cost. Not sure that it is needed for the triangle inequality.
    weights = np.asarray(graph.edge_labels())
    min_dist = np.min(weights[np.nonzero(weights)]) 
 
    for edge in graph.edges():
        if edge[2] is None:
            graph.set_edge_label(edge[0], edge[1], min_dist)
    return graph
    
def create_load(P, D, min_load= 2, max_load = 2):
    '''
    #Introduces load constraints randomly given an interval of min, max load to consider
    :P: set of pickup nodes
    :D: set of delivery nodes
    '''
    rnd = current_randstate().python_random()
    load = [0]
    for i in P:
        load.append(rnd.randrange(min_load, max_load + 1))
    for i in D:
        load.append(-load[i-len(P)])
    load.append(0)
    return load
    
def choose_data(data, day = 20, month = "Dec", year = 2017):
    """
    :data: data in the form of dataframe
    """
    new_data =  data.loc[(data["day_num"] == day) & (data["month"] == month) &  (data["year"] == year)]
    return new_data
    
def check_triangle_inequality(g):
    """
    Check if the weights f the graph satisfy the triable inequality
    In order to check this, we replace the cost between two nodes with their shortest path and see if the costs have changed. If not then the initial costs already satisfied the triangle inequality. 
    :g: sagemath default graph currently c_graph
    :returns: True if triangle inequality is satisfied, False otherwise
    """

    #calculate all pairs shortest path 
    dist, pred = g.shortest_path_all_pairs(by_weight=True, algorithm="Dijkstra_Boost")
    edges_before = g.edges()
    #create graph G' which is guaranteed to satisfy the triangle inequality
    for edge in g.edges():
        u = edge[0]
        v = edge[1]
        w = dist[u][v]
        g.set_edge_label(u, v, w)
    edges_after = g.edges()
    for edge_i in edges_before:
        for edge_j in edges_before:
            if edge_i == edge_j and g.edge_label(edge_i[0], edge_i[1]) != g.edge_label(edge_j[0], edge_j[1]):
                print("Triangle inequality is not satisfied")
                return False
    print("Triangle inequality is satisfied!")
    return True
    
def create_graph(daily_data):
    """
    This scripts "main function" it gets the data dataframe and creates the graph and the data_dict with all the information needed to run Branch and Bound
    """
    vertices = range(1,daily_data.shape[0]+1)
    mapping_dict_num_to_code = {}
    mapping_code_to_num = {}
    P = []
    D = []
    n = daily_data.shape[0]
    i = 1
    time_windows = {}
    vertices = []
    edges = []
    edges_codes = []
    vertices_encoded = ['s']
    for index, row in daily_data.iterrows():
        mapping_dict_num_to_code[i] = row["departure_node"]
        time_windows[i] = row["seconds"]
        try:
            mapping_code_to_num[row["departure_node"]].append(i)
        except:
            mapping_code_to_num[row["departure_node"]] = [i]
            
        vertices_encoded.append(row["departure_node"])
        P.append(i)
        mapping_dict_num_to_code[i + n] = row["arrival_node"]
        
        try:
            mapping_code_to_num[row["arrival_node"]].append(n+i)
        except:
            mapping_code_to_num[row["arrival_node"]]= [n+i]
            
            
        vertices_encoded.append(row["arrival_node"])
        D.append(n+i)
        i += 1
      
    vertices = range(2*n+2)
    vertices_encoded.extend('t')
    vertices_encoded = set(vertices_encoded)
    g =  digraphs.Complete(2*n+2)
    
    #delete not needed edges
    g = delete_edges(g, P, D, n)
    g = set_edges_cost(vertices_encoded, mapping_code_to_num, n, g)
    print(check_triangle_inequality(g))
    
    return g, P, D, time_windows, mapping_dict_num_to_code
    
    

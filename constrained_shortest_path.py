import pandas as pd
import pickle
from itertools import chain
import time
from sage.all import *
import sys
rnd = current_randstate().python_random()
from operator import itemgetter
import read_data_script as read_data 
from collections import deque
import local_search_new as ls
import numpy as np

class Label:
    '''
    This class is used to construct the Labels for the Labelling algorithm used to solve the constrained shortest path problem
    '''
    def __init__(self ,number_of_node, arrival_time, load, cost, V_parent, O_parent, P, D, parent_node, n, dij_cost, path_start_time, prev_path_times, path_prev, path_loads_prev, path_reduced_costs_prev, path_costs_prev):
        self.id = number_of_node 
        self.t = arrival_time #arrival time in Label L 
        self.l = load # cummulative load at this node using the a path
        self.c = cost # the cost on this solution is the distance
        self.V = Set(V_parent) # V is the set of requests served by this path that may or may not have already been served
        self.construct_V(P, D)
        self.O = Set(O_parent) #O is the set of open requests. Requests that have started but not delivered yet
        self.construct_O(P, D, n)
        self.p = parent_node # previous node in the path
        self.d = dij_cost
        self.path_times = prev_path_times
        self.store_times()
        self.path_starting_time = path_start_time
        self.set_path_starting_time()
        self.path_ids = path_prev
        self.path_ids.append(number_of_node)
        self.path_loads = path_loads_prev
        self.path_reduced_costs = path_reduced_costs_prev
        self.path_costs = path_costs_prev
        self.path_loads[self.id] = self.l
        self.path_reduced_costs[self.id] =self.d
        self.path_costs[self.id] = self.c  

    def deepcopy(label_obj):
        self.id = label_obj.id
        self.t = label_obj.t
        self.l = label_obj.l
        self.c = label_obj.c
        self.V = copy(label_obj.V)# V is the set of requests served by this path that may or may not have already been served
        self.O = copy(label_obj.O)
        self.p = parent_node # previous node in the path
        self.d = dij_cost
        self.path_times = copy(label_obj.path_times)
        self.path_starting_time = copy(label_obj.path_starting_time)
        self.path_ids = copy(label_obj.path_ids)
        self.path_loads = copy(label_obj.path_loads)
        self.path_reduced_costs = copy(label_obj.path_reduced_costs)
        self.path_costs = copy(label_obj.path_costs)
        self.path_loads[self.id] = copy(label_obj.path_loads[self.id])
        self.path_reduced_costs[self.id] = copy(label_obj.path_reduced_costs[self.id])
        self.path_costs[self.id] = copy(label_obj.path_costs[self.id])
        return 
            
    def set_path_starting_time(self):
        if self.p == 0:
            self.path_starting_time = self.t

    def print_data(self):
        print("id",self.id)
        print("arrival time",self.t)
        print("load",self.l)
        print("cost",self.c)
        print("self.O",self.O)
        print("parent",self.p)
        print("dij",self.d)
        print("self.path_times",self.path_times)
        print("self.set_path_starting_time()",self.path_starting_time)
        print("self.path_ids",self.path_ids) 
        
    def store_times(self):
        self.path_times[self.id] = self.t
        
    def construct_V(self, P, D):
        if self.id in P:
            self.V = self.V.union(Set([self.id]))
        
    
    def construct_O(self, P, D, n):

        if self.id in P:
            self.O = self.O.union(Set([self.id]))
        elif (self.id in D): #and self.id-n in self.O:
            self.O = self.O.difference(Set([self.id-n]))
           
class PricingSubproblems():
    def __init__(self, g, data_dict, max_time_in_veh = 1800):
        self.g = g
        self.n = data_dict["num_of_requests"]
        self.load_per_req = data_dict["load"]
        self.time_window = data_dict["time_windows"]
        self.speed = data_dict["vehicle_speed"]
        self.P = data_dict["P"]
        self.D = data_dict["D"]
        self.vehicle_capacity = data_dict["vehicle_capacity"]
        #self.veh_max_waiting_time = data_dict["vehicle_maximum_waiting_time"]
        #self.passenger_max_waiting_time = data_dict["client_maximum_waiting_time"]
        self.planning_horizon  = data_dict["short_planning_horizon"]
        self.min_time_win = data_dict["min_time_window"]
        self.max_time_in_vehicle = max_time_in_veh
        self.num_of = 0 
        self.reduced_cost_dict = {}
        self.transformed_reduced_cost = {}
        self.theta ={}
        self.data_dict = data_dict
    
    def travel_time(self, from_node, to_node):
        """
        Gets the travel times between two locations.
        :from_node: int
        :to_node: int
        :Return: float corresponding to seconds
        """
 
        travel_time = self.g.edge_label(from_node,to_node)
        return travel_time 
    
    def get_path_lns(self, path):
        """
        Takes a path in Label form and returns the path in node ids, the cost in distance and the reduced_cost. 
        :path: Label 
        :Return: list of nodes in path, cost(float), reduced cost(float)
        """
        path_id_list = []
        cost_list = []
        red_cost = []
        if path is None:
            return None
            
        for node in path:
            path_id_list.append(node.id)
            cost_list.append(node.c)
            red_cost.append(node.d)
        
        return path_id_list, cost_list, red_cost
    
    def get_path(self, path):
        """
        Takes a path in Label form and returns the path in node ids, the cost in distance and the reduced_cost. 
        :path: Label 
        :Return: list of nodes in path, cost(float), reduced cost(float)
        """
        if path is None:
            return None
        return path.path_ids, path.c, path.d
        
    def get_paths(self, result):
        """
        result: list of paths with Labels
        Return: the corresponding path in node_id, costs in destance and reduced cost
        """
        all_paths = []
        costs = []
        d_costs = []
        for path in result:
            all_paths.append(path.path_ids)
            costs.append(path.c)
            d_costs.append(path.d)
        return all_paths, costs, d_costs
    
    def convert_path_to_labels(self, paths_list, g):
        #used when applying local search
        paths_lab_list = []
        for path in paths_list:
            r_num = 0
            U = []
            s = 0
            cost = 0
            path_lab = []
            V_parent = []
            O_parent = [] 
            parent_node = -1
            source_time = self.time_window[s][0]
            prev_path_times ={}
            L = Label(s, source_time, self.load_per_req[s], cost, V_parent, O_parent, self.P, self.D, parent_node, self.n, 0, self.min_time_win, prev_path_times, [], {}, {}, {}) 
            path_lab.append(L)
            prev_node = s
            for i in range(len(path)-1):
                current_node = path[i+1]
                cost = self.g.edge_label(prev_node,current_node ) + L.c #cummulative cost
                dij_cost =  self.g.edge_label(prev_node,current_node ) + L.c
                t_j = self.travel_time(prev_node, current_node) + L.t
                cummulative_load = L.l + self.load_per_req[current_node]
                L_new = Label(current_node, t_j, cummulative_load, cost, L.V, L.O, self.P, self.D, prev_node, self.n, dij_cost, self.min_time_win, copy(L.path_times), copy(L.path_ids),copy(L.path_loads), copy(L.path_reduced_costs), copy(L.path_costs))
                path_lab.append(L_new)
                prev_node = current_node
                L = L_new
            paths_lab_list.append(path_lab[-1])
                
        return paths_lab_list
    
    
    def get_random_nodes(self,nodes_sorted, num_of_neighbors_to_expand, weight_low_ranks = 1.0, weight_high_ranks = 2.0):
        #randomly choose nodes from the list "nodes _sorted" by assigning probabilities on the first and the second half of the list
        #according to the inputs weight_low_ranks, weight_high_ranks. 
        if len(nodes_sorted) <= num_of_neighbors_to_expand:
            num_of_neighbors_to_expand = len(nodes_sorted)
        prob =  [weight_high_ranks]*(len(nodes_sorted)//2) + [weight_low_ranks]*(len(nodes_sorted) - len(nodes_sorted)//2 )
        prob /= np.sum(prob)
        number = np.random.choice(nodes_sorted, num_of_neighbors_to_expand-1, p=prob,replace = False)
        return list(number)   
        
    def check_dominance(self, L, labels_ending_at):
        '''
        To use this function triangle inequality should be satisfied both on distances and times.
        In the paper the dominance rule refers to the cost and time constrained. After some experimentation on introducing the capacity constraints 
        we concluded to the dominance rule below. The capacity constrained "l" can not be set as <= cause then suboptimal solutions are found. Thus we introduce it  as <.
        U: list of paths 
        L: Label of the node to be added on the path
        :labels_ending_at: dictionary containing the node id as key and as value the list of Labels ending in that node.
        :return: True if there is a path that dominates L, False otherwise
        '''
        checking_node = L.id
        if checking_node == 2*self.n +1:
            return False
        for L_i in labels_ending_at[checking_node]:
            if L_i.d <= L.d and L_i.V.issubset(L.V)  and L_i.O.issubset(L.O)  and L_i.t <= L.t and L_i.l < L.l:
                self.num_of += 1
                return True    
        return False
    
    def calculate_theta(self, node1, node2, duals):
        for j in self.P:
            max_theta = -1000000
            for k in range(1,2 * self.n +1):
                if j == k:
                    continue
                
                for i in range(1,2 * self.n +1):
                    if i == k or i==j  or i == self.n +j or k == self.n +j:
                        continue
                    try:
                        red_cost_1= self.calculate_reduced_cost_new(i, k, duals)
                        red_cost_2 = self.calculate_reduced_cost_new(i, self.n + j, duals)
                        red_cost_3 = self.calculate_reduced_cost_new(self.n+j, k, duals)
                        temp =  red_cost_1 - (red_cost_2 + red_cost_3)
                    except:
                        continue
                    if temp > max_theta:
                        max_theta = temp
                        self.theta[j] =  temp
                        
    def calculate_transformed_reduced_cost(self, node1, node2, duals):
        #calculates transformed reduced cost for edges starting from pickup nodes
        if self.theta == {}:
            self.calculate_theta(node1, node2, duals)
        
        if node1 in self.P:
            self.transformed_reduced_cost[(node1,node2)] = self.calculate_reduced_cost(node1,node2, duals) - self.theta[node1]
            return self.transformed_reduced_cost[(node1,node2)]
        else:
            self.transformed_reduced_cost[(node1,node2)] = self.calculate_reduced_cost(node1,node2, duals) + self.theta[node1- self.n]
            return self.transformed_reduced_cost[(node1,node2)]
    
    def calculate_reduced_cost_new(self, node1, node2, duals):
        mi = []
        if len(duals) > 2*self.n +2:
            for i in range(2*self.n + 2, len(duals)):
                mi.append(duals[i])
        else:
            mi.append(0.0)

        if node1 in self.P:
            dij = self.g.edge_label(node1, node2) - duals[node1]
            for i in mi:
                dij = dij - i
            self.reduced_cost_dict[(node1, node2)] = dij
        else:
            dij = self.g.edge_label(node1, node2)
            for i in mi:
                dij = dij -  i
            self.reduced_cost_dict[(node1, node2)] = dij
        return dij
    
    
    def calculate_reduced_cost(self, node1, node2, duals):
        
        if len(duals) > 2*self.n +2:
            
            mi = duals[2*self.n +2]
        else:
            mi = 0.0

        if node1 in self.P:
            dij = self.g.edge_label(node1, node2) - duals[node1] - mi
            self.reduced_cost_dict[(node1, node2)] = dij
        else:
            dij = self.g.edge_label(node1, node2) -  mi
            self.reduced_cost_dict[(node1, node2)] = dij
        return dij
     
     
    def get_reduced_cost(self, node1, node2, duals):
        if node1 == 0:
            try:
                return self.reduced_cost_dict[(node1, node2)]
            except:
                return self.calculate_reduced_cost_new(node1, node2, duals)
        else:
            try:
                return self.transformed_reduced_cost[(node1,node2)]
            except:
                return self.calculate_transformed_reduced_cost(node1, node2, duals)
            
    """   
    def get_reduced_cost(self, node1, node2, duals):
        '''
        Calculates the reduced cost using the rule on the paper we are implementing
        :node1: int
        :node2: int
        :duals: list of floats(dual variables)
        :Return: float, reduced cost from node1 to node2
        '''
        #In case no additional constraints have been introduced calculate reduced cost
        if len(duals) <= 2*self.n +2:
            try:
                return self.reduced_cost_dict[(node1,node2)]
            except KeyError:
                dij = self.calculate_reduced_cost(node1, node2, duals)
                dij_new = self.calculate_reduced_cost_new(node1, node2, duals)
                
                return dij_new
        else: #this is the case where we have to transform the reduced cost for the triangle inequality to hold only for the pickup nodes though
            if node1 in self.D or node1 == 0 or node1 == 2*self.n +1:
                try:
                    return self.reduced_cost_dict[(node1,node2)]
                except KeyError:
                    #dij = self.calculate_reduced_cost(node1, node2, duals)
                    dij_new = self.calculate_reduced_cost_new(node1, node2, duals)
                    return dij_new
            else:
                try:
                    return self.transformed_reduced_cost[(node1, node2)]
                except KeyError:
                    dij  = self.calculate_transformed_reduced_cost(node1, node2, duals)
                    return dij 
    """    
                   
        
    def deepcopy(self, obj_list):
        """
        Deep copy of a list of objects
        obj_list: list of objects
        :Returns: a list of the same objects with different reference on the database
        """
        new_list = []
        for obj in obj_list:
            new_obj = copy(obj)
            new_list.append(new_obj)
        return new_list
         
            
    def label_elimination(self, L):
        #implementation of the label elimination idea proposed in the paper. It doesn't apply in our data and there for it is not used.
        open_requests = list(L.O)
        result = []
        open_deliveries = set([i+self.n for i in open_requests])
        U = deque()
        U.append([L])
        
        while U:
            nodes_in_path = []
            path = U.popleft() #U is a queue first in first out!
            
            if len(path) == len(open_deliveries) +1:
                return True
                
            for p in path:
                nodes_in_path.append(p.id)
            
            
            L = path[-1]
            for delivery_node in open_deliveries:
                if L.id ==  delivery_node:
                    continue

                time_node = L.t + self.travel_time(L.id, delivery_node)
                try:
                    path_starting_time = path[1].t
                except:
                    path_starting_time = path[0].t
                if time_node <= self.time_window[delivery_node][1] and delivery_node not in nodes_in_path and time_node <= self.planning_horizon + path_starting_time:
                    L_new = Label(delivery_node, time_node, 0, 0, L.V, L.O, self.P, self.D, i, self.n, 0, self.min_time_win, copy(L.path_times),copy(L.path_ids),copy(L.path_loads), copy(L.path_reduced_costs), copy(L.path_costs))
                    new_path = list(path)
                    new_path.append(L_new)
                    U.append(new_path)
                    #result.append(new_path)

        return False
              
    def get_sorted_list_on_reduced_cost(self, node, duals, open_deliveries):
        """
        Return a sorted list of the incident nodes of node
        :node: int,  node id of which we sort its neighbours
        :duals: list of dual variables
        :open_deliveries: list of ints corresponding to delivery node_ids
        """

        neighbors = []
        neighbors_red_cost = []
        temp = []
        delivery_nodes = []
        for inc_node in self.g.edges_incident(node):
            if inc_node[1] in self.D and inc_node[1] in open_deliveries :
                delivery_nodes.append(inc_node[1])
            elif inc_node[1] == 2*self.n+1:
                continue
            elif inc_node[1] in self.P:
                neighbors.append(inc_node[1])
                neighbors_red_cost.append(self.get_reduced_cost(node, inc_node[1], duals)) #maybe there will be problem.
        if neighbors_red_cost != []: 
            temp = zip(neighbors_red_cost, neighbors)
            temp.sort()
        return temp, delivery_nodes
    
    def calculate_time_node(self, node1, node2, L):
        """
        Calculates time from node1 to node2 and returns the time after adding on the path corresponding to label L the time from node1 to node2
        :node1: int
        :node2: int
        :L: Label object 
        """
        if  node1 != 0:
            time_node = self.travel_time(node1, node2) + L.t
        else:
            time_node = self.time_window[node2][0]
        return time_node
    
    
    def initialize_U_all_requests(self, duals, labels_ending_at):
        """
        Create first labels corresponding to paths of the type (0, i) for every request i
        :duals: list of dual variables
        :labels_ending_at: dictionary with keys the id of the last node of the path and values a list of Labels(therefore paths) found in previous heuristics that failed to find a negative reduced cost path
        return: list containing labels
        """
        depot_id = 0 
        cost = 0
        U = deque() 
        dij_cost = 0
        parent_node = -1
        source_time = self.min_time_win  #time of the first request made on the set of requests that we have
        L = Label(depot_id, source_time, self.load_per_req[depot_id], cost, [], [], self.P, self.D, parent_node, self.n, dij_cost, self.min_time_win, {}, [], {}, {}, {}) 
        U.append(L)
        for req in self.P:
            i = L.id
            time_node = self.calculate_time_node(i, req, L)
            cummulative_load = L.l + self.load_per_req[req] 
            U, labels_ending_at = self.expand_path(U, cummulative_load, time_node, labels_ending_at, L, i, req,  duals)
         
        return U
    
    def initialize_U(self, duals):
        """
        Create first label corresponding to a path including only the starting node.
        :duals: list of dual variables
        return: List containing one label
        """
        depot_id = 0 
        cost = 0
        U = deque() 
        dij_cost = 0
        parent_node = -1
        source_time = self.min_time_win  #time of the first request made on the set of requests that we have
        L = Label(depot_id, source_time, self.load_per_req[depot_id], cost, [], [], self.P, self.D, parent_node, self.n, dij_cost, self.min_time_win, {}, [],{},{},{}) 
        U.append(L)
        return U
           
    def expand_path(self, U, cummulative_load, time_node, labels_ending_at, L, i, current_node,  duals):
        """
        Check the sp1 algorithm condition and the dominance criteria, expand the path and store the label on the labels_ending_at dict 
        :U: list of paths to where the new label will be added
        :cummulative_load: int load on the current_node.
        :time_node:float, time at the current_node
        :labels_ending_at: dictionary with keys the id of the last node of the path and values a list of Labels(therefore paths) found in previous heuristics that failed to find a negative reduced cost path
        :L: Label object of the path we are trying to expand
        :i: int, previous node_id on the path
        :current_node: int, node on which we are expanding the path 
        :duals: list of dual variables
        """
        cost = self.g.edge_label(i,current_node) + L.c 
      
        dij = self.get_reduced_cost(i, current_node, duals)
        dij_cost = dij + L.d
        
        if 0 < current_node and current_node <= self.n and (current_node not in L.V):
            L_new = Label(current_node, time_node, cummulative_load, cost, L.V, L.O, self.P, self.D, i, self.n, dij_cost, L.path_starting_time, copy(L.path_times), copy(L.path_ids),copy(L.path_loads), copy(L.path_reduced_costs), copy(L.path_costs)) 
            if self.check_dominance(L_new,labels_ending_at) == False:
                labels_ending_at[current_node].append(L_new)
                U.append(L_new)
        elif self.n < current_node and current_node <= 2*self.n+1 and ((current_node-self.n) in L.O):
            L_new = Label(current_node, time_node, cummulative_load, cost, L.V, L.O, self.P, self.D, i, self.n, dij_cost, L.path_starting_time, copy(L.path_times), copy(L.path_ids), copy(L.path_loads), copy(L.path_reduced_costs), copy(L.path_costs))
            if self.check_dominance(L_new, labels_ending_at) == False:
                labels_ending_at[current_node].append(L_new)
                U.append(L_new)
        elif current_node == 2*self.n+1 and L.O.is_empty():
            L_new = Label(current_node, time_node, cummulative_load, cost, L.V, L.O, self.P, self.D, i, self.n, dij_cost, L.path_starting_time, copy(L.path_times), copy(L.path_ids),copy(L.path_loads), copy(L.path_reduced_costs), copy(L.path_costs))
            if self.check_dominance(L_new, labels_ending_at) == False:
                labels_ending_at[current_node].append(L_new)
                U.append(L_new)    
        return U, labels_ending_at
    
    def expand_path_no_dom(self, U, cummulative_load, time_node, labels_ending_at, L, i, current_node,  duals):
        """
        Check the sp1 algorithm conditions and if there is no violation expand the path and store the label on the labels_ending_at dict 
        :U: list of paths to where the new label will be added
        :cummulative_load: int load on the current_node.
        :time_node:float, time at the current_node
        :labels_ending_at: dictionary with keys the id of the last node of the path and values a list of Labels(therefore paths) found in previous heuristics that failed to find a negative reduced cost path
        :L: Label object of the path we are trying to expand
        :i: int, previous node_id on the path
        :current_node: int, node on which we are expanding the path 
        :duals: list of dual variables
        """
        cost = self.g.edge_label(i,current_node) + L.c 
        dij = self.get_reduced_cost(i, current_node, duals)
        dij_cost = dij + L.d
        
        if 0 < current_node and current_node <= self.n and (current_node not in L.V):
            L_new = Label(current_node, time_node, cummulative_load, cost, L.V, L.O, self.P, self.D, i, self.n, dij_cost, L.path_starting_time, copy(L.path_times), copy(L.path_ids), copy(L.path_loads),  copy(L.path_reduced_costs), copy(L.path_costs))   
            labels_ending_at[current_node].append(L_new)
            U.append(L_new)
        elif self.n < current_node and current_node <= 2*self.n+1 and ((current_node-self.n) in L.O):
            L_new = Label(current_node, time_node, cummulative_load, cost, L.V, L.O, self.P, self.D, i, self.n, dij_cost, L.path_starting_time, copy(L.path_times), copy(L.path_ids),copy(L.path_loads),  copy(L.path_reduced_costs), copy(L.path_costs)) 
            labels_ending_at[current_node].append(L_new)
            U.append(L_new)
        elif current_node == 2*self.n+1 and L.O.is_empty():
            L_new = Label(current_node, time_node, cummulative_load, cost, L.V, L.O, self.P, self.D, i, self.n, dij_cost, L.path_starting_time, copy(L.path_times), copy(L.path_ids),copy(L.path_loads),  copy(L.path_reduced_costs), copy(L.path_costs)) 
            labels_ending_at[current_node].append(L_new)
            U.append(L_new)    
        return U, labels_ending_at
    
    def calculate_H1(self, discovered_paths, duals, num_of_neighbors_to_expand, labels_ending_at, weight_low_ranks= 1.0, weight_high_ranks= 3, num_paths = 10, timeout = 5):
        """
        H1 does the same as SP1 but we expand the paths only in "num_od_neighbors_to_expand" nodes.
        :discovered_paths: dictionary of paths in  id form that are already in the solution space (Omega) of the master problem.
                           The keys are the frozensets of nodes in the path and the values a list of paths containing these nodes. Dictionary form is used to accelerate search.
        :duals: list of dual variables obtained by the linear programming problem 
        :labels_ending_at: dictionary with keys the id of the last node of the path and values a list of Labels(therefore paths) found in previous heuristics that failed to find a negative reduced cost path
        :num_of_paths: maximum number of paths to return
        :Return: list of labels with negative reduced cost( at most paths_per_itter labels), labels_ending_at dict to be used in next heuristic if result list was empty.
        """
        #Initialize variables 
        starting_time = time.time()
        result = []
        #self.reduced_cost_dict = {}
        #self.transformed_reduced_cost = {}
        U = self.initialize_U(duals)
        
        
        #Introduce the paths found from the previous heuristics in this iteration U to expand, excepting the paths ending in the target node. 
        if labels_ending_at == {}:
            for vert in self.g.vertices():
                labels_ending_at[vert] = []
        else:
            for k,v in labels_ending_at.items():
                if   k ==2*self.n + 1: 
                    labels_ending_at[k] = []
                for l in v:
                    if l != 2*self.n + 1:
                        U.append(l)
              
        while U: 
            current_time = time.time()
            L = U.pop()
            i = L.id
            open_deliveries = set([j+self.n for j in list(L.O)])
            
            if current_time - starting_time >= 2000: #if h1 runs more than 2000 seconds return and run sp1 because probably the optimal solution in the master problem is already found is already found.
                return result, labels_ending_at
            
            #returning_condition: if h1 exceeds 30 seconds and at least one negative cost path is found return what is already found.
            if len(result) > 0 and current_time - starting_time > timeout:
                return result, labels_ending_at
                    
            #return check
            if i == 2*self.n+1:
                if L.d < 0: 
                    current_path , temp_cost , temp_dcost = L.path_ids, L.c, L.d
                    try:
                        temp_paths = discovered_paths[frozenset(current_path)]
                    except:
                        result.append(L)
                        if len(result) > num_paths: 
                            return result, labels_ending_at
                        continue
                    if current_path not in temp_paths:
                        result.append(L)
                        if len(result) >= num_paths:
                            return result, labels_ending_at
                        continue
                else:
                    continue
                    
            if self.g.edges_incident(i) != []:
                #choosing num_of_neighbors_to_expand nodes by sorting the list of neighbors on reduced cost
                pickup_inc = []
                delivery_inc = []
                request_nodes, delivery_nodes = self.get_sorted_list_on_reduced_cost(i, duals, open_deliveries)
                
                if i > self.n and i != 2* self.n +1:
                    delivery_nodes.append(2* self.n +1)
                    
                if request_nodes == []:
                    inc_nodes  = delivery_nodes
                else:
                    distances, incident_nodes_sorted = zip(*request_nodes)
                    request_nodes_chosen = self.get_random_nodes(incident_nodes_sorted, num_of_neighbors_to_expand, weight_low_ranks, weight_high_ranks)
                    if delivery_nodes ==[]:
                        inc_nodes  = request_nodes_chosen
                    else:
                        inc_nodes  = delivery_nodes + request_nodes_chosen
                        shuffle(inc_nodes)
                #expand current path towards selected incident_nodes
                for j in range(len(inc_nodes)): #get all incident nodes from i
                    current_node = inc_nodes[j]
                    if current_node in L.path_times: #skip if neighbor node is already in the path 
                        continue
                
                    if current_node in self.D and current_node not in open_deliveries: #skip if current code is in delivery nodes but the request node hasn't been visited
                        continue
                    elif current_node in self.P and L.l == self.vehicle_capacity:
                        continue
                    
                
                    time_node = self.calculate_time_node(i, current_node, L)
                    cummulative_load = L.l + self.load_per_req[current_node] 
                    
                    if (cummulative_load <= self.vehicle_capacity and time_node <= self.time_window[current_node][1] and \
                    time_node >= self.time_window[current_node][0] and time_node <= self.planning_horizon + L.path_starting_time):
                        U, labels_ending_at = self.expand_path(U, cummulative_load, time_node, labels_ending_at, L, i, current_node, duals)
                
        return result, labels_ending_at 

    def calculate_H3(self, discovered_paths, duals, use_all_requests = False, num_of_neighbors_to_expand = 2):
        """
        Starting from a route containing only [0,i] for every request i we expand the path choosing the num_of_neighbors_to_expand incident nodes with the lowest reduced cost
        :discovered_paths: dictionary storing paths already in master problem set of paths with keys a frozenset of the nodes in the paths.
        :duals: list of dual variables
        :use_all_requests: used to discriminate between H3 and H3_all. H3_all starts with the sets [0,i) for all requests i while H3  starts from (0,j) where j is randomly chosen between the requests. 
        :num_of_neighbors_to_expand: num of neighbors to expand in each iteration
        :Return: list of labels with negative reduced cost( at most paths_per_itter labels), labels_ending_at dict to be used in next heuristic if result list was empty.
        """
        #Initialize variables
        starting_time = time.time()
        paths_per_itter = 8
        #self.reduced_cost_dict = {}
        #self.transformed_reduced_cost = {}
        result = []

        #initialize dictionary storing the Labels using the node_id on the path as key 
        labels_ending_at = {}
        for vert in self.g.vertices():
            labels_ending_at[vert] = []
        
        #initialize U depending on whether this is called by H3 or H3_all
        if use_all_requests:
            U = self.initialize_U_all_requests(duals, labels_ending_at)
        else:
            U = self. initialize_U(duals)
            
        while U: 
            current_time = time.time()
            L = U.popleft()
            i = L.id
            
            open_deliveries = set([j + self.n for j in list(L.O)])

            #return check
            if i == 2*self.n+1:
                if len(result) >0 and current_time - starting_time > 5:
                    return result, labels_ending_at
                
                if L.d < 0: 
                    current_path , temp_cost , temp_dcost = self.get_path(L)
                    try:
                        temp_paths = discovered_paths[frozenset(current_path)]
                    except:
                        result.append(L)
                        if len(result) >= paths_per_itter:
                            return result, labels_ending_at
                        else:
                            continue
                        return path
                    if current_path not in temp_paths:
                        result.append(L)
                        if len(result) > paths_per_itter:
                            return result, labels_ending_at
                        else:
                            continue
                        return result
                else:
                    continue
                    
            if self.g.edges_incident(i) != []:
                pickup_inc = []
                delivery_nodes = []
                incident_nodes, delivery_nodes = self.get_sorted_list_on_reduced_cost(i, duals, open_deliveries)
                if incident_nodes != []: 
                    distances, incident_nodes_sorted = zip(*incident_nodes)
                     
                    for iterator1 in range(len(incident_nodes_sorted)):
                    
                        time_node = self.calculate_time_node(i, incident_nodes_sorted[iterator1], L)
                        
                        if time_node <= self.time_window[incident_nodes_sorted[iterator1]][1] and incident_nodes_sorted[iterator1] in self.P and len(pickup_inc) <= num_of_neighbors_to_expand and incident_nodes_sorted[iterator1] not in L.path_times:
                            pickup_inc.append(incident_nodes_sorted[iterator1])                   
                            if len(pickup_inc)== num_of_neighbors_to_expand:
                                    break
                                   
                    if i > self.n and i != 2* self.n +1:
                        delivery_nodes.append(2* self.n +1)
                    inc_nodes  = delivery_nodes + pickup_inc 
                    shuffle(inc_nodes)
                else:
                    if i > self.n and i != 2* self.n +1:
                        delivery_nodes.append(2* self.n +1)
                    inc_nodes  = delivery_nodes
                
                #expand current path towards selected incident_nodes
                for j in range(len(inc_nodes)): #get all incident nodes from i
                    current_node = inc_nodes[j]
                        
                    if current_node in L.path_times: #skip if neighbor node is already in the path 
                        continue
                
                    if current_node in self.D and current_node not in open_deliveries: #skip if current code is in delivery nodes but the request node hasn't been visited
                        continue
                    elif current_node in self.P and L.l == self.vehicle_capacity:
                        continue
                        
                    time_node = self.calculate_time_node(i, current_node, L)
                    cummulative_load = L.l + self.load_per_req[current_node] 
                   
                    if (cummulative_load <= self.vehicle_capacity and time_node <= self.time_window[current_node][1] and \
                    time_node >= self.time_window[current_node][0] and time_node <= self.planning_horizon + L.path_starting_time):
                        U, labels_ending_at = self.expand_path(U, cummulative_load, time_node, labels_ending_at, L, i, current_node, duals)
               
        return result, labels_ending_at  
    
    
    def calculate_SP1(self, discovered_paths, duals, labels_ending_at, timeout = 3600):
        """
        Sp1 is the elementary path problem used to find negative reduced cost paths to add to the master problem
        :discovered_paths: dictionary of paths in  id form that are already in the solution space (Omega) of the master problem.
                           The keys are the frozensets of nodes in the path and the values a list of paths containing these nodes. Dictionary form is used to accelerate search.
        :duals: list of dual variables obtained by the linear programming problem 
        :labels_ending_at: dictionary with keys the id of the last node of the path and values a list of Labels(therefore paths) found in previous heuristics that failed to find a negative reduced cost path
        :Rerutn: a path with negative reduced cost or [] if none was found
        """
        #Initialize data structures
        result = [] 
        #self.reduced_cost_dict = {}
        #self.transformed_reduced_cost = {}
        self.num_of = 0
        U = self.initialize_U(duals)
        sp1_start_time = time.time()
        #Introduce the paths found from the previous heuristics in this iteration U to expand, excepting the paths ending in the target node. 
        if labels_ending_at == {}:
            for vert in self.g.vertices():
                labels_ending_at[vert] = []
        else:
            for k,v in labels_ending_at.items():
                
                if  k == 2*self.n + 1: # k >= self.n//2 or
                    labels_ending_at[k] = []
                
                for l in v:
                    if l != 2*self.n + 1:
                        U.append(l)
              
        #SP1 main loop
        while U: 
            L = U.pop()
            i = L.id
            sp1_current_time = time.time()
            
            #if running time exceeds timeout return [] and thus settle with a suboptimal solution
            if sp1_current_time - sp1_start_time:
                return result
                
            #check if negative reduced cost path is found
            if i == 2*self.n+1:
                if L.d < 0: 
                    current_path , temp_cost , temp_dcost = self.get_path(L)
                    try:
                        temp_paths = discovered_paths[frozenset(current_path)]
                    except:
                        return L
                    if current_path not in temp_paths:
                        return L
                else:
                    continue
                
            
            #We only need to check if the path can be expanded to delivery nodes of requests already on the set L.O(the set of request started on the path but not finished)
            open_deliveries = set([j+self.n for j in list(L.O)])
            
            #expand path on every incident edge  
            for current_node in self.g.neighbor_out_iterator(i):
                if current_node in L.path_times: #skip if neighbor node is already in the path 
                    continue
                
                if current_node in self.D and current_node not in open_deliveries:#skip if current code is in delivery nodes but the request node hasn't been visited
                    continue
                
                
                time_node = self.calculate_time_node(i, current_node, L) #calculates the arrival time on the current_node
                cummulative_load = L.l + self.load_per_req[current_node] #calculates the load of the vehicle after visiting current_node
                
                #checking capacity and time constraints
                if (cummulative_load <= self.vehicle_capacity and time_node <= self.time_window[current_node][1] and 
                    time_node >= self.time_window[current_node][0] and time_node <= self.planning_horizon + L.path_starting_time):
                    U, labels_ending_at = self.expand_path(U, cummulative_load, time_node, labels_ending_at, L, i, current_node, duals)
        
        return result
        
    def calculate_SP1_no_dom(self, discovered_paths, duals, labels_ending_at):
        """
        Sp1 is the elementary path problem used to find negative reduced cost paths to add to the master problem. 
        This implementation doesn't use dominance criteria and we only use it for testing whether the complete sp1 algorithm finds the optimum solution.
        """
        #Initialize data structures
        result = [] 
        #self.reduced_cost_dict = {}
        #self.transformed_reduced_cost = {}
        self.num_of = 0
        U = self.initialize_U(duals)
        
        #Introduce the paths found from the previous heuristics in this iteration U to expand, excepting the paths ending in the target node. 
        if labels_ending_at == {}:
            for vert in self.g.vertices():
                labels_ending_at[vert] = []
        else:
            for k,v in labels_ending_at.items():
                if  k == 2*self.n + 1: # k >= self.n//2 or
                    labels_ending_at[k] = []
                
                for l in v:
                    if l != 2*self.n + 1:
                        U.append(l)
              
        #SP1 main loop
        while U: 
            L = U.pop()
            i = L.id
            
            
            #check if negative reduced cost path is found
            if i == 2*self.n+1:
                if L.d < 0: 
                    current_path , temp_cost , temp_dcost = self.get_path(L)
                    try:
                        temp_paths = discovered_paths[frozenset(current_path)]
                    except:
                        return L
                    if current_path not in temp_paths:
                        return L
                else:
                    continue
                
            
            #We only need to check if the path can be expanded to delivery nodes of requests already on the set L.O(the set of request started on the path but not finished)
            open_deliveries = set([j+self.n for j in list(L.O)])
            
            #expand path on every incident edge  
            for current_node in self.g.neighbor_out_iterator(i):
                
                if current_node in L.path_times: #skip if neighbor node is already in the path 
                    continue
                
                if current_node in self.D and current_node not in open_deliveries:#skip if current code is in delivery nodes but the request node hasn't been visited
                    continue
                
                
                time_node = self.calculate_time_node(i, current_node, L) #calculates the arrival time on the current_node
                cummulative_load = L.l + self.load_per_req[current_node] #calculates the load of the vehicle after visiting current_node
                
                #checking capacity and time constraints
                if (cummulative_load <= self.vehicle_capacity and time_node <= self.time_window[current_node][1] and 
                    time_node >= self.time_window[current_node][0] and time_node <= self.planning_horizon + L.path_starting_time):
                    U, labels_ending_at = self.expand_path_no_dom(U, cummulative_load, time_node, labels_ending_at, L, i, current_node, duals)
        
        return result
    

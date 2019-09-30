from sage.all import *
import time
import numpy as np
from constrained_shortest_path_changed import Label, PricingSubproblems

class LNS(PricingSubproblems):
    
    def deepcopy(self, obj_list):
  
        new_obj = Label.deepcopy(obj)
      
        return new_obj
    
    
    def removeNodes(self, path, available_requests_to_remove, duals):#ok
        if available_requests_to_remove ==[]:
            return None, None, None, None
            
        rnd = current_randstate().python_random()
        random_index = rnd.randrange(0,len(available_requests_to_remove))
        i = available_requests_to_remove[random_index]
        for node_ind in range(len(path.path_ids)):
            if path.path_ids[node_ind] == i:
                pickup_removed_index = node_ind
                path.path_ids= copy(path.path_ids)
                L_id = path.path_ids.pop(node_ind)
                request_removed = L_id
                path.path_loads = copy(path.path_loads)
                path.path_loads.pop(L_id, None)
                path.path_times = copy(path.path_times)
                path.path_times.pop(L_id,None)
                
                path.path_costs = copy(path.path_costs)
                path.path_costs.pop(L_id,None)
                
                path.path_reduced_costs = copy(path.path_reduced_costs)
                path.path_reduced_costs.pop(L_id,None)
                path.V = copy(path.V)
                path.V.difference(Set([L_id]))
                
                path.O = copy(path.O)
                path.O.difference(Set([L_id]))
                
            if path.path_ids[node_ind]== self.n+i:
                delivery_removed_index = node_ind
                L_id = path.path_ids.pop(node_ind)
                path.V.difference(Set([L_id]))
                path.O.difference(Set([L_id]))
                path.path_loads.pop(L_id, None)
                path.path_times.pop(L_id,None)
                path.path_costs.pop(L_id,None)
                path.path_reduced_costs.pop(L_id,None)
                break
                
        return path, i, pickup_removed_index, delivery_removed_index
    
    def insert_delivery(self, path, delivery_inds_available, delivery_node_insert, paths_visited, duals,  set_of_paths_hash_table):
        #we will insert the delivery node in every possible position of the path so that the corresponding path is the min cost and is not in the paths visited
        rnd = current_randstate().python_random()
        min_cost_path = np.inf
        final_path = None
        final_ind  = None

        while delivery_inds_available != []:
            random_position = rnd.randrange(0, len(delivery_inds_available))
            delivery_index = delivery_inds_available[random_position]
 
            

            current_node = path.path_ids[delivery_index-1]
            node = path.path_ids[delivery_index]
            time_node = self.calculate_time_node(current_node, node, path) #calculates the arrival time on the current_node

            
            cummulative_load = path.path_loads[current_node] + self.load_per_req[node] #calculates the load of the vehicle after visiting current_node
                
          
            if cummulative_load <= self.vehicle_capacity and time_node <= self.time_window[delivery_node_insert][1] and time_node >= self.time_window[delivery_node_insert][0]:
                cost = self.g.edge_label(current_node,node) + path.path_costs[current_node] 
                dij = self.get_reduced_cost(current_node, node, duals)
                dij_cost = dij + path.path_reduced_costs[current_node]

                          
                path_new = self.update(deepcopy(path), delivery_node_insert, dij_cost, cost, time_node,cummulative_load, delivery_index + 1, duals)
                
                new_path_act = path_new.path_ids
                try:
                    paths_omega_same_length = set_of_paths_hash_table[frozenset(new_path_act)]
                except:
                    paths_omega_same_length = []
                
                
                
                if new_path_act in paths_visited or new_path_act in paths_omega_same_length:
                    delivery_inds_available.remove(delivery_index )                    
                    continue
                
                    
                delivery_inds_available.remove(delivery_index)
                
                if path_new.d < min_cost_path and path_new.path_ids not in paths_omega_same_length:
                    final_path = path_new
                    min_cost_path = path_new.d
                    final_ind = delivery_index
            else:
                delivery_inds_available.remove(delivery_index )
                continue
       
        return final_path
        
    def calculate_time_node(self, node1, node2, path):
        """
        Calculates time from node1 to node2 and returns the time after adding on the path corresponding to label L the time from node1 to node2
        :node1: int
        :node2: int
        :L: Label object 
        """
        if  node1 != 0:
            time_node = self.travel_time(node1, node2) + path.path_times[node1]
        else:
            time_node = self.time_window[node2][0]
        return time_node
    
    def randomizedInsert_new(self, path, pickup_removed_index, del_removed_index, paths_visited, request_removed, duals,  set_of_paths_hash_table):
        # Insert request
        request_nodes, delivery_nodes = self.get_sorted_list_on_reduced_cost(path.path_ids[pickup_removed_index-1], duals, [])
        path_new  = None
        distances, incident_nodes_sorted = zip(*request_nodes)
        incident_nodes_sorted  = [inc_node for inc_node in incident_nodes_sorted if inc_node not in path.path_ids]
        insert_flag = False
       
        while insert_flag == False:
            request_nodes_chosen = self.get_random_nodes(incident_nodes_sorted, min(len(incident_nodes_sorted), 4))
       
            try:
                request_nodes_chosen.remove(request_removed)
            except:
                pass
                
            ind = pickup_removed_index-1
            for req in request_nodes_chosen:
                current_node = path.path_ids[pickup_removed_index-1]
                time_node = self.calculate_time_node(current_node, req, path) #calculates the arrival time on the current_node
                cummulative_load = path.path_loads[current_node] + self.load_per_req[req] #calculates the load of the vehicle after visiting current_node
 
                #checking capacity and time constraints
                if current_node == 0 or (cummulative_load <= self.vehicle_capacity and time_node <= self.time_window[req][1] and time_node >= self.time_window[req][0]):
                    try:
                        next_node = path.path_ids[ind+1]
                        self.g.edge_label(req, next_node)
                    except:
                        continue
                    cost = self.g.edge_label(current_node,req) + path.path_costs[current_node] 
                    dij = self.get_reduced_cost(current_node, req, duals)
                    dij_cost = dij + path.path_reduced_costs[current_node]
                    path_new = self.update(path, req, dij_cost, cost, time_node,cummulative_load, ind+1, duals)
                    insert_flag = True
                    break
            
            if path_new  == None:
                return None
                
            delivery_inds_available = range(pickup_removed_index , len(path_new.path_ids)-1)
         
            if delivery_inds_available == []:
                delivery_inds_available = [pickup_removed_index]
            
            path = self.insert_delivery(path_new, delivery_inds_available, req + self.n, paths_visited, duals,  set_of_paths_hash_table)   
        return path
   
    def update(self, path, new_node_id, dij_cost, cost, time_node,cummulative_load, ind, duals): 
 
        path.path_ids.insert(ind, new_node_id)
        path.path_reduced_costs[new_node_id] = dij_cost 
        path.path_costs[new_node_id] = cost 
        path.path_loads[new_node_id] = cummulative_load      
        path.path_times[new_node_id] = time_node
        path.V.union(Set([new_node_id]))
        path.O.union(Set([new_node_id]))
 
        for i in range(ind, len(path.path_ids)-1):
            node =  path.path_ids[i+1]
            
            try:
                previous_node = path.path_ids[i]
                path.path_reduced_costs[node] = path.path_reduced_costs[previous_node] + self.get_reduced_cost(previous_node, node, duals)
                path.path_costs[node] = path.path_costs[previous_node] + self.g.edge_label(previous_node, node)
                path.path_loads[node] = path.path_loads[previous_node] + self.load_per_req[node]
                path.path_times[node] = path.path_times[previous_node] + self.calculate_time_node(previous_node, node, path)
            except:
                continue
        return path
        
    def LNS(self, path_in, sigma, duals,  set_of_paths_hash_table):
        unsigned_oo = np.inf
        request_removed_lst = []
        paths_visited = []
        path = copy(path_in)
        available_requests_to_remove = [r for r in path.path_ids if r in self.P]
        initial_path = path.path_ids
        
        for i in range(sigma):
            f = unsigned_oo
            while path.d < f:
                f = path.d
                path_new = copy(path)
                path_new, request_removed, pickup_removed_index, del_removed_index = self.removeNodes(path_new, available_requests_to_remove, duals)
                
                if path_new is None:
                    return None
                
                
                path_new_full  = self.randomizedInsert_new(path_new,pickup_removed_index, del_removed_index , paths_visited, request_removed, duals,  set_of_paths_hash_table)
                
              
                if path_new_full is None:
                    available_requests_to_remove.remove(request_removed)
                    continue
                
                
                paths_visited.append(path_new_full.path_ids)
              
                try:
                    paths_omega_same_length = set_of_paths_hash_table[frozenset(path_new_full.path_ids)]
                except:
                    paths_omega_same_length = []
                if path_new_full.d <= f: 
                    
                    if path_new_full.d < 0 and path_new_full != initial_path and path_new_full.path_ids not in paths_omega_same_length:
                        return path_new_full
                        
                        
                    path = path_new_full
                    available_requests_to_remove = [r for r in path.path_ids if r in self.P]
                else:
                    available_requests_to_remove.remove(request_removed)
      
       
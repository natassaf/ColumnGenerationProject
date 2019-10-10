import pandas as pd
import pickle
from itertools import chain
import time
from sage.all import *
import sys
import constrained_shortest_path as ss
import read_data_script as read_data 
from LNS import LNS 
from IPython.display import display
from yappi import get_func_stats, COLUMNS_FUNCSTATS, COLUMNS_THREADSTATS
import branch_and_bound_pending as bnb

#General Functions
################################

def save_data(data, url):
    with open(url, "wb") as file:
        data_dict = pickle.dump(data,file)

def get_paths_solution(res_variables, omega):
    """
    :res_variables: RMP solver result. Dictionary of number of variable as key and resulting value as value
    :omega: list of paths in solution space of the RMP (restricted master problem)
    Return: list of paths in solution.
    """
    paths_in_solution = []
    for key, value in res_variables.items():
        if value != 0 and omega[key]  :
            paths_in_solution.append(omega[key])
    return paths_in_solution

def print_results(res_variables, omega, r_cost):
    num_of_vehicles_used = 0
    
    path_list = []
    total_cost = 0
    print("Path is:")
    index_list = []
    for key, value in res_variables.items():
        if value != 0:
            num_of_vehicles_used += 1
            path_list.append(omega[key])
            index_list.append(key)
            print(omega[key],value)
            print("objective func_per_path",r_cost[key])
            total_cost += r_cost[key]
     
    #remove from total cost the cost from node 0  to every node and the cost from any node to node 2*n + 1
    total_cost = total_cost - 10000 * num_of_vehicles_used - 1000 * num_of_vehicles_used
    print("num_of_vehicles_used", num_of_vehicles_used)
    print("total cost", total_cost)
    
    return index_list
    
def create_a(set_of_paths, vs):
    """
    Given the current set of paths and the vertices of the graph return the A matrix used in the constraints indicating how many time node i appears in path j 
    :set_of_paths: paths is RMP solution space
    :vs: list of vertices 
    Return: Matrix a. Sagemath automatically use an effective way of storing sparse matrices
    """
    a = matrix(len(vs),len(set_of_paths))
    for r in range(len(set_of_paths)):
        a[:,r] = 0
        for i in vs:
            if i in set_of_paths[r]:
                a[i,r] += 1
    return a

def save_in_dict(g, data_dict, num_of_paths_in_sol, res_variables, res_obj, omega, r_cost, n, num_set_of_paths, count_iter_same_obj, day, set_of_paths_hash_table, duals):
    res_dict = {}
    res_dict["num_paths_in_solution"] = num_of_paths_in_sol
    res_dict["res_variables"] = res_variables
    res_dict["res_obj"] = res_obj
    res_dict["omega"] = omega
    res_dict["r_cost"] = r_cost
    res_dict["n"] = n
    res_dict["num_set_of_paths"] = num_set_of_paths 
    res_dict["num_of_redundant_iterations"] = count_iter_same_obj   
    res_dict["duals"] = duals
    res_dict["time"] =  ending_time - starting_time
    g.save("C:/Users/nataz/Desktop/Github_repositories/PDPTW/code/results/graph_pdptw")
    with open('C:/Users/nataz/Desktop/Github_repositories/PDPTW/code/results/set_of_paths_hash_table'+ str(n) , 'wb') as handle:
        pickle.dump(set_of_paths_hash_table, handle)
    with open('C:/Users/nataz/Desktop/Github_repositories/PDPTW/code/results/results_'+ str(n) , 'wb') as handle:
        pickle.dump(res_dict, handle)
    with open('C:/Users/nataz/Desktop/Github_repositories/PDPTW/code/results/data_dict_'+ str(n), 'wb') as handle:
        pickle.dump(data_dict, handle)
        
def get_paths_in_sol(res_variables, omega_labels, omega, r_cost):
    '''
    Retrieves the paths in label form that are included in the current solution
    '''
    path_list = []
    num_of_paths_in_sol = 0
    for key, value in res_variables.items():
        if value != 0:
            path_list.append(omega_labels[key])
            num_of_paths_in_sol += 1
    return path_list, num_of_paths_in_sol   

################Initial Path Functions################

def insert_single_request_paths(g, data_dict, P, D):
    """
    :P: set of pickup nodes
    :D: set of delivery nodes
    return: list of single request paths, list of corresponding costs
    """
    single_route_paths = []
    single_route_costs = []
    for req in P:
        single_route_paths.append([0, req, req + data_dict["num_of_requests"], 2* data_dict["num_of_requests"] + 1])
        single_route_costs.append(g.edge_label(0,req) + g.edge_label(req, req + data_dict["num_of_requests"]) + g.edge_label(req + data_dict["num_of_requests"],2* data_dict["num_of_requests"] + 1))
    return single_route_paths, single_route_costs
    
def create_init_path(g, data_dict, pricing_problem_obj):
    """
    Creates a initial path
    """
    path = []
    c = []
    set_of_paths = []
    set_of_costs = []
    set_of_paths, set_of_costs = insert_single_request_paths(g, data_dict, data_dict["P"], data_dict["D"])
    paths_labels = pricing_problem_obj.convert_path_to_labels(set_of_paths, g)

    return set_of_paths, set_of_costs, paths_labels
    
################Column Generation Functions ################  
    
def solve_simplex(a, r_cost, vs, set_of_paths):
    """
    :p: MixedIntegerLinearProgram object
    :a: how many times row(node) i appears in column(path) j
    :r_cost: distances matrix used to define the objective function 
    :vs: list of vertices
    :set_of_paths: list of paths used as variables
    return: dictionary of variables and the values of the solution, cost of solution, list of dual variables
    """
    
    #create variable
    p =  MixedIntegerLinearProgram(maximization = False)
    p.solver_parameter("simplex_or_intopt", "simplex_only") 

    y = p.new_variable(integer = False, nonnegative=True)
    
    #create objective function
    p.set_objective(p.sum(r_cost[r] * y[r] for r in range(len(set_of_paths))))
    vs_size = len(vs)
    for i in vs[1:vs_size-1]:
        if i != 0 and i != len(vs)-1:
            p.add_constraint(p.sum(a[i,r] * y[r] for r in range(len(set_of_paths)))==1) 

    #solve the master problem
    p.solve()
    
    d = p.get_backend()
    
    duals = []
    for i in vs[:vs_size-2]:
        duals.append(round(d.get_row_dual(i),2))
    duals.insert(0, 0)
    duals.append(0)
    
    return p.get_values(y), p.get_objective_value(), duals

def run_H3(set_of_paths_hash_table, duals, obj):
    """
    In this heuristic we 
    :set_of_paths_hash_table: paths already in RMP in dictionary with frozenset of nodes as keys
    "obj: object of constraint_shortest_path script
    
    """
    #run H3 
    h3_start_time = time.time()
    h3_res, labels_ending_at  = obj.calculate_H3(set_of_paths_hash_table, duals)
    h3_end_time = time.time()
    if h3_res != []:
        #all_res.extend(h3_res)
        h3_paths, h3_costs, h3_reduced_costs = obj.get_paths(h3_res)
        print("time h3",h3_end_time - h3_start_time)
        return h3_paths, h3_costs, h3_reduced_costs , h3_res, labels_ending_at
    return [],[],[],[], labels_ending_at
    
def run_H3_all(set_of_paths_hash_table, duals, obj):
    #run H3
    h3_start_time = time.time()
    h3_res, labels_ending_at  = obj.calculate_H3(set_of_paths_hash_table, duals, True)
    h3_end_time = time.time()
    
    if h3_res != []:
        print("h3 time", h3_end_time - h3_start_time)
        h3_paths, h3_costs, h3_reduced_costs = obj.get_paths(h3_res)
        return h3_paths, h3_costs, h3_reduced_costs , h3_res, labels_ending_at
    return [],[],[],[], labels_ending_at

def run_H1(set_of_paths_hash_table, duals, labels_ending_at, num_of_arcs, obj,timeout = 5, weight_low_ranks= 1.0, weight_high_ranks = 3.0):
    #run H1
    start_h1 = time.time()
    h1_res, labels_ending_at = obj.calculate_H1(set_of_paths_hash_table, duals, num_of_arcs, labels_ending_at, weight_low_ranks, weight_high_ranks, timeout)
    end_h1 = time.time()
    
    if h1_res != []:
        print("h1 time", end_h1 - start_h1)
        h1_paths, h1_costs, h1_reduced_cost,  = obj.get_paths(h1_res)
        return h1_paths, h1_costs, h1_reduced_cost, h1_res,labels_ending_at
    return [],[],[],[], labels_ending_at

def run_sp1(set_of_paths_hash_table, duals,labels_ending_at, obj, timeout = 3600):
    start_time_sp1 = time.time()
    sp1_res = obj.calculate_SP1(set_of_paths_hash_table, duals, labels_ending_at, timeout = 3600)
    end_time_sp1 = time.time()

    if sp1_res != []:
        print("calculate_SP1",end_time_sp1 - start_time_sp1,)
        sp1_paths, sp1_costs, sp1_reduced_cost = obj.get_path(sp1_res)
        return [sp1_paths], [sp1_costs], [sp1_reduced_cost], [sp1_res]
    else:
        return [], [], [], []


def run_sp1_no_dom(set_of_paths_hash_table, duals,labels_ending_at, obj):
    start_time_sp1 = time.time()
    sp1_res = obj.calculate_SP1_no_dom(set_of_paths_hash_table, duals,labels_ending_at)
    end_time_sp1 = time.time()

    if sp1_res != []:
        print("calculate_SP1",end_time_sp1 - start_time_sp1)
        sp1_paths, sp1_costs, sp1_reduced_cost = obj.get_path(sp1_res)
        return [sp1_paths], [sp1_costs], [sp1_reduced_cost], [sp1_res]
    else:
        return [], [], [], []

        
def run_LNS(paths_solution_lab, sigma, duals,  set_of_paths_hash_table, obj):
    ls_obj = LNS(obj.g, obj.data_dict)
    start_time_lns = time.time()
    for path_lab in paths_solution_lab:
        new_path = ls_obj.LNS(path_lab, sigma, duals,  set_of_paths_hash_table)
    end_time_lns = time.time()
    print("LNS time:", end_time_lns - start_time_lns)
    return new_path

      
def run_all_heuristics(obj, n, duals, a, paths_solution_lab, set_of_paths_hash_table):
    '''
    Solves heuristics (H3, H4, H1) and exact algorithm sp1 zith the objective to find q negqtive reduced cost path
    '''
    #initialize variables
    all_paths = []
    all_costs = []
    all_reduced_costs = []
    all_res = []
    h1_res = []
    paths_found_labels = []
    iter = 0
    obj.reduced_cost_dict = {}
    obj.transformed_reduced_cost = {}
    #H3 Heuristic
    print("H3 heuristic is running one request")
    h3_paths, h3_costs, h3_reduced_costs, h3_res, labels_ending_at = run_H3(set_of_paths_hash_table, duals,obj)
    if h3_paths != []:
        print("Path found!")
        return h3_paths, h3_costs, h3_reduced_costs , h3_res
    print("H3 heuristic is running all request")
    
    h3_paths, h3_costs, h3_reduced_costs, h3_res, labels_ending_at = run_H3_all(set_of_paths_hash_table, duals, obj)
    if h3_paths != []:
        print("Path found!")
        return h3_paths, h3_costs, h3_reduced_costs , h3_res
    
    sigma = 5
    print("LNS is running for paths in solution and sigma =",sigma)
    lns_res = run_LNS(paths_solution_lab, sigma, duals,  set_of_paths_hash_table, obj)
   
    if lns_res !=None:
        lns_paths, lns_costs, lns_reduced_cost = obj.get_path(lns_res)
        print("Path found!")
        return [lns_paths], [lns_costs], [lns_reduced_cost], [lns_res]
    

    num_of_arcs = n//4 + 2
    print("H1 heuristic is running small")
    h1_paths, h1_costs, h1_reduced_cost, h1_res, labels_ending_at = run_H1(set_of_paths_hash_table, duals, labels_ending_at, num_of_arcs, obj)
    if h1_paths != []:
        print("Path found!")
        return h1_paths, h1_costs, h1_reduced_cost, h1_res
  
    num_of_arcs = n//4 + 3
    print("H1 heuristic is running low rank")
    weight_low = 2.0
    weight_high = 1.0
    h1_paths, h1_costs, h1_reduced_cost, h1_res, labels_ending_at = run_H1(set_of_paths_hash_table, duals, labels_ending_at, num_of_arcs, obj, weight_low, weight_high )
    if h1_paths != []:
        return h1_paths, h1_costs, h1_reduced_cost, h1_res
       
    print("H1 heuristic is running large")
    num_of_arcs = num_of_arcs + n//4 
    timeout = 30
    h1_paths, h1_costs, h1_reduced_cost, h1_res, labels_ending_at = run_H1(set_of_paths_hash_table, duals, labels_ending_at, num_of_arcs, obj,timeout)
    if h1_paths != []:
        print("Path found!")
        return h1_paths, h1_costs, h1_reduced_cost, h1_res
    #labels_ending_at = {}
    
    print("SP1 is running")
    sp1_paths, sp1_costs, sp1_reduced_cost, sp1_res = run_sp1(set_of_paths_hash_table, duals,labels_ending_at,obj, timeout = 3600)
    if sp1_paths != []:
        print("Negative reduced cost found!")
        return sp1_paths, sp1_costs, sp1_reduced_cost, sp1_res
    
    '''
    labels_ending_at = {}   
    print("SP1 is running_no dom")
    sp1_paths, sp1_costs, sp1_reduced_cost, sp1_res = run_sp1_no_dom(set_of_paths_hash_table, duals,labels_ending_at,obj)
    '''
    return sp1_paths, sp1_costs, sp1_reduced_cost, sp1_res
    
def column_generation_solver(g, data_dict, a, r_cost, set_of_paths, pricing_problem_obj, set_of_paths_hash_table, omega_in_labels):
    """
    Finds a integer or fractional solution
    """
    start_of_master_problem = time.time()
    #sigma = (pricing_problem_obj.n) // 2 
    vs = g.vertices()
   
    pricing_result_paths = set_of_paths #initial set  of paths to get a 1st solution
    n = data_dict["num_of_requests"]
    
    #solve the relaxed master problem using simplex (fractional solutions are allowed)
    res_variables, res_obj, duals = solve_simplex(a, r_cost, vs, set_of_paths)
    
    
    res_obj_old = 0
    count_iter_same_obj = 0
    paths_solution_lab, num_of_paths_in_sol = get_paths_in_sol(res_variables, omega_in_labels, set_of_paths, r_cost)

    #keep solving the pricing problem as long as we  find negative reduced cost paths 
    while pricing_result_paths != []: 
        pricing_result_paths, pricing_result_cost, pricing_result_reduced_cost = [], [], []
        
        paths_in_solution = get_paths_solution(res_variables, set_of_paths) #used in local search
        
        #run the heuristics in the following order (H3, H4, H1, SP1)
        start_heuristics = time.time()
        pricing_result_paths, pricing_result_cost, pricing_result_reduced_cost, all_res  = run_all_heuristics(pricing_problem_obj ,n,  duals, a, paths_solution_lab, set_of_paths_hash_table) 
        end_heuristics = time.time()
        
        #add the path found in the set of paths, calculate matrix a and resolve the master problem
        if pricing_result_paths != []:
            
            current_of_master_problem = time.time()
            count_iter_same_obj += 1
            omega_in_labels.extend(all_res)
            set_of_paths.extend(pricing_result_paths)
            #introduce the new path in the hash table            
            for p in pricing_result_paths:
                try:
                    if p not in set_of_paths_hash_table[frozenset(p)]:
                        set_of_paths_hash_table[frozenset(p)].append(p)
                except:
                    set_of_paths_hash_table[frozenset(p)]= [p]
            r_cost.extend(pricing_result_cost)
            a = create_a(set_of_paths, vs)
        res_variables, res_obj, duals = solve_simplex(a, r_cost, vs, set_of_paths)
      
        res_obj_old =  res_obj               
        paths_solution_lab, num_of_paths_in_sol = get_paths_in_sol(res_variables, omega_in_labels, set_of_paths, r_cost)
        print("Number of paths in the master problem", len(set_of_paths))

    return res_variables, res_obj, set_of_paths, r_cost, count_iter_same_obj, num_of_paths_in_sol, set_of_paths_hash_table, duals
    
def column_generation_main(g, data_dict, P, D, day, n):
    """
    Create an initial path and call the column_generation_solver function
    :g: A sage graph 
    :data_dict: dictionary containing the data mentioned in the report
    :P: pickup nodes
    :D: delivery nodes
    :day: 
    :n: number of requests
    :return: fractional or initial solution - > res_variables: dictionary of variables and values, res_obj:total cost of solution, 
           omega: paths in master problem, r_cost: cost of each path in the rmp, count_iter_same_obj: number of iterations with the same obj function, num_of_paths_in_sol:number of paths used in solution which is num of vehicles, set_of_paths_hash_table: discovered paths in dictionary form to accelerate the pace.  
    """
    path_list = []
    set_of_paths = []
    set_of_paths_hash_table = {}
    global starting_time
    global ending_time
    
    #create object of the pricing problem solvers class and local search class
    pricing_problem_obj = ss.PricingSubproblems(g, data_dict)
    
    #create as initial paths the paths containing one request and the corresponding delivery
    initial_paths, r_cost, omega_labels = create_init_path(g, data_dict, pricing_problem_obj)
    
    
    set_of_paths.extend(initial_paths)
    set_of_paths_hash_table[frozenset([0,2*n+1])]= [0,2*n+1]
    #create hash table that will accelerate the check on whether a path exists already in the master problem columns.
    for p in initial_paths:
        try:
            set_of_paths_hash_table[frozenset(p)].append(p)
        except:
            set_of_paths_hash_table[frozenset(p)] = [p]

    #create sparse matrix to store matrix a of the paper. It indicates how many times a variable is contained in a given column(path)
    a = matrix(RR,g.order(),len(set_of_paths),sparse = True)
    a = create_a(set_of_paths, g.vertices())
    
    #solve the problem using simplex and column generation
    res_variables, res_obj, omega, r_cost, count_iter_same_obj, num_of_paths_in_sol, set_of_paths_hash_table, duals = column_generation_solver(g, data_dict, a, r_cost, set_of_paths, pricing_problem_obj, set_of_paths_hash_table,omega_labels)
    
    
    #save results
    ending_time = time.time()
    save_in_dict(g, data_dict, num_of_paths_in_sol, res_variables, res_obj, omega, r_cost, n, len(omega), count_iter_same_obj, day, set_of_paths_hash_table, duals)
    return res_variables, res_obj, omega, r_cost, count_iter_same_obj, num_of_paths_in_sol, set_of_paths_hash_table, pricing_problem_obj
    
#####Branch and Bound Functions#####
def solve_simplex_branch(a, r_cost, vs, set_of_paths, u, bound_type):
    """
    :p: MixedIntegerLinearProgram object
    :a: how many times row(node) i appears in column(path) j
    :r_cost: distances matrix used to define the objective function 
    """
    
    #create variable
    p =  MixedIntegerLinearProgram(maximization = False)
    p.solver_parameter("simplex_or_intopt", "simplex_only") 
    #p.solver_parameter("primal_v_dual", "GLP_DUAL")
    y = p.new_variable(integer = False, nonnegative=True)
    #create objective function
    
    p.set_objective(p.sum(r_cost[r] * y[r] for r in range(len(set_of_paths))))
    vs_size = len(vs)
    num_of_added_consraints = 0
    for i in vs[1:vs_size-1]:
        if i != 0 and i != len(vs)-1:
            p.add_constraint(p.sum(a[i,r] * y[r] for r in range(len(set_of_paths)))==1) 
    if bound_type == "ceil":
        p.add_constraint(p.sum(y[r] for r in range(len(set_of_paths))) >= ceil(u))
        num_of_added_consraints += 1
    elif bound_type == "floor":
        p.add_constraint(p.sum(y[r] for r in range(len(set_of_paths))) <= floor(u))
        num_of_added_consraints += 1
        
    #solve the master problem
    p.solve()
    d = p.get_backend()
    
    duals = []
    for i in vs[:vs_size-2]:
        duals.append(round(d.get_row_dual(i),1))
    for j in range(1,num_of_added_consraints+1):
        duals.append(round(d.get_row_dual(i + j),1))
    #print(d.is_variable_binary(0))
    
    duals.insert(0, 0)
    duals.append(0)
    return p.get_values(y), p.get_objective_value(), duals

def correct_rounding_error_sol_paths(res_variables):
    u = 0
    e = 0.0001
    for key, value in res_variables.items():
        if value != 0:
            if value <= e:
                res_variables[key] = 0
            else:
                if abs(res_variables[key]  - int(res_variables[key]))<0.01:
                    res_variables[key] = round(res_variables[key])
                else:
                    res_variables[key] = round(res_variables[key], 2)
                u += value
                

    return u, res_variables     

    
def create_month_to_string_dict():
    months_dict = {}
    months = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "July", "Aug", "Sep", "Oct", "Nov", "Dec"]
    for i in range(12):
        months_dict[i] = months[i]
    return months_dict
   
def check_branch_and_bound(res_variables, omega):
    path_list = []
    values_of_solution = []
    branch_and_bound_flag = False #no branch and bound needed
    for key, value in res_variables.items():
        if value != 0:
            path_list.append(omega[key])
            values_of_solution.append(value)
            if value != 1:
                return True
    return False
 
def print_actual_path(res_variables, omega, mapping_dict_num_to_code,n):
    path_list = []
    
    print("The actual stops to follow are:")
    index_list = []
    for key, value in res_variables.items():
        node_temp = []
        if value != 0:
            for node in omega[key]:
                if node != 0 and node != 2*n+1:
                    
                    node_temp.append(mapping_dict_num_to_code[node])
            path_list.append(node_temp)
    print(path_list)
 
################ Main Function ################
if __name__ == "__main__":
    global starting_time 
    global ending_time
        
    requests_data = read_data.readcsv("C:/Users/nataz/Desktop/Github_repositories/PDPTW/requests")#sys.argv[1]
    month_dict = create_month_to_string_dict()

    #Get input from user
    user_input = raw_input("Please insert the date of requests you wish to run the algorithm at in the format day/month/year: (ex: 21/12/2017) ")
    n_selection = raw_input("Please select he number of requests you want to optimize: (ex:12)")
    input = user_input.split("/")
    input[1] = month_dict[int(input[1])-1]

    daily_data = read_data.choose_data(requests_data, int(input[0]) ,input[1],int(input[2]))
    if daily_data.empty:
        print("No data found for this date")
        sys.exit()

    #time is counted in seconds
    T = 60 * 8 * 60 #planning horizon  
    W = 30*60 #time window interval in seconds
    n = int(n_selection)

    
    #create data needed 
    g, P, D, time_windows, mapping_dict_num_to_code = read_data.create_graph(daily_data[:n])
    data_dict = read_data.create_data_dict(g, P, D, time_windows, n, T, W)
    g = read_data.graph_preprocessing(g, n, P, D, data_dict)


    #Run Column Generation
    starting_time = time.time()
    res_variables, res_obj, omega, r_cost, count_iter_same_obj, num_of_paths_in_sol, set_of_paths_hash_table, pricing_problem_obj  = column_generation_main(g, data_dict, P, D, input[0], n)
    end = time.time()
    print("Time to run Column Generation: ",end - starting_time)

    #Discard the paths that are the result of rounding errors while solving simplex
    #u, res_variables = correct_rounding_error_sol_paths(res_variables)

    #Print Column Generation Result
    print_results(res_variables, omega, r_cost)

    #Check if branch and bound is needed
    bnb_flag = check_branch_and_bound(res_variables, omega)
    if bnb_flag:
        print("Starting Branch and Bound")
        bnb.run_branch_and_bound(g, data_dict, res_variables, res_obj,  r_cost, omega, set_of_paths_hash_table, pricing_problem_obj) 
    else:
        print_actual_path(res_variables, omega, mapping_dict_num_to_code, n)
        


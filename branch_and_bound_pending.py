import pickle
import yappi
from sage.all import *
#import master_problem_final as master
import time 
import constrained_shortest_path_changed as ss
from collections import defaultdict

class SolutionInfo:
    def __init__(self, res_variables,objective_func_cost, a, r_cost, set_of_paths, set_of_paths_hash_table, num_of_requests, branch_type = "branching_rule_1"):
        self.matrix_a = a
        self.solver_result_dict = res_variables
        self.r_cost = r_cost
        self.set_of_paths = set_of_paths
        self.set_of_paths_hash_table = set_of_paths_hash_table
        self.n = num_of_requests
        self.u = self.calculate_u()
        self.objective_func_cost = objective_func_cost
        self.next_branching_method = branch_type
        
    def calculate_u(self):
        #calculates the u variable of branch and bound level 1 rule in Desrochers et al(1992) paper/section 5  
        u = 0
        for key, value in self.solver_result_dict.items():
            if value != 0:
                u += value
        return u
 
 
 
def deepcopy(sol_info_object):
    
    matrix_a = copy(sol_info_object.matrix_a)
    solver_result_dict = copy(sol_info_object.solver_result_dict)
    r_cost = copy(sol_info_object.r_cost)
    set_of_paths = copy(sol_info_object.set_of_paths)
    set_of_paths_hash_table = copy(sol_info_object.set_of_paths_hash_table)
    n = copy(sol_info_object.n)
    u = copy(sol_info_object.u)
    objective_func_cost = copy(sol_info_object.objective_func_cost)
    next_branching_method = copy(sol_info_object.next_branching_method)
    return SolutionInfo(solver_result_dict, objective_func_cost, matrix_a, r_cost, set_of_paths, set_of_paths_hash_table, n)
    
def create_a(set_of_paths, n):
    """
    Given the current set of paths and the vertices of the graph return the A matrix used in the constraints indicating how many time node i appears in path j 
    """
    vs = range(2* n +2)
    a = matrix(len(vs),len(set_of_paths))
    for r in range(len(set_of_paths)):
        a[:,r] = 0
        for i in vs:
            if i in set_of_paths[r]:
                a[i,r] += 1
    return a
    
def calculate_u(solver_result_dict):
    #calculates the u variable of branch and bound level 1 rule in Desrochers et al(1992) paper/section 5  
    u = 0
    for key, value in solver_result_dict.items():
        if value != 0:
            u += value
    return u
   

def solve_simplex_branch_lev2(sol_info, bound_type):
    """
    :p: MixedIntegerLinearProgram object
    :a: how many times row(node) i appears in column(path) j
    :r_cost: distances matrix used to define the objective function 
    """
  
    vs = range(2* sol_info.n +2)
    #create variable
    p =  MixedIntegerLinearProgram(maximization = False)
    p.solver_parameter("simplex_or_intopt", "simplex_only") 
    #p.solver_parameter("primal_v_dual", "GLP_DUAL")
    y = p.new_variable(integer = False, nonnegative=True)
    #create objective function
    
    p.set_objective(p.sum(sol_info.r_cost[r] * y[r] for r in range(len(sol_info.set_of_paths))))
    vs_size = len(vs)
    num_of_added_consraints = 0
    for i in vs[1:vs_size-1]:
        if i != 0 and i != len(vs)-1:
            p.add_constraint(p.sum(sol_info.matrix_a[i,r] * y[r] for r in range(len(sol_info.set_of_paths)))==1) 
    if bound_type == "ceil":
        print("lev2 bound ceil",ceil(sol_info.objective_func_cost))
        p.add_constraint(p.sum(y[r]* sol_info.r_cost[r] for r in range(len(sol_info.set_of_paths))) > ceil(sol_info.objective_func_cost))
        p.add_constraint(p.sum(y[r] for r in range(len(sol_info.set_of_paths))) >= ceil(sol_info.u))
        num_of_added_consraints += 2
    elif bound_type == "floor":
        print("lev2 bound floor",ceil(sol_info.objective_func_cost))
        p.add_constraint(p.sum(y[r] for r in range(len(sol_info.set_of_paths))) >= ceil(sol_info.u))
        p.add_constraint(p.sum(y[r] * sol_info.r_cost[r] for r in range(len(sol_info.set_of_paths))) < floor(sol_info.objective_func_cost))
        num_of_added_consraints += 2
        
    #solve the master problem
    
    p.solve()
    
    d = p.get_backend()
    
    duals = []
    for i in vs[:vs_size-2]:
        #print("i",i, "dual", d.get_row_dual(i))
        duals.append(d.get_row_dual(i))
    for j in range(1,num_of_added_consraints+1):
        duals.append(d.get_row_dual(i + j))
    #print(d.is_variable_binary(0))
    
    duals.insert(0, 0)
    duals.append(0)
    print("duals_lev2",duals)
    #p.show()
    #print(set_of_paths)

    return p.get_values(y), p.get_objective_value(), duals
    
def solve_simplex_branch_lev1(sol_info, bound_type):
    """
    :p: MixedIntegerLinearProgram object
    :a: how many times row(node) i appears in column(path) j
    :r_cost: distances matrix used to define the objective function 
    """
  
    vs = range(2* sol_info.n +2)
    #create variable
    p =  MixedIntegerLinearProgram(maximization = False)
    p.solver_parameter("simplex_or_intopt", "simplex_only") 
    #p.solver_parameter("primal_v_dual", "GLP_DUAL")
    y = p.new_variable(integer = False, nonnegative=True)
    #create objective function
    
    p.set_objective(p.sum(sol_info.r_cost[r] * y[r] for r in range(len(sol_info.set_of_paths))))
    vs_size = len(vs)
    num_of_added_consraints = 0
    for i in vs[1:vs_size-1]:
        if i != 0 and i != len(vs)-1:
            p.add_constraint(p.sum(sol_info.matrix_a[i,r] * y[r] for r in range(len(sol_info.set_of_paths)))==1) 
    if bound_type == "ceil":
        print("current bound u",ceil(sol_info.u))
        p.add_constraint(p.sum(y[r] for r in range(len(sol_info.set_of_paths))) >= ceil(sol_info.u))
        num_of_added_consraints += 1
    elif bound_type == "floor":
        print("current bound u",floor(sol_info.u))
        p.add_constraint(p.sum(y[r] for r in range(len(sol_info.set_of_paths))) <= floor(sol_info.u))
        num_of_added_consraints += 1
        
    #solve the master problem
    
    p.solve()
    
    d = p.get_backend()
    
    duals = []
    for i in vs[:vs_size-2]:
        #print("i",i, "dual", d.get_row_dual(i))
        duals.append(d.get_row_dual(i))
    for j in range(1,num_of_added_consraints+1):
        duals.append(d.get_row_dual(i + j))
    #print(d.is_variable_binary(0))
    
    duals.insert(0, 0)
    duals.append(0)
    print("duals_lev1",duals)
    #p.show()
    #print(set_of_paths)

    return p.get_values(y), p.get_objective_value(), duals
        
        
#na ftiaksw klash gia na exw opote thelw ta g, datadict P, D    
 #DEN DOULEUEI AYTO
 

def generate_columns(sol_info, subproblem_obj, duals):
    start_heuristics = time.time()
    #pricing_problem_obj = ss.PricingSubproblems(g, data_dict, P, D)
    
    pricing_result_paths, pricing_result_cost, pricing_result_reduced_cost, all_res  = master.run_all_heuristics(subproblem_obj ,n, [], duals, sol_info.matrix_a, [],subproblem_obj, sol_info.set_of_paths_hash_table) 
    end_heuristics = time.time()
    
    #print("pricing_result_paths, pricing_result_reduced_cost",pricing_result_paths, pricing_result_reduced_cost)
    #print(pricing_result_paths in set_of_paths)
    
    #add the path found in the set of paths, calculate matrix a and resolve the master problem
    if pricing_result_paths != []:
        current_of_master_problem = time.time()
        count_iter_same_obj += 1
        omega_in_labels.extend(all_res)
        set_of_paths.extend(pricing_result_paths)

        #print("set_of_paths",len(set_of_paths))
        
        #introduce the new path in the hash table
        for p in pricing_result_paths:
            try:
                if p not in set_of_paths_hash_table[frozenset(p)]:
                    set_of_paths_hash_table[frozenset(p)].append(p)
            except:
                set_of_paths_hash_table[frozenset(p)]= [p]
        r_cost.extend(pricing_result_cost)
        a = create_a(set_of_paths, vs)
        return set_of_paths, r_cost, a, set_of_paths_hash_table, True
    else:
        return set_of_paths, r_cost, a, set_of_paths_hash_table, False
        
def choose_branching_rule(sol_info):
    #First branching strategy is to add a constraint on the number of vehicles used. 
    #Branch_1: sum of vehicles <= floor(u), 
    #Branch 2:  sum of vehicles <= floor(u) 
    #where u is the sum of vehicles in the current solution.
    u =calculate_u(sol_info.solver_result_dict)
    if u != int(u):
        return "branching_rule_1"
    elif sol_info.objective_func_cost != int(sol_info.objective_func_cost):
        return "branching_rule_2"
    else:
        return "branching_rule_3"
    
def check_sol_integrality(solver_result_dict):
    for key, value in solver_result_dict.items():
        if value != 0 and value != 1:
            return False
    return True

def run_branch_and_bound(g, data_dict, res_variables, objective_func_cost,  r_cost, set_of_paths, set_of_paths_hash_table, subproblem_obj):
    #Store the fractional solution info
    num_of_requests = data_dict["num_of_requests"]
    a = create_a(set_of_paths, num_of_requests)
    result = []
    #create first sol_info object with branching strategy 1
    fractional_sol_info = SolutionInfo(res_variables, objective_func_cost, a, r_cost, set_of_paths, set_of_paths_hash_table, num_of_requests, "branching_rule_1")
    u, res_variables = correct_rounding_error_sol_paths(res_variables)
    #master.print_results(res_variables, fractional_sol_info.set_of_paths, fractional_sol_info.r_cost)
    bnb_queue = []
    bnb_queue.append(fractional_sol_info)
    while bnb_queue:
        sol_info = bnb_queue.pop()
        #branch_type = choose_branching_rule(sol_info)
        print(sol_info.next_branching_method)
        
        if sol_info.next_branching_method == "branching_rule_1":
            bound_type_1 = "ceil"
            try:
                #solving master problem in branch 1
                res_variables, res_obj, duals_1 = solve_simplex_branch_lev1(sol_info, bound_type_1)
                u, res_variables = correct_rounding_error_sol_paths(res_variables)
                master.print_results(res_variables, sol_info.set_of_paths, sol_info.r_cost)
                #save solution
                if check_sol_integrality(res_variables):
                    result.append({"res_variables":res_variables,"res_obj": res_obj, "set_of_paths": sol_info.set_of_paths,"r_cost": sol_info.r_cost })
                    lowest_int_bound = res_obj
                        
                #generate new columns
                    set_of_paths, r_cost, matrix_a, set_of_paths_hash_table, found_new_paths = generate_columns(sol_info, subproblem_obj,duals_1)
                    #if no new columns can be added then append in bnb queue and check branching rule 2
                    if not found_new_paths:
                        sol_info_new = deepcopy(sol_info)
                        sol_info_new.next_branching_method = "branching_rule_2"
                        bnb_queue.append(sol_info_new)
                    else:
                        sol_info_new = SolutionInfo(res_variables, res_obj, matrix_a , r_cost, set_of_paths, set_of_paths_hash_table, num_of_requests, "branching_rule_1" )
                        bnb_queue.append(sol_info_new)
                else: 
                    sol_info_new = SolutionInfo(res_variables, res_obj, matrix_a , r_cost, set_of_paths, set_of_paths_hash_table, num_of_requests, "branching_rule_2" )
                    bnb_queue.append(sol_info_new)
            except:
                print("No feasible solution found")
                sol_info_new = deepcopy(sol_info)
                sol_info_new.next_branching_method = "branching_rule_2"
                bnb_queue.append(sol_info_new)
                
                
            bound_type_2 = "floor"
            
            try:
                print("branching_here",sol_info.next_branching_method)
                res_variables, res_obj, duals_2 = solve_simplex_branch_lev1(sol_info, bound_type_2)
                u, res_variables = correct_rounding_error_sol_paths(res_variables)
                master.print_results(res_variables, sol_info.set_of_paths, sol_info.r_cost)
                if check_sol_integrality(res_variables):
                    
                    #if res_obj < lowest_int_bound:
                    result.append({"res_variables":res_variables,"res_obj": res_obj, "set_of_paths": sol_info.set_of_paths,"r_cost": sol_info.r_cost })
                    lowest_int_bound = res_obj
                        
                    set_of_paths, r_cost, matrix_a, set_of_paths_hash_table, found_new_paths = generate_columns(sol_info, subproblem_obj,duals_2)
                    if not found_new_paths :
                        sol_info_new = deepcopy(sol_info)
                        sol_info_new.next_branching_method = "branching_rule_2"
                        bnb_queue.append(sol_info_new)
                    else:
                        sol_info_new = SolutionInfo(res_variables, res_obj, matrix_a , r_cost, set_of_paths, set_of_paths_hash_table, num_of_requests, "branching_rule_1" )
                        bnb_queue.append(sol_info_new)
                else:
                    sol_info_new = SolutionInfo(res_variables, res_obj, matrix_a , r_cost, set_of_paths, set_of_paths_hash_table, num_of_requests, "branching_rule_2" )
                    bnb_queue.append(sol_info_new)
            except:
                print("No feasible solution found")
                sol_info_new = deepcopy(sol_info)
                sol_info_new.next_branching_method = "branching_rule_2"
                bnb_queue.append(sol_info_new)
                
        elif sol_info.next_branching_method == "branching_rule_2":
            print("branching_rule_2")
            if sol_info.next_branching_method == "branching_rule_2":
                bound_type_1 = "ceil"
                try:
                    #solving master problem in branch 1
                    res_variables, res_obj, duals_1 = solve_simplex_branch_lev2(sol_info, bound_type_1)
                    u, res_variables = correct_rounding_error_sol_paths(res_variables)
                    master.print_results(res_variables, sol_info.set_of_paths, sol_info.r_cost)
                    #save solution
                    if check_sol_integrality(res_variables):
                        result.append({"res_variables":res_variables,"res_obj": res_obj, "set_of_paths": sol_info.set_of_paths,"r_cost": sol_info.r_cost })
                        lowest_int_bound = res_obj
                            
                    #generate new columns
                    set_of_paths, r_cost, matrix_a, set_of_paths_hash_table, found_new_paths = generate_columns(sol_info, subproblem_obj,duals_1)
                    #if no new columns can be added then append in bnb queue and check branching rule 2
                    if not found_new_paths:
                        sol_info_new = deepcopy(sol_info)
                        sol_info_new.next_branching_method = "branching_rule_3"
                        bnb_queue.append(sol_info_new)
                    else:
                        sol_info_new = SolutionInfo(res_variables, res_obj, matrix_a , r_cost, set_of_paths, set_of_paths_hash_table, num_of_requests, "branching_rule_2" )
                        bnb_queue.append(sol_info_new)
                except:
                    print("No feasible solution found")
                    sol_info_new = deepcopy(sol_info)
                    sol_info_new.next_branching_method = "branching_rule_3"
                    bnb_queue.append(sol_info_new)
                    
                    
                bound_type_2 = "floor"
                print("branching_rule_2")
                try:
                    print("branching_here",sol_info.next_branching_method)
                    res_variables, res_obj, duals_2 = solve_simplex_branch_lev2(sol_info, bound_type_2)
                    u, res_variables = correct_rounding_error_sol_paths(res_variables)
                    print(bound_type, "solution\n")
                    master.print_results(res_variables, sol_info.set_of_paths, sol_info.r_cost)
                    if check_sol_integrality(res_variables):
                        result.append({"res_variables":res_variables,"res_obj": res_obj, "set_of_paths": sol_info.set_of_paths,"r_cost": sol_info.r_cost })
                        #if res_obj < lowest_int_bound:
                        lowest_int_bound = res_obj
                    set_of_paths, r_cost, matrix_a, set_of_paths_hash_table, found_new_paths = generate_columns(sol_info, subproblem_obj,duals_2)
                    if not found_new_paths:
                        sol_info_new = deepcopy(sol_info)
                        sol_info_new.next_branching_method = "branching_rule_3"
                        bnb_queue.append(sol_info_new)
                    else:
                        sol_info_new = SolutionInfo(res_variables, res_obj, matrix_a , r_cost, set_of_paths, set_of_paths_hash_table, num_of_requests, "branching_rule_2" )
                        bnb_queue.append(sol_info_new)
                except:
                    print("No feasible solution found")
                    sol_info_new = deepcopy(sol_info)
                    sol_info_new.next_branching_method = "branching_rule_3"
                    bnb_queue.append(sol_info_new)
                    
        elif sol_info.next_branching_method == "branching_rule_3":
            print("Need to run branching rule 3 which is not yet implemented")
            return None
 
    
def bound2_heuristic(g, set_of_paths, res_variables, P, D):
    #create a dictionary with the arc cost using as cost the value of the solution path in the master problem
    branch_cost_dict = {}
    arcs_cost =  defaultdict(lambda:0)
    for key, value in res_variables.items():
        if value != 0:
            for i in range(1, len(set_of_paths[key])):
                arc = (set_of_paths[key][i-1], set_of_paths[key][i])  
                arcs_cost[arc] = value
                
                try:
                    branch_cost_dict[set_of_paths[key][i-1]] += round(value,2)
                except:
                    branch_cost_dict[set_of_paths[key][i-1]] = 0
                
                
    #print("branch_cost_dict",branch_cost_dict)
  
    
    result = []
    branch_set = [0]
 
    result.append((branch_set, 0))
    U = copy(result)
    
    while U:
        branch_set = U.pop()

        if len(branch_set[0]) > 4:
            return result
        nodes_in_set = branch_set[0]
        
        prev_node_cost = branch_set[1]
        
        for current_node in g.vertices():
            branch_set_new = copy(branch_set[0])
            if current_node in branch_set_new:
                continue
            
            branch_set_new.append(current_node)
            try:
                new_set_cost = prev_node_cost + branch_cost_dict[(current_node)] - sum(arcs_cost[(current_node, node)] for node in nodes_in_set)
                if len(branch_set_new)>2:
                    result.append((branch_set_new, new_set_cost) )
                U.append((branch_set_new, new_set_cost) )
            except:
                if len(branch_set_new)>2:
                    result.append((branch_set_new, prev_node_cost -  sum(arcs_cost[(current_node, node)] for node in nodes_in_set)))
                U.append((branch_set_new, new_set_cost) )
        
            
    return result
    
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
             
def read_results():
    with open('C:/Users/nataz/Desktop/Github_repositories/PDPTW/code/results/7/results_7','r') as file:
        res_dict =  pickle.load(file)
    with open('C:/Users/nataz/Desktop/Github_repositories/PDPTW/code/results/7/data_dict_7','r') as file:
        data_dict =  pickle.load(file) 
    with open('C:/Users/nataz/Desktop/Github_repositories/PDPTW/code/results/7/set_of_paths_hash_table7','r') as file:
        set_of_paths_hash_table =  pickle.load(file) 
    g = load("C:/Users/nataz/Desktop/Github_repositories/PDPTW/code/results/7/graph_pdptw")

    print(res_dict.keys())
    omega = res_dict["omega"]
    r_cost =res_dict["r_cost"]
    res_variables = res_dict["res_variables"]
    objective_func_cost = res_dict["res_obj"] 
    #set_of_paths_hash_table = res_dict["set_of_paths_hash_table"]
    print("n", res_dict["n"]) 
    print("objective function", res_dict["res_obj"])
    print("res_variables",res_dict["res_variables"])
    print("duals",res_dict["duals"])
    result = bound2_heuristic(g,omega, res_variables, data_dict["P"], data_dict["D"])
    print(type(result))
    print(result)
    best_set = []
    max_distance_from_integer= 0
    for res in result:
        if abs(round(res[1])-res[1])>= max_distance_from_integer:
            max_distance_from_integer = abs(round(res[1])-res[1])
            best_set.append(res[0])
    print(best_set)
    #for key, value in result.items():
    #    print(len(result[key]))
    branching_set = best_set[-1]
    #branch using the one of the sets obtained
    subproblem_obj = ss.PricingSubproblems(g, data_dict)
    u, res_variables = correct_rounding_error_sol_paths(res_variables)
    run_branch_and_bound(g, data_dict, res_variables, objective_func_cost,  r_cost ,omega,  set_of_paths_hash_table, subproblem_obj)


 

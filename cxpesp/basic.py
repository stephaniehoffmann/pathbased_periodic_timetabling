# 09.02.2023
# Stephanie Riedmüller

from pyscipopt import Model, quicksum, SCIP_PARAMSETTING            # type: ignore
import networkx as nx
from collections import defaultdict



def expanded_graph(L, T, weights):

    lb_driving = {}
    ub_driving = {}
    if 'duration' in L[0].keys():
        for line in L.keys():
            for i, (v,w) in enumerate(L[line]['edges']):
                lb_driving[line, (v,w)] = L[line]['duration'][i]
                ub_driving[line, (v,w)] = L[line]['duration'][i]
    else:
        for line in L.keys():
            for (v,w) in L[line]['edges']:
                lb_driving[line, (v,w)] = 1
                ub_driving[line, (v,w)] = 1
    lb_waiting = 0
    ub_waiting = 2
    lb_transfer = 1
    ub_transfer = T
    lb_turnaround = 2
    ub_turnaround = T + 1
    C = defaultdict(nx.DiGraph)
    A = defaultdict(nx.DiGraph)
    for line in L.keys():
    # driving arcs
        for (v,w) in L[line]['edges']:
            for t1 in range(T):
                for t2 in range(T):
                    duration = ((t2 - t1 - lb_driving[line, (v,w)]) % T) + lb_driving[line, (v,w)]
                    if duration >= lb_driving[line, (v,w)] and duration <= ub_driving[line, (v,w)]:
                        C[line].add_edge((v, line, 'dep', t1, '+'), (w, line, 'arr', t2, '+'), 
                                        type='driving', duration=duration, weight=1)
                        C[line].add_edge((w, line, 'dep', t1, '-'), (v, line, 'arr', t2, '-'), 
                                        type='driving', duration=duration, weight=1)
    # waiting and turnaround arcs
        for v in L[line]['nodes']:
            for t1 in range(T):
                for t2 in range(T):
                    if (v, line, 'arr', t1, '+') in C[line].nodes():
                        if (v, line, 'dep', t2, '+') in C[line].nodes():
                            duration = ((t2 - t1 - lb_waiting) % T) + lb_waiting
                            if duration >= lb_waiting and duration <= ub_waiting:
                                C[line].add_edge((v, line, 'arr', t1, '+'), (v, line, 'dep', t2, '+'), 
                                                type='waiting', duration=duration, weight=1)
                                C[line].add_edge((v, line, 'arr', t1, '-'), (v, line, 'dep', t2, '-'), 
                                                type='waiting', duration=duration, weight=1)
                        else:
                            duration = ((t2 - t1 - lb_turnaround) % T) + lb_turnaround
                            if duration >= lb_turnaround and duration <= ub_turnaround:
                                C[line].add_edge((v, line, 'arr', t1, '+'), (v, line, 'dep', t2, '-'), 
                                                type='turnaround', duration=duration, weight=1)
                    else:
                        duration = ((t2 - t1 - lb_turnaround) % T) + lb_turnaround
                        if duration >= lb_turnaround and duration <= ub_turnaround:
                            C[line].add_edge((v, line, 'arr', t1, '-'), (v, line, 'dep', t2, '+'), 
                                            type='turnaround', duration=duration, weight=1)
    # transfer arcs
    for v in list(set([node for i in L.keys() for node in L[i]['nodes']])):
        for l1 in L.keys():
            for l2 in L.keys():
                if l1 != l2 and v in L[l1]['nodes'] and v in L[l2]['nodes']:
                    for t1 in range(T):
                        for t2 in range(T):
                            duration = ((t2 - t1 - lb_transfer) % T) + lb_transfer
                            if duration >= lb_transfer and duration <= ub_transfer:
                                if (v, l1, 'arr', t1, '+') in C[l1].nodes():
                                    if (v, l2, 'dep', t2, '+') in C[l2].nodes():
                                        A[l1, l2, v].add_edge((v, l1, 'arr', t1, '+'), (v, l2, 'dep', t2, '+'), 
                                                            type='transfer', duration=duration, weight=weights[(l1, l2, v)])
                                    if (v, l2, 'dep', t2, '-') in C[l2].nodes():
                                        A[l1, l2, v].add_edge((v, l1, 'arr', t1, '+'), (v, l2, 'dep', t2, '-'), 
                                                            type='transfer', duration=duration, weight=weights[(l1, l2, v)])
                                if (v, l1, 'arr', t1, '-') in C[l1].nodes():
                                    if (v, l2, 'dep', t2, '+') in C[l2].nodes():
                                        A[l1, l2, v].add_edge((v, l1, 'arr', t1, '-'), (v, l2, 'dep', t2, '+'), 
                                                            type='transfer', duration=duration, weight=weights[(l1, l2, v)])
                                    if (v, l2, 'dep', t2, '-') in C[l2].nodes():
                                        A[l1, l2, v].add_edge((v, l1, 'arr', t1, '-'), (v, l2, 'dep', t2, '-'), 
                                                            type='transfer', duration=duration, weight=weights[(l1, l2, v)])


    return C, A


def cxpesp (C, A, L, T, print_pesp=False, writeLP=False, vtype='C'):
    model = Model('cXPESP')
    model.hideOutput()
    model.setPresolve(SCIP_PARAMSETTING.OFF)
    model.setHeuristics(SCIP_PARAMSETTING.OFF)
    model.disablePropagation()

    x = {}
    z = {}
    cycles_all = cycles_in_graph (C, L, T)
    for line in L.keys():
        for c in range(cycles_all[line]['number_of_cycles']):
            x[line, c] = model.addVar(vtype=vtype, name=f'x_cxpesp_{line}_{c}', lb=0)

        model.addCons(quicksum(x[line, c] for c in range(cycles_all[line]['number_of_cycles'])) == 1)
    

    transfer_duration = {}
    transfer_weight = {}
    i = 0
    for a in A:
        transfer_duration[a] = nx.get_edge_attributes(A[a], 'duration')
        transfer_weight[a] = nx.get_edge_attributes(A[a], 'weight')
        for arc in A[a].edges():
            name = i
            i += 1
            z[arc] = model.addVar(vtype=vtype, name=f'z_{name}_{arc}', lb=0)
   
    for l1, l2, v in A.keys():
        for vt in A[l1,l2,v].nodes():
            for direction in ['+', '-']:
                if vt[2] == 'arr' and [w for w in A[l1,l2,v][vt] if w[4] == direction] != []:
                    model.addCons(quicksum(x[l1, c] for c in range(cycles_all[l1]['number_of_cycles']) if vt in cycles_all[l1][c]['nodes']) 
                                - quicksum(z[vt, w] for w in A[l1,l2,v][vt] if w[4] == direction) == 0)
                if vt[2] == 'dep' and [u for u in A[l1,l2,v].predecessors(vt) if u[4] == direction] != []:
                    model.addCons(quicksum(x[l2, c] for c in range(cycles_all[l2]['number_of_cycles']) if vt in cycles_all[l2][c]['nodes']) 
                                - quicksum(z[u, vt] for u in A[l1,l2,v].predecessors(vt) if u[4] == direction) == 0)


    model.setObjective(quicksum(x[line, c] * cycles_all[line][c]['duration'] * cycles_all[line][c]['weight'] for line, c in x.keys()) 
                        + quicksum(z[arc] * transfer_duration[a][arc] * transfer_weight[a][arc] for a in A for arc in A[a].edges()), 
                        sense='minimize')

    model.optimize()
    if writeLP:
        model.writeProblem(filename='written_lps/model.lp')

    obj_val = model.getObjVal()
    sol = model.getBestSol()
    #model.printBestSol()
    
    if print_pesp:
    	x_pesp, z_pesp = cxpesp_to_pesp (T, C, A, L, cycles_all, sol, x, z)
        print('\n' + 20*'=' + '\n Projection onto PESP: \n' + 20*'=' + '\n')
        print_results (x_pesp, 'x_pesp')
        print_results (z_pesp, 'z_pesp')

    nc, nt = number_of_variables(x, z)
    n_cons = len([constraint for constraint in model.getConss() if str(constraint).startswith('c')])

    return obj_val, nc, nt, n_cons


def cxpesp_to_pesp (T, C, A, L, cycles_all, cxpesp, x, z):
    x_pesp = defaultdict(int)
    z_pesp = defaultdict(int)

    for line in L:
        duration = nx.get_edge_attributes(C[line], 'duration')
        for c in range(cycles_all[line]['number_of_cycles']):
            for ((v1, line, s1, t1, d1), (v2, line, s2, t2, d2)) in cycles_all[line][c]['edges']:
                x_pesp[(v1, line, s1, d1), (v2, line, s2, d2)] += duration[((v1, line, s1, t1, d1), (v2, line, s2, t2, d2))] * cxpesp[x[line, c]]
    
    s = ('arr', 'dep')
    for l1, l2, v in A:
        duration = nx.get_edge_attributes(A[l1, l2, v], 'duration')
        for ((v, l1, s1, t1, d1), (v, l2, s2, t2, d2)) in A[l1, l2, v].edges():   
            if ((v, l1, s1, t1, d1), (v, l2, s2, t2, d2)) in z.keys():
                z_pesp[(v, l1, s1, d1), (v, l2, s2, d2)] += duration[((v, l1, s1, t1, d1), (v, l2, s2, t2, d2))] * cxpesp[z[(v, l1, s1, t1, d1), (v, l2, s2, t2, d2)]]

    return x_pesp, z_pesp



def cycle_attributes (line, cycle, C):
    edge_duration = nx.get_edge_attributes(C[line], 'duration')
    edge_weights = nx.get_edge_attributes(C[line], 'weight')
    duration = 0
    weight = 0
    for i in range(len(cycle)):
        duration += edge_duration[cycle[i-1], cycle[i]]
        weight += edge_weights[cycle[i-1], cycle[i]]
    weight = weight / len(cycle)
    return duration, weight


def rekursive_paths (graph, previous_paths, current_length, max_length):
    new_paths = []
    for path in previous_paths:
        for w in graph[path[-1]]:          
            new_paths.append(path + [w])
    if current_length + 1 < max_length:
        new_paths = rekursive_paths(graph, new_paths, current_length + 1, max_length)
    return new_paths


def find_simple_cycles (graph, nodes, max_length):
    all_simple_cycles = []
    start = [[v] for v in nodes]
    all_simple_paths = rekursive_paths(graph, start, 1, max_length + 1)
    for path in all_simple_paths:
        if path[0] == path[-1]:
            path.pop()
            all_simple_cycles.append(path)
    return all_simple_cycles


def cycles_in_graph (C, L, T):
    cycles_all = {}
    for line in L.keys():
        cycles_line = find_simple_cycles(C[line], [(L[line]['nodes'][0], line, 'dep', t, '+') for t in range(T) if (L[line]['nodes'][0], line, 'dep', t, '+') in C[line].nodes()] , 4 * (len(L[line]['nodes']) - 1))
        cycles_all[line] = {}
        cycles_all[line]['number_of_cycles'] = len(cycles_line)
        for c in range(cycles_all[line]['number_of_cycles']):          
            duration, weight = cycle_attributes (line, cycles_line[c], C)
            cycles_all[line][c] = {}
            cycles_all[line][c]['nodes'] = cycles_line[c]
            cycles_all[line][c]['edges'] = [(cycles_all[line][c]['nodes'][i-1], cycles_all[line][c]['nodes'][i]) for i in range(len(cycles_all[line][c]['nodes']))]
            cycles_all[line][c]['duration'] = duration
            cycles_all[line][c]['weight'] = weight
    return cycles_all


def one_cycle_per_line (C, L):
    cycles_all = {}
    for line in L.keys():
        path = [(L[line]['nodes'][0], line, 'dep', 0, '+')]
        for i in range(4 * (len(L[line]['nodes']) - 1) -1):
            path.append(list(C[line][path[-1]])[0])
        cycles_all[line] = {}
        cycles_all[line][0] = {}
        cycles_all[line][0]['nodes'] = path
        cycles_all[line][0]['edges'] = [(cycles_all[line][0]['nodes'][i-1], cycles_all[line][0]['nodes'][i]) for i in range(len(cycles_all[line][0]['nodes']))]
        cycles_all[line][0]['duration'] = cycle_attribute  (cycles_all[line][0]['edges'], C, line, attribute='duration')
        cycles_all[line][0]['weight'] = cycle_attribute  (cycles_all[line][0]['edges'], C, line, attribute='weight')
    return cycles_all


def cycle_attribute (cycle, C, line, attribute):
    edge_attribute = nx.get_edge_attributes(C[line], attribute)
    res = 0
    for arc in cycle:
        res += edge_attribute[arc]
    if attribute == 'weight':
        res = res / len(cycle)
    return res


def print_results (result_dict, name):
    print(f'Values for {name}: ')
    for key in result_dict.keys():
        print(f'{name}_{key} = {result_dict[key]}')


def number_of_variables(x, z):
    number_of_cycles = len(x.keys())
    number_transfer_arcs = len(z.keys())
    return number_of_cycles, number_transfer_arcs


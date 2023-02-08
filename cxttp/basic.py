# 09.02.2023
# Stephanie Riedmüller

from pyscipopt import Model, quicksum , SCIP_PARAMSETTING              # type: ignore
import networkx as nx
from collections import defaultdict
from cxpesp.basic import cycles_in_graph, number_of_variables


def expanded_graph_pr(L, T, weights):

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
    A = nx.DiGraph()
    G = nx.DiGraph()
    OD_arcs = nx.DiGraph()
    for line in L.keys():
    # driving arcs
        for (v,w) in L[line]['edges']:
            for t1 in range(T):
                for t2 in range(T):
                    duration = ((t2 - t1 - lb_driving[line, (v,w)]) % T) + lb_driving[line, (v,w)]
                    if duration >= lb_driving[line, (v,w)] and duration <= ub_driving[line, (v,w)]:
                        C[line].add_edge((v, line, 'dep', t1, '+'), (w, line, 'arr', t2, '+'), duration=duration, weight=1)
                        C[line].add_edge((w, line, 'dep', t1, '-'), (v, line, 'arr', t2, '-'), duration=duration, weight=1)
                        G.add_edge((v, line, 'dep', t1, '+'), (w, line, 'arr', t2, '+'), duration=duration, weight=1)
                        G.add_edge((w, line, 'dep', t1, '-'), (v, line, 'arr', t2, '-'), duration=duration, weight=1)
    # waiting and turnaround arcs
        for v in L[line]['nodes']:
            for t1 in range(T):
                for t2 in range(T):
                    if (v, line, 'arr', t1, '+') in C[line].nodes():
                        if (v, line, 'dep', t2, '+') in C[line].nodes():
                            duration = ((t2 - t1 - lb_waiting) % T) + lb_waiting
                            if duration >= lb_waiting and duration <= ub_waiting:
                                C[line].add_edge((v, line, 'arr', t1, '+'), (v, line, 'dep', t2, '+'), duration=duration, weight=1)
                                C[line].add_edge((v, line, 'arr', t1, '-'), (v, line, 'dep', t2, '-'), duration=duration, weight=1)
                                G.add_edge((v, line, 'arr', t1, '+'), (v, line, 'dep', t2, '+'), duration=duration, weight=1)
                                G.add_edge((v, line, 'arr', t1, '-'), (v, line, 'dep', t2, '-'), duration=duration, weight=1)
                        else:
                            duration = ((t2 - t1 - lb_turnaround) % T) + lb_turnaround
                            if duration >= lb_turnaround and duration <= ub_turnaround:
                                C[line].add_edge((v, line, 'arr', t1, '+'), (v, line, 'dep', t2, '-'), duration=duration, weight=1)
                                G.add_edge((v, line, 'arr', t1, '+'), (v, line, 'dep', t2, '-'), duration=duration, weight=1)
                    else:
                        duration = ((t2 - t1 - lb_turnaround) % T) + lb_turnaround
                        if duration >= lb_turnaround and duration <= ub_turnaround:
                            C[line].add_edge((v, line, 'arr', t1, '-'), (v, line, 'dep', t2, '+'), duration=duration, weight=1)
                            G.add_edge((v, line, 'arr', t1, '-'), (v, line, 'dep', t2, '+'), duration=duration, weight=1)
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
                                        G.add_edge((v, l1, 'arr', t1, '+'), (v, l2, 'dep', t2, '+'), duration=duration, weight=weights[(l1, l2, v)])
                                        A.add_edge((v, l1, 'arr', t1, '+'), (v, l2, 'dep', t2, '+'), duration=duration, weight=weights[(l1, l2, v)])
                                    if (v, l2, 'dep', t2, '-') in C[l2].nodes():
                                        G.add_edge((v, l1, 'arr', t1, '+'), (v, l2, 'dep', t2, '-'), duration=duration, weight=weights[(l1, l2, v)])
                                        A.add_edge((v, l1, 'arr', t1, '+'), (v, l2, 'dep', t2, '-'), duration=duration, weight=weights[(l1, l2, v)])
                                if (v, l1, 'arr', t1, '-') in C[l1].nodes():
                                    if (v, l2, 'dep', t2, '+') in C[l2].nodes():
                                        G.add_edge((v, l1, 'arr', t1, '-'), (v, l2, 'dep', t2, '+'), duration=duration, weight=weights[(l1, l2, v)])
                                        A.add_edge((v, l1, 'arr', t1, '-'), (v, l2, 'dep', t2, '+'), duration=duration, weight=weights[(l1, l2, v)])
                                    if (v, l2, 'dep', t2, '-') in C[l2].nodes():
                                        G.add_edge((v, l1, 'arr', t1, '-'), (v, l2, 'dep', t2, '-'), duration=duration, weight=weights[(l1, l2, v)])
                                        A.add_edge((v, l1, 'arr', t1, '-'), (v, l2, 'dep', t2, '-'), duration=duration, weight=weights[(l1, l2, v)])
    # arcs for OD-pairs
    for line in L.keys():
        for v in L[line]['nodes']:
            for t in range(T):
                for d in ['+', '-']:
                    if (v, line, 'dep', t, d) in C[line].nodes():
                        G.add_edge((v, 'dep'), (v, line, 'dep', t, d), duration=0, weight=0)
                        OD_arcs.add_edge((v, 'dep'), (v, line, 'dep', t, d), duration=0, weight=0)
                    if (v, line, 'arr', t, d) in C[line].nodes():
                        G.add_edge((v, line, 'arr', t, d), (v, 'arr'), duration=0, weight=0)
                        OD_arcs.add_edge((v, line, 'arr', t, d), (v, 'arr'), duration=0, weight=0)

    return C, A, G, OD_arcs


def cxttp (G, C, ODpairs, L, T, cutoff, writeLP=False, pesp=False, vtype='C'):
    model = Model('cXPESP_pr')
    model.hideOutput()
    model.setPresolve(SCIP_PARAMSETTING.OFF)
    model.setHeuristics(SCIP_PARAMSETTING.OFF)
    model.disablePropagation()
    x = dict()
    cycles_all = cycles_in_graph (C, L, T)
    for line in L.keys():
        for c in range(cycles_all[line]['number_of_cycles']):
            cc = cycles_all[line][c]['duration']
            x[line, c] = model.addVar(vtype=vtype, name=f'x_cxpesp_{line}_{c}_{cc}', lb=0, obj=cycles_all[line][c]['duration'] * cycles_all[line][c]['weight'])

        model.addCons(quicksum(x[line, c] for c in range(cycles_all[line]['number_of_cycles'])) == 1)
    
    paths = dict()
    y = dict()
    for s, t in ODpairs.keys():
        paths[s, t] = find_st_paths(G, s, t, cutoff)
        for p in paths[s, t].keys():
            pp = paths[s, t][p]['nodes']
            y[s, t, p] = model.addVar(vtype='C', name=f'y_cxpesp_{s}_{t}_{p}_{pp}', lb=0, obj = ODpairs[s, t] * paths[s, t][p]['duration'])

        model.addCons(quicksum(y[s, t, i] for i in paths[s, t].keys()) == 1)


    for line in C.keys():     
        for arc in C[line].edges():
            for s, t in ODpairs.keys():
                model.addCons(
                    quicksum(x[line, c] for c in range(cycles_all[line]['number_of_cycles']) if arc in cycles_all[line][c]['edges'])
                    - quicksum(y[s, t, p] for p in paths[s, t].keys() if arc in paths[s, t][p]['arcs']) 
                    >= 0)


    model.optimize()
    if writeLP:
        model.writeProblem(filename='written_lps/model.lp')

    obj_val = model.getObjVal()
    #model.printBestSol()

    nc, nt = number_of_variables(x, y)

    n_cons = len([constraint for constraint in model.getConss() if str(constraint).startswith('c')])

    if pesp:
        sol = model.getBestSol()
        print('Projection onto PESP')
        x_pesp = projection_onto_pesp (C, L, cycles_all, sol, x)
        for key, val in x_pesp.items():
            print(key, val)

    return obj_val, nc, nt, n_cons


def find_st_paths(G, s, t, cutoff):
    st_paths = dict()
    # The cutoff must be a reasonable passenger path length depending on the instance
    simple_paths = nx.all_simple_paths(G=G, source=(s,'dep'), target=(t,'arr'), cutoff=cutoff)
    for i, path in enumerate(simple_paths):
        st_paths[i] = dict()
        st_paths[i]['nodes'] = path
        st_paths[i]['arcs'] = [(path[i], path[i+1]) for i in range(len(path) - 1)]
        arc_durations = nx.get_edge_attributes(G, 'duration')
        st_paths[i]['duration'] = sum([arc_durations[arc] for arc in st_paths[i]['arcs']])
    return st_paths


def projection_onto_pesp (C, L, cycles_all, cxpesp, x):
    x_pesp = defaultdict(int)

    for line in L:
        duration = nx.get_edge_attributes(C[line], 'duration')
        for c in range(cycles_all[line]['number_of_cycles']):
            for ((v1, line, s1, t1, d1), (v2, line, s2, t2, d2)) in cycles_all[line][c]['edges']:
                x_pesp[(v1, line, s1, d1), (v2, line, s2, d2)] += duration[((v1, line, s1, t1, d1), (v2, line, s2, t2, d2))] * cxpesp[x[line, c]]
    
    return x_pesp



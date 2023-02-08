# 09.02.2023
# Stephanie Riedmüller

from pyscipopt import Model, quicksum               # type: ignore
import networkx as nx
from collections import defaultdict

def expanded_graph_pesp(L, T, weights):
    lb_driving = {}
    ub_driving = {}
    if "duration" in L[0].keys():
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
            C[line].add_edge((v, line, 'dep', '+'), (w, line, 'arr', '+'), 
                            type='driving', weight=1, lb = lb_driving[line, (v,w)], ub = ub_driving[line, (v,w)])
            C[line].add_edge((w, line, 'dep', '-'), (v, line, 'arr', '-'), 
                            type='driving', weight=1, lb = lb_driving[line, (v,w)], ub = ub_driving[line, (v,w)])
    # waiting and turnaround arcs
        for v in L[line]['nodes']:
            for d1 in ['+', '-']:
                for d2 in ['+', '-']:
                    if (v, line, 'arr', d1) in C[line].nodes():
                        if (v, line, 'dep', d2) in C[line].nodes():
                            if d1 == d2:
                                arc_type = 'waiting'
                                lb = lb_waiting
                                ub = ub_waiting
                                C[line].add_edge((v, line, 'arr', d1), (v, line, 'dep', d2), 
                                            type=arc_type, weight=1, lb=lb, ub=ub)
                            elif (v, line, 'arr', d2) not in C[line].nodes():
                                arc_type = 'turnaround'
                                lb = lb_turnaround
                                ub = ub_turnaround
                                C[line].add_edge((v, line, 'arr', d1), (v, line, 'dep', d2), 
                                            type=arc_type, weight=1, lb=lb, ub=ub)
    # transfer arcs
    for v in list(set([node for i in L.keys() for node in L[i]['nodes']])):
        for l1 in L.keys():
            for l2 in L.keys():
                if l1 != l2 and v in L[l1]['nodes'] and v in L[l2]['nodes']:
                    if (v, l1, 'arr', '+') in C[l1].nodes():
                        if (v, l2, 'dep', '+') in C[l2].nodes():
                            A[l1, l2, v].add_edge((v, l1, 'arr', '+'), (v, l2, 'dep', '+'), 
                                                type='transfer', weight=weights[(l1, l2, v)], lb=lb_transfer, ub=ub_transfer)
                        if (v, l2, 'dep', '-') in C[l2].nodes():
                            A[l1, l2, v].add_edge((v, l1, 'arr', '+'), (v, l2, 'dep', '-'), 
                                                type='transfer', weight=weights[(l1, l2, v)], lb=lb_transfer, ub=ub_transfer)
                    if (v, l1, 'arr', '-') in C[l1].nodes():
                        if (v, l2, 'dep', '+') in C[l2].nodes():
                            A[l1, l2, v].add_edge((v, l1, 'arr', '-'), (v, l2, 'dep', '+'), 
                                                type='transfer', weight=weights[(l1, l2, v)], lb=lb_transfer, ub=ub_transfer)
                        if (v, l2, 'dep', '-') in C[l2].nodes():
                            A[l1, l2, v].add_edge((v, l1, 'arr', '-'), (v, l2, 'dep', '-'), 
                                                type='transfer', weight=weights[(l1, l2, v)], lb=lb_transfer, ub=ub_transfer)

    return C, A


def pesp (C, A, L, T, vtype='I'):
    model = Model('PESP')
    model.hideOutput()

    x = {}
    z = {}
    pi = {}
    p = {}
    eps = 0.5

    for line in L.keys():
        lb = nx.get_edge_attributes(C[line], 'lb')
        ub = nx.get_edge_attributes(C[line], 'ub')
        for e in C[line].edges():
            x[e] = model.addVar(vtype=vtype, name=f'x_pesp_{e}', lb=lb[e], ub=ub[e])
            p[e] = model.addVar(vtype=vtype, name=f'p_{e}', lb=0)
        for v in C[line].nodes():
            pi[v] = model.addVar(vtype='C', name=f'pi_{v}', lb=0, ub=T - eps)
        for v, w in C[line].edges():
            model.addCons(pi[w] - pi[v] + (T * p[v, w]) - x[v, w] == 0)
    
    for a in A:
        lb = nx.get_edge_attributes(A[a], 'lb')
        ub = nx.get_edge_attributes(A[a], 'ub')
        for e in A[a].edges():
            z[e] = model.addVar(vtype=vtype, name=f'z_{e}', lb=lb[e], ub=ub[e])
            p[e] = model.addVar(vtype=vtype, name=f'p_{e}', lb=0)
            v, w = e
            model.addCons(pi[w] - pi[v] + (T * p[v, w]) - z[v, w] == 0)

    arc_weights = {}
    transfer_weight = {}
    edge_type = {}
    for line in L:
        arc_weights[line] = nx.get_edge_attributes(C[line], 'weight')
        edge_type[line] = nx.get_edge_attributes(C[line], 'type')
    for a in A:
        transfer_weight[a] = nx.get_edge_attributes(A[a], 'weight')
 
    model.setObjective(quicksum(x[e] * arc_weights[line][e]  for line in L for e in C[line].edges()) 
                    + quicksum(z[e] * transfer_weight[a][e]  for a in A for e in A[a].edges()), 
                    sense='minimize')


    model.optimize()
    obj_val = model.getObjVal()

    n_cons = len([constraint for constraint in model.getConss() if str(constraint).startswith('c')])
    n_vars = len(x.keys()) + len(z.keys())

    return obj_val, n_cons, n_vars


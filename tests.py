# 09.02.2023
# Stephanie Riedmüller

from collections import defaultdict
from pyscipopt import Model              # type: ignore
from time import time
import pandas as pd
import json
import random

from pesp.pesp import expanded_graph_pesp, pesp
from cxpesp.basic import expanded_graph, cxpesp
from cxpesp.column_generation import cxpesp_cg
from cxttp.basic import expanded_graph_pr, cxttp
from cxttp.column_generation import cxttp_cg
from xttp.xttp import xttp


def read_lines(filepath):
    L = {}
    i = 0
    with open(filepath, 'r') as f:
        line = f.readline().replace('\n', ' ').split(';')
        while line != ['']:
            if not line[0].startswith('#'):
                L[i] = {}
                L[i]['nodes'] = []              
                for v in line:
                    L[i]['nodes'].append(int(v))
                L[i]['edges'] = [(L[i]['nodes'][j], L[i]['nodes'][j+1]) for j in range(len(L[i]['nodes'])-1)]
                i += 1
            line = f.readline().replace('\n', ' ').split(';')
    return L


def read_json(filepath):
    with open(filepath, 'r') as f:
        L = json.load(f)
        keys = [x for x in range(len(L.keys()))]
        L = dict(zip(keys, list(L.values())))
        for line in L.keys():
            nodes = L[line]['nodes']
            L[line]['edges'] = [(nodes[j], nodes[j+1]) for j in range(len(nodes)-1)]
    return L


def read_problem_from_file(name_file):
    model = Model('read')
    model.readProblem(name_file)
    model.optimize()
    sol = model.getBestSol()
    for var in model.getVars():
        if sol[var] > 0:
            print(f'{var.name}:{sol[var]}')
    print(f'ObjVal: {model.getObjVal()}')
    return 


def random_weights (activities, min_weight, max_weight):
    weights = dict()
    for activity in activities:
            weights[activity] = random.randint(min_weight, max_weight)
    return weights


def activities (L):
    res = list()
    for v in list(set([node for i in L.keys() for node in L[i]['nodes']])):
        for l1 in L.keys():
            for l2 in L.keys():
                if l1 != l2 and v in L[l1]['nodes'] and v in L[l2]['nodes']:
                    res.append((l1,l2,v))
    return res


def generate_weights (T, L, randomized_weights=False):
    weights = defaultdict(lambda: 1, {})

    C, A, G, OD_arcs = expanded_graph_pr(L, T, weights=weights)
    C, G = break_symmetry_cxttp(C, G, T)

    stations = list(set([node for i in L.keys() for node in L[i]['nodes']]))
    ODpairs = {(s, t) : random.randint(1,5) for s in stations for t in stations if s != t}

    if randomized_weights == True:
        weights = random_weights (activities(L), 1, 5)
    
    else:
        _ , _ , _ , _ , y  = cxttp_cg(T=2, L=L, C=C, G=G, A=A, OD_arcs=OD_arcs, ODpairs=ODpairs, vtype='C', output_y=True)

        relevant_activities = activities(L)
        
        for (l1, l2, v) in relevant_activities:
            weight = 1
            for (s,t) in ODpairs.keys():
                for path in y[s,t]:
                    for (u1, u2) in  path[1]:
                        if u1[0] == v and u2[0] == v and u1[1] == l1 and u2[1] == l2:
                            weight += path[0]
            weights[(l1, l2, v)] = weight
        
    return weights, ODpairs


def break_symmetry_cxpesp (C, T):
    v, k, s, t, d = list(C[0].nodes())[0]
    C[0].remove_nodes_from([(v, k, s, timestep, d) for timestep in range(1, T)])
    return C


def break_symmetry_cxttp (C, G, T):
    v, k, s, t, d = list(C[0].nodes())[0]
    C[0].remove_nodes_from([(v, k, s, timestep, d) for timestep in range(1, T)])
    G.remove_nodes_from([(v, k, s, timestep, d) for timestep in range(1, T)])
    return C, G


# ------------


def test_pesp(T_lst, L, vtype, name, weight):
    df = pd.DataFrame(columns=['Obj_PESP','cons', 'vars', 'Time'])
    for T in T_lst:
        C_pesp, A_pesp = expanded_graph_pesp(L, T, weights=weight)
        t1 = time()
        obj_val, n_cons, n_vars = pesp (C_pesp, A_pesp, L, T, vtype=vtype)
        t2 = time() 
        t3 = round(t2-t1, 3)
        df.loc[T, :] = [obj_val, n_cons, n_vars, t3 ]
        df.to_csv(f'results/{name}.csv')
    print(df)
    return df


def test_cxpesp(T_lst, L, vtype, name, weight):
    df = pd.DataFrame(columns=['Obj_cXPESP', '#cycles', '#transfers', '#constraints', 'time_cXPESP'])
    for T in T_lst:
        C, A = expanded_graph(L, T, weights=weight)
        C = break_symmetry_cxpesp(C, T)
        t1 = time()
        obj_val_cxpesp, nc, nt, n_cons = cxpesp (C, A, L, T,vtype=vtype)
        t2 = time()
        t3  = round(t2-t1, 3)
        df.loc[T, :] = [obj_val_cxpesp, nc, nt, n_cons, t3 ]
        df.to_csv(f'results/{name}.csv')
    print(df)
    return df


def test_cxpesp_cg(T_lst, L, vtype, name, weight):
    df = pd.DataFrame(columns=['Obj_cg', '#cycles', '#transfers', '#constraints', 'time'])
    for T in T_lst:
        C, A = expanded_graph(L, T, weights=weight)
        C = break_symmetry_cxpesp(C, T)
        t1 = time()
        obj_val_cxpesp_cg_pricer, nc_cg_pricer, nz_cg_pricer, n_cons2  = cxpesp_cg(T, L, C, A, vtype=vtype)
        t2 = time()
        t3  = round(t2-t1, 3)
        df.loc[T, :] = [obj_val_cxpesp_cg_pricer, nc_cg_pricer, nz_cg_pricer, n_cons2, t3 ]
        df.to_csv(f'results/{name}.csv')
    print(df)
    return df


def test_xttp(T_lst, L, ODpairs, vtype, cutoff, name, symmetrybreak, weight):
    df = pd.DataFrame(columns=['Obj_XTTP', 'n_vars_x', 'n_vars_y', 'n_cons', 'time'])
    for T in T_lst:        
        C, A, G, OD_arcs = expanded_graph_pr(L, T, weights=weight)
        if symmetrybreak == True:
            C, G = break_symmetry_cxttp(C, G, T)
        t1 = time()
        obj_val, n_cons, n_vars_x, n_vars_y = xttp (G=G, C=C, L=L, ODpairs=ODpairs, vtype=vtype, cutoff=cutoff)
        t2 = time()
        t3  = round(t2-t1, 3)
        df.loc[T, :] = [obj_val, n_vars_x, n_vars_y, n_cons, t3 ]
        df.to_csv(f'results/{name}.csv')
    print(df)
    return df


def test_cxttp(T_lst, L, ODpairs, vtype, cutoff, name, weight):
    df = pd.DataFrame(columns=['Obj_cXTTP', '#cycles', '#y', '#constraints', 'time'])
    for T in T_lst:
        C, A, G, OD_arcs = expanded_graph_pr(L, T, weights=weight)
        C, G = break_symmetry_cxttp(C, G, T)
        t1 = time()
        obj_val_cxttp, nc, nt, n_cons = cxttp (G=G, C=C, ODpairs=ODpairs, L=L, T=T, vtype=vtype, cutoff=cutoff)
        t2 = time()
        t3  = round(t2-t1, 3)
        df.loc[T, :] = [obj_val_cxttp, nc, nt, n_cons, t3 ]
        df.to_csv(f'results/{name}.csv')
    print(df)
    return df


def test_cxttp_cg(T_lst, L, ODpairs, vtype, name, weight):
    df = pd.DataFrame(columns=['Obj_cg', '#cycles', '#paths', '#constraints', 'time'])
    for T in T_lst:
        C, A, G, OD_arcs = expanded_graph_pr(L, T, weights=weight)
        C, G = break_symmetry_cxttp(C, G, T)
        t1 = time()
        obj_val, nc, nz, n_cons  = cxttp_cg(T=T, L=L, C=C, G=G, A=A, OD_arcs=OD_arcs, ODpairs=ODpairs, vtype=vtype)
        t2 = time()
        t3  = round(t2-t1, 3)
        df.loc[T, :] = [obj_val, nc, nz, n_cons, t3 ]
        df.to_csv(f'results/{name}.csv')
    print(df)
    return df



###########################################################


L = read_lines('instances/input_2cross.txt')
#L = read_json('instances/input_3berlin.json')
#L = read_json('instances/input_berlin.json')

weight , ODpairs = generate_weights (T=2, L=L, randomized_weights=True, name='random')
#weight = {'arcs': defaultdict(lambda: 1, {}), 'transfers': defaultdict(lambda: 1, {})}

period = [5,10,15,20,30,40,50,60]

df_pesp_IP = test_pesp(period, L, vtype='I', weight=weight)
df_pesp_LP = test_pesp(period, L, vtype='C', weight=weight)
df_cxpesp_LP = test_cxpesp(period, L, vtype='C', weight=weight)
df_cxpesp_cg = test_cxpesp_cg(period, L, vtype='C', weight=weight)

cutoff = 6      # suitable maximal passenger path length, dependent on instance
df_xttp_IP = test_xttp(period, L, ODpairs, vtype='I', cutoff=cutoff, symmetrybreak=True, weight=weight)
df_xttp_LP_1 = test_xttp(period, L, ODpairs, vtype='C', cutoff=cutoff, symmetrybreak=True, weight=weight)
df_xttp_LP_2 = test_xttp(period, L, ODpairs, vtype='C', cutoff=cutoff, symmetrybreak=False, weight=weight)
df_cxttp_LP = test_cxttp(period, L, ODpairs, vtype='C', cutoff=cutoff, weight=weight)
df_cxttp_cg = test_cxttp_cg(period, L, ODpairs, vtype='C', weight=weight)




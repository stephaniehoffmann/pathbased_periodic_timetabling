# 09.02.2023
# Stephanie Riedmüller

from pyscipopt import Model, quicksum, SCIP_PARAMSETTING              # type: ignore
import networkx as nx
from cxttp.basic import find_st_paths



def xttp (G, C, L, ODpairs, vtype='I', cutoff=10):
    model = Model('XTTP')
    model.hideOutput()
    model.setPresolve(SCIP_PARAMSETTING.OFF)
    model.setHeuristics(SCIP_PARAMSETTING.OFF)
    model.disablePropagation()

    x = dict()
    y = dict()
    paths = dict()

    for line in L.keys():
        arc_durations = nx.get_edge_attributes(C[line], 'duration')
        arc_weight = nx.get_edge_attributes(C[line], 'weight')
        for e in C[line].edges():
            x[e] = model.addVar(vtype=vtype, name=f'x_xttp_{e}', lb=0, ub=1, obj=arc_durations[e] * arc_weight[e])

        
        for v, u in L[line]['edges']:
            model.addCons(quicksum(x[e] for e in C[line].edges() if e[0][0] == v and e[1][0] == u) == 1)
            model.addCons(quicksum(x[e] for e in C[line].edges() if e[0][0] == u and e[1][0] == v) == 1)


        for v in C[line].nodes():
            model.addCons(quicksum(x[v,u] for u in C[line][v]) - quicksum(x[u,v] for u in C[line].predecessors(v)) == 0)

    
    for s, t in ODpairs.keys():
        paths[s, t] = find_st_paths(G, s, t, cutoff=cutoff)
        for p in paths[s, t].keys():
            y[s, t, p] = model.addVar(vtype='C', name=f'y_cxpesp_{s}_{t}_{p}', lb=0, obj = ODpairs[s, t] * paths[s, t][p]['duration'])

        model.addCons(quicksum(y[s, t, i] for i in paths[s, t].keys()) == 1)

    for (s,t) in ODpairs:
        for line in L.keys():
            for e in C[line].edges():
                relevant_paths = [p for p in paths[s, t].keys() if e in paths[s, t][p]['arcs']]
                if len(relevant_paths) > 0:
                    model.addCons(x[e] - quicksum(y[s, t, p] for p in relevant_paths) >= 0)
    

    model.optimize()
    obj_val = model.getObjVal()
    #model.printBestSol()

    n_cons = len([constraint for constraint in model.getConss() if str(constraint).startswith('c')])
    n_vars_x = len(x.keys())
    n_vars_y = len(y.keys())

    return obj_val, n_cons, n_vars_x, n_vars_y


# 09.02.2023
# Stephanie Riedmüller

from pyscipopt import Model, Pricer, SCIP_RESULT, quicksum       # type: ignore
from cxpesp.basic import cycle_attribute, one_cycle_per_line
import networkx as nx
from copy import deepcopy
from collections import defaultdict


class CyclePricer(Pricer):

    # The reduced cost function for the variable pricer
    def pricerredcost(self):


        # Retrieving the dual solutions
        dualSolutions = {'mu' : list(), 'pi' : defaultdict(int), 'tau' : defaultdict(int)}
        for idx, cons in enumerate(self.data['cons']['partitionCons']):
            dualSolutions['mu'].append(self.model.getDualsolLinear(cons))
        for idx, cons in enumerate(self.data['cons']['transferConsOut']):       
            dualSolutions['pi'][self.data['transfer_nodes']['out'][idx]] = self.model.getDualsolLinear(cons)             
        for idx, cons in enumerate(self.data['cons']['transferConsIn']):       
            dualSolutions['tau'][self.data['transfer_nodes']['in'][idx]] = self.model.getDualsolLinear(cons)
        
        # Finding a new cycle variable
        for line in self.data['lines'].keys():
            tmp_graph = deepcopy(self.data['graph']['C'][line])
            duration = nx.get_edge_attributes(self.data['graph']['C'][line], 'duration')

            v1, v2 = self.data['lines'][line]['edges'][0]
            for t in range(self.data['T']):
                v1t = (v1, line, 'dep', t, '+')
                if v1t in self.data['graph']['C'][line].nodes():
                    arcs = self.data['graph']['C'][line][(v1, line, 'dep', t, '+')]
                    for v2t in arcs:
                        _, l, s, t2, d = v2t
                        tmp_graph.remove_edge(v1t, v2t ) 
                        tmp_graph.add_edge(v1t, (f'artificial_{v2}', l, s, t2, d))
                        duration[v1t, (f'artificial_{v2}', l, s, t2, d)] = duration[v1t, v2t]


            reduced_costs = dict()
            line_weight = self.data['cycles'][line][0]['weight']
            for u1, u2 in tmp_graph.edges():
                rc = duration[u1, u2] * line_weight
                for line2 in self.data['lines'].keys():
                    if line != line2:
                        rc = (rc    - dualSolutions['pi'][(u1, '+', line, line2)] 
                                    - dualSolutions['pi'][(u1, '-', line, line2)] 
                                    - dualSolutions['tau'][(u2, '+',line2, line)] 
                                    - dualSolutions['tau'][(u2, '-', line2, line)] )
                reduced_costs[(u1, u2)] = {"reduced_costs" : rc}
        
            nx.set_edge_attributes(tmp_graph, values=reduced_costs)           
            

            for t in range(self.data['T']):
                source = (v2, l, s, t, d)
                target = (f'artificial_{v2}', l, s, t, d)
                if source in tmp_graph.nodes() and target in tmp_graph.nodes():
                    path_nodes = nx.shortest_path(tmp_graph, source=source, target=target, weight="reduced_costs", method="bellman-ford")
                    path_arcs = [(path_nodes[i], path_nodes[i+1]) for i in range(len(path_nodes) - 1)]
                    reduced_path_cost = 0
                    for arc in path_arcs:
                        reduced_path_cost += reduced_costs[arc]["reduced_costs"]
            
                    eps = 1e-08
                    if dualSolutions['mu'][line] > reduced_path_cost + eps and len(path_nodes) > 0:
                        path_nodes.pop()
                        cycle = [(path_nodes[i-1], path_nodes[i]) for i in range(len(path_nodes))]
                        currentNumVar = len(self.data['var']['cycleVars'][line])

                        self.data['cycles'][line].append({  'edges' : cycle, 
                                                            'nodes' : path_nodes,
                                                            'duration' : cycle_attribute(cycle, self.data['graph']['C'], line, attribute='duration'), 
                                                            'weight' : cycle_attribute(cycle, self.data['graph']['C'], line, attribute='weight')})

                        newVar = self.model.addVar(vtype=self.data['vtype'], name=f'Cycle_{line}_{currentNumVar}', lb=0, 
                                                    obj= self.data['cycles'][line][currentNumVar]['duration'] * self.data['cycles'][line][currentNumVar]['weight'], 
                                                    pricedVar = True)

                        self.data['var']['cycleVars'][line].append(newVar)

                        self.model.addConsCoeff(self.data['cons']['partitionCons'][line], newVar, 1.0)

                        for i, c in enumerate(self.data['cons']['transferConsOut']):
                            if self.data['transfer_nodes']['out'][i][0] in self.data['cycles'][line][currentNumVar]['nodes']:
                                self.model.addConsCoeff(c, newVar, 1.0)
                    
                        for i, c in enumerate(self.data['cons']['transferConsIn']):
                            if self.data['transfer_nodes']['in'][i][0] in self.data['cycles'][line][currentNumVar]['nodes']:
                                self.model.addConsCoeff(c, newVar, 1.0)

    
        return {'result':SCIP_RESULT.SUCCESS}

    # The initialisation function for the variable pricer to retrieve the transformed constraints of the problem
    def pricerinit(self):
        for cons_type in ['partitionCons', 'transferConsOut', 'transferConsIn']:
            for i, c in enumerate(self.data['cons'][cons_type]):
                self.data['cons'][cons_type][i] = self.model.getTransformedCons(c)


def cxpesp_cg(T, L, C, A, vtype='C'):
    model = Model("cXPESP")

    model.setPresolve(0)
    model.hideOutput()

    # creating a pricer
    pricer = CyclePricer()
    model.includePricer(pricer, "CyclePricer", "Pricer to identify new cycles in line")

    # adding the initial variables to the restricted master problem
    cycles_subset = one_cycle_per_line(C, L)
    
    cycleVars = defaultdict(list)
    cycles = defaultdict(list)

    for line in L.keys():
        cycleVars[line].append(model.addVar(vtype=vtype, name=f'Cycle_{line}_{0}', lb=0, obj= cycles_subset[line][0]['duration'] * cycles_subset[line][0]['weight'])) 
        cycles[line].append(cycles_subset[line][0])

    transferVars = defaultdict(dict)
    transfer_duration = dict()
    transfer_weight = dict()

    for a in A.keys():
        transfer_duration[a] = nx.get_edge_attributes(A[a], 'duration')
        transfer_weight[a] = nx.get_edge_attributes(A[a], 'weight')
        for j, arc in enumerate(A[a].edges()):
            transferVars[a][arc] = model.addVar(vtype=vtype, name=f'Transfer_{a}_{j}', lb=0, obj= transfer_duration[a][arc] * transfer_weight[a][arc])


    # adding constraints
    partitionCons = list()

    for line in L.keys():
        partitionCons.append(model.addCons(
            quicksum(cycleVars[line][i] for i in range(len(cycleVars[line]))) == 1, 
            name=f'Partition_{line}', separate = False, modifiable = True))
        
    
    transferCons = {'out' : list(), 'in' : list()}
    nodeTransferCons = {'out' : list(), 'in' : list()}

    for a in A.keys():
        (l1, l2, v) = a     
        for vt in A[l1,l2,v].nodes():
            for direction in ['+', '-']:
                if vt[2] == 'arr':
                    t_lst = [w for w in A[l1,l2,v][vt] if w[4] == direction]
                    if t_lst:
                        nodeTransferCons['out'].append((vt, direction, l1, l2))
                        transferCons['out'].append(model.addCons(quicksum(cycleVars[l1][i] for i in [0] if vt in cycles[l1][i]['nodes']) 
                                    - quicksum(transferVars[a][vt, w] for w in t_lst) == 0, 
                                    name=f'Coupling_out_{vt}_{direction}', separate = False, modifiable = True))
                if vt[2] == 'dep':
                    t_lst = [u for u in A[l1,l2,v].predecessors(vt) if u[4] == direction]
                    if t_lst:
                        nodeTransferCons['in'].append((vt, direction, l1, l2))
                        transferCons['in'].append(model.addCons(quicksum(cycleVars[l2][i] for i in [0] if vt in cycles[l2][i]['nodes']) 
                                    - quicksum(transferVars[a][u, vt] for u in t_lst) == 0, 
                                    name=f'Coupling_in_{vt}_{direction}', separate = False, modifiable = True))


    # Setting the pricer_data for use in the init and redcost functions
    pricer.data = dict()
    pricer.data['var'] = {'cycleVars' : cycleVars, 'transferVars' : transferVars }
    pricer.data['cons'] = {'partitionCons' : partitionCons, 'transferConsOut' : transferCons['out'], 'transferConsIn' : transferCons['in'] }
    pricer.data['transfer_nodes'] = {'out': nodeTransferCons['out'], 'in': nodeTransferCons['in']}
    pricer.data['graph'] = {'C' : C, 'A' : A}
    pricer.data['cycles'] = cycles
    pricer.data['lines'] = L
    pricer.data['T'] = T
    pricer.data['vtype'] = vtype

    model.optimize()
    #model.printBestSol()

    n_cycles = sum([len(pricer.data['var']['cycleVars'][line]) for line in L.keys()] )
    n_transfers = sum([len(pricer.data['var']['transferVars'][a]) for a in A.keys() ])
    n_cons = sum([len(pricer.data['cons'][cons_type]) for cons_type in pricer.data['cons'].keys()])

    return model.getObjVal(), n_cycles, n_transfers, n_cons




# 09.02.2023
# Stephanie Riedmüller


from pyscipopt import Model, Pricer, SCIP_RESULT, SCIP_PARAMSETTING, quicksum       # type: ignore
from cxpesp.basic import cycle_attribute, one_cycle_per_line
import networkx as nx
from copy import deepcopy
from numpy import infty
from collections import defaultdict


class CyclePricer(Pricer):

    # The reduced cost function for the variable pricer
    def pricerredcost(self):

        # against rounding issues
        eps = 1e-08

        # Retrieving the dual solutions
        dualSolutions = {'mu' : defaultdict(int), 'pi' : defaultdict(int), 'tau' : defaultdict(int)}
        for r, cons in self.data['cons']['partitionCons'].items():
            dualSolutions['mu'][r] = self.model.getDualsolLinear(cons)
        for (s,t), cons in self.data['cons']['passengerflowCons'].items():       
            dualSolutions['pi'][s, t] = self.model.getDualsolLinear(cons)             
        for (arc, s, t), cons in self.data['cons']['couplingCons'].items():       
            dualSolutions['tau'][arc, s, t] = self.model.getDualsolLinear(cons)

        # finding new cycle variables
        for line in self.data['lines'].keys():
            tmp_graph = deepcopy(self.data['graph']['C'][line])
            duration = nx.get_edge_attributes(self.data['graph']['C'][line], 'duration')

            v1, v2 = self.data['lines'][line]['edges'][0]
            for t in range(self.data['T']):
                v1t = (v1, line, 'dep', t, '+')
                if v1t in self.data['graph']['C'][line].nodes():
                    arcs = self.data['graph']['C'][line][(v1, line, 'dep', t, '+')]
                    for v2t in arcs:
                        _, l, m, t2, d = v2t
                        tmp_graph.remove_edge(v1t, v2t ) 
                        tmp_graph.add_edge(v1t, (f'artificial_{v2}', l, m, t2, d))
                        duration[v1t, (f'artificial_{v2}', l, m, t2, d)] = duration[v1t, v2t]


            reduced_costs = dict()
            line_weight = self.data['cycles'][line][0]['weight']
            for arc in tmp_graph.edges():
                red_cost = duration[arc] * line_weight
                for s, t in self.data['ODpairs']:
                    red_cost -= dualSolutions['tau'][arc, s, t]
                reduced_costs[arc] = {"reduced_costs" : red_cost}
        
            nx.set_edge_attributes(tmp_graph, values=reduced_costs)           
            

            for t in range(self.data['T']):
                source = (v2, l, m, t, d)
                target = (f'artificial_{v2}', l, m, t, d)
                if source in tmp_graph.nodes() and target in tmp_graph.nodes():
                    path_nodes = nx.shortest_path(tmp_graph, source=source, target=target, weight="reduced_costs", method="bellman-ford")
                    path_arcs = [(path_nodes[i], path_nodes[i+1]) for i in range(len(path_nodes) - 1)]
                    reduced_path_cost = 0
                    for arc in path_arcs:
                        reduced_path_cost += reduced_costs[arc]["reduced_costs"]
            
                    
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

                    
                        for arc in cycle:
                            for source, target in self.data['ODpairs']:
                                if (arc, source, target) in self.data['cons']['couplingCons'].keys():
                                    c = self.data['cons']['couplingCons'][arc, source, target]
                                    self.model.addConsCoeff(c, newVar, 1.0)

      
        # finding new passenger path variables
        for s,t in self.data['ODpairs']:

            reduced_costs2 = dict()
            costs = nx.get_edge_attributes(self.data['graph']['G'], 'duration')
            for arc in self.data['graph']['G'].edges():
                reduced_costs2[arc] = {"reduced_costs2" : (self.data['ODpairs'][s,t] * costs[arc]) + dualSolutions['tau'][arc, s, t] + eps}
            
            nx.set_edge_attributes(self.data['graph']['G'], values=reduced_costs2)


            st_path_nodes = nx.shortest_path(self.data['graph']['G'], source=(s, 'dep'), target=(t, 'arr'), weight="reduced_costs2", method="dijkstra")
            st_path_arcs = [(st_path_nodes[i], st_path_nodes[i+1]) for i in range(len(st_path_nodes) - 1)]
            reduced_st_path_cost = 0
            st_path_duration = 0
            for arc in st_path_arcs:
                reduced_st_path_cost += reduced_costs2[arc]["reduced_costs2"]
                st_path_duration += costs[arc]


            if dualSolutions['pi'][s, t] > reduced_st_path_cost + eps  and len(st_path_nodes) > 0:
                
                currentNumVar = len(self.data['var']['pathVars'][s, t].keys())

                self.data['paths'][s, t][currentNumVar] = { 'arcs' : st_path_arcs, 
                                                            'nodes' : st_path_nodes,
                                                            'duration' : st_path_duration}
                

                newpathVar = self.model.addVar(vtype='C', name=f'Path_{s}_{t}_{currentNumVar}', lb=0, 
                                            obj= self.data['ODpairs'][s, t] * st_path_duration, 
                                            pricedVar = True)

                self.data['var']['pathVars'][s, t][currentNumVar] = newpathVar

                self.model.addConsCoeff(self.data['cons']['passengerflowCons'][s, t], newpathVar, 1.0)

                for arc in st_path_arcs:
                    if (arc, s, t) in self.data['cons']['couplingCons'].keys():
                        c = self.data['cons']['couplingCons'][arc, s, t]
                        self.model.addConsCoeff(c, newpathVar, -1.0)
   
        return {'result':SCIP_RESULT.SUCCESS}

    # The initialisation function for the variable pricer to retrieve the transformed constraints of the problem
    def pricerinit(self):
        for cons_type in ['partitionCons', 'passengerflowCons', 'couplingCons']:
            for i, c in self.data['cons'][cons_type].items():
                self.data['cons'][cons_type][i] = self.model.getTransformedCons(c)


def cxttp_cg(T, L, C, G, A, OD_arcs, ODpairs, vtype='C', output_y=False):
    model = Model("cXPESP_pr")

    model.setPresolve(0)
    model.hideOutput()

    # creating a pricer
    pricer = CyclePricer()
    model.includePricer(pricer, "CyclePricer", "Pricer to identify new cycles in line")

    # adding the initial variables
    cycles_subset = one_cycle_per_line(C, L)
    
    cycleVars = defaultdict(list)
    cycles = defaultdict(list)

    for line in L.keys():
        cycleVars[line].append(model.addVar(vtype=vtype, name=f'Cycle_{line}_{0}', lb=0, obj= cycles_subset[line][0]['duration'] * cycles_subset[line][0]['weight']))
        cycles[line].append(cycles_subset[line][0])


    pathVars = defaultdict(dict)
    paths = defaultdict(dict)

    # construct reduced graph to find path variables for a feasible restricted master problem
    G_red = nx.DiGraph()
    for line in L.keys():
        for cycle in cycles[line]:
            G_red.add_edges_from(cycle['edges'])
    G_red.add_edges_from(A.edges())
    G_red.add_edges_from(OD_arcs.edges())

    costs = nx.get_edge_attributes(G, 'duration')
    for s, t in ODpairs.keys():
        shortest_st_path = nx.shortest_path(G_red, source=(s,'dep'), target=(t, 'arr'), weight="duration", method="bellman-ford")
        st_path_arcs = [(shortest_st_path[i], shortest_st_path[i+1]) for i in range(len(shortest_st_path) - 1)]
        st_path_duration = 0
        for arc in st_path_arcs:
            st_path_duration += costs[arc]
        paths[s, t][0] = {  'nodes' : shortest_st_path, 
                            'arcs': st_path_arcs,
                            'duration': st_path_duration}
        pathVars[s, t][0] = model.addVar(vtype='C', name=f'Path_{s}_{t}_0', lb=0, obj = ODpairs[s, t] * paths[s, t][0]['duration'])

    
    # avoid empty rows
    slack = model.addVar(vtype='C', name=f'slack', lb=0, ub=0, obj= 0) 


    # adding constraints
    partitionCons = dict()
    passengerflowCons = dict()

    for line in L.keys():
        partitionCons[line] = model.addCons(cycleVars[line][0] == 1, 
            name=f'Partition_{line}', separate = False, modifiable = True)
    for s, t in ODpairs.keys():
        passengerflowCons[s, t] = model.addCons(pathVars[s, t][0]  == 1, 
            name=f'PassengerFlow_{s}_{t}', separate = False, modifiable = True)
        
    
    couplingCons = dict()

    for line in C.keys():
        num = 0     
        for arc in C[line].edges():
            for s, t in ODpairs.keys():
                num +=1
                couplingCons[arc, s, t] = model.addCons(slack + quicksum(cycleVars[line][c] for c in [0] if arc in cycles[line][0]['edges']) 
                    - quicksum(pathVars[s, t][p] for p in [0] if arc in paths[s, t][p]['arcs']) >= 0, 
                    name=f'couplingCons_{num}_{s}_{t}', separate = False, modifiable = True)


    # Setting the pricer_data for use in the init and redcost functions
    pricer.data = dict()
    pricer.data['var'] = {'cycleVars' : cycleVars, 'pathVars' : pathVars }
    pricer.data['cons'] = { 'partitionCons' : partitionCons, 
                            'passengerflowCons' : passengerflowCons, 
                            'couplingCons' : couplingCons }
    pricer.data['graph'] = {'C' : C, 'G' : G}
    pricer.data['cycles'] = cycles
    pricer.data['paths'] = paths
    pricer.data['lines'] = L
    pricer.data['T'] = T
    pricer.data['ODpairs'] = ODpairs
    pricer.data['vtype'] = vtype

    model.optimize()
    #model.printBestSol()

    n_cycles = sum([len(pricer.data['var']['cycleVars'][line]) for line in L.keys()] )
    n_paths = sum([len(list(pricer.data['var']['pathVars'][s, t].keys())) for s, t in ODpairs])
    n_cons = sum([len(pricer.data['cons'][cons_type]) for cons_type in pricer.data['cons'].keys()])

    
    if output_y:
        y = dict()   
        sol = model.getBestSol()
        for (s,t) in ODpairs.keys():
            y[s,t] = list()
            for i in pathVars[s, t].keys():
                y[s,t].append((sol[pathVars[s, t][i]], paths[s, t][i]['arcs']))
        return model.getObjVal(), n_cycles, n_paths, n_cons, y
    else:           
        return model.getObjVal(), n_cycles, n_paths, n_cons


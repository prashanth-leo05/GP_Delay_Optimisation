from collections import deque
import os
import networkx as nx
import re
import glob
import gpkit as gp
import sys
import tqdm
#import matplotlib.pyplot as plt

#sys.setrecursionlimit(10000)

class Netlist:
    def __init__(self,file=None,parse=True):
        self.file=file
        self.DAG=nx.DiGraph()
        if parse:
          self.parse()
          
    def __repr__(self):
        return f"{self.DAG}"
    
    def gates_adjacency(self):
        adj_list=dict()
        for node in self.DAG:
            if 'type' in self.DAG.nodes[node].keys():
                adj_list[gate:=node]=[]
                for net in self.DAG.successors(gate):
                  adj_list[gate].extend(self.DAG.successors(net))
        return adj_list

    def longest_path(self):
        ''' Find the path with the largest depth and 
        skip every other node starting from second 
        node as they are not gates '''
        path = nx.dag_longest_path(self.DAG)[1::2]
        path = [f"{self.DAG.nodes[node]['type']}_{node}" for node in path]
        return "=>".join(path)

    def parse(self):
        parser_re = re.compile("""([A-Z]+\d)_  #matches the gate type
                                  (\d+)        #matches the gate index
                                  \s
                                  \(([^\)]*)\) #matches input and output nets of a gate""",
                               re.VERBOSE)  
        with open(self.file,'r') as f:
            text=f.read()

        for match in re.finditer(parser_re,text):
            gate_type, gate_index = match.group(1, 2)
            out, *ins = [net.strip() for net in match.group(3).split(',')]
            for net in ins:
                self.DAG.add_edge(net, gate_index)
            self.DAG.add_edge(gate_index, out)
            self.DAG.nodes[gate_index]['type'] = gate_type
            self.DAG.nodes[gate_index]['size']=gp.Variable(f"x_{gate_index}")

class GPproblem(Netlist):
  min_size=1            #minimum allowable size of each gate
  max_size=100          #maximum allowable size of each gate
  C_L=1000              #load capacitance at all the outputs
  C_max=50              #maximum capacitance at all the inputs
  #logical effort
  g={"NOT": lambda n: 1 , "NAND":lambda n: (int(n)+2)/3, 
     "NOR": lambda n: (2*int(n)+1)/3}

  def __init__(self,file,mode="Timing",model=False,T_spec=None):
    #self.netlist=netlist
    super().__init__(file)
    self.mode=mode
    self.T_spec=T_spec
    if model:
      self.model()

  def __repr__(self):
    return f"GPproblem for {self.netlist}" 

  def model(self):
    constraints=[]
    DAG=self.DAG
    gates_adjacency=self.gates_adjacency()
    nodes_list=DAG.nodes
    T_spec=gp.Variable("T_wall")
    objective=T_spec
    tol=0.000001 # input arrival time tolerance
    for gate in gates_adjacency:
      #maximum and minimum allowable sizes for each gates
      constraints.append(nodes_list[gate]["size"]>=self.min_size)
      constraints.append(nodes_list[gate]["size"]<=self.max_size)
    # different objective if mode is not Timing
    if self.mode=="Area":
      objective=0
      T_spec=self.T_spec
      for gate in gates_adjacency:
        objective += nodes_list[gate]["size"]
    ''' Visit nodes in topologically sorted order so that by the 
    time we visit a node all the arrival times of the node are 
    already calculated.
    '''
    for node in nx.topological_sort(DAG):
      if not DAG.in_degree(PI:=node):
        nodes_list[PI]["arrival_time"]=gp.Variable(f"a_{PI}")
        constraints.append(nodes_list[PI]["arrival_time"] == tol)
        for gate in DAG.successors(PI):
          gate_type=nodes_list[gate]["type"][:-1]
          gate_size=nodes_list[gate]["size"]
          n_inputs=nodes_list[gate]["type"][-1]
          constraints.append(self.g[gate_type](n_inputs)*gate_size <= 
                             self.C_max)

      elif (gate:=node) in gates_adjacency:
        fanout_net=list(DAG.successors(gate))[0]
        gate_type=nodes_list[gate]["type"][:-1]
        n_inputs= nodes_list[gate]["type"][-1]
        gate_size=nodes_list[gate]["size"]
        nodes_list[fanout_net]["arrival_time"]=gp.Variable(f"a_{fanout_net}")
        parasitic_effort = int(n_inputs)
        stage_effort=0
        for fanout_gate in DAG.successors(fanout_net):
          fanout_gate_size=nodes_list[fanout_gate]["size"]
          fanout_gate_type=nodes_list[fanout_gate]["type"][:-1]
          n_inputs=nodes_list[fanout_gate]["type"][-1]
          stage_effort+=fanout_gate_size*self.g[fanout_gate_type](n_inputs)
        if not stage_effort:
          stage_effort+= self.C_L
          constraints.append(nodes_list[PO:=fanout_net]["arrival_time"]<=
                             T_spec)
        stage_effort /= gate_size
        d=stage_effort + parasitic_effort

        for fanin_net in DAG.predecessors(gate):
          constraints.append(d+nodes_list[fanin_net]["arrival_time"]<=
                             nodes_list[fanout_net]["arrival_time"])
      else:
        continue
    return  objective,constraints

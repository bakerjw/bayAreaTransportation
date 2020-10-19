#Author: Mahalia Miller
#Date: Jan. 21, 2013

#this code does iterative travel assignment. 
import sys, util, pdb, math
import networkx as nx

iteration_vals = [0.4, 0.3, 0.2, 0.1] #assign od vals in this amount per iteration. These are recommeded values from the Nature paper, http://www.nature.com/srep/2012/121220/srep01001/pdf/srep01001.pdf

class ITA:
  def __init__(self, G, demand):
    #G is a networkx graph, it can be damaged. Edges need to have a free flow travel time, a capacity, a variable called flow (which we'll change to keep track of flows assigned to each link), and a variable called t_a (which we'll change based on the flows)
    self.G = G
    self.demand = demand
    self.total_demand = 0
    for o in demand.keys():
        for d in demand[o].keys():
            self.total_demand += demand[o][d]

  def assign(self):
    #does 4 iterations in which it assigns od, updates t_a. find paths to minimize travel time and we record route for each od pair
    total_flow = 0
    total_dem = 0
    for i in range(4): #do 4 iterations
      for origin in self.demand.keys(): #origin is an actual A (one end of an actual edge in the graph)
        #find the shortest paths from this origin to each destination
        # print origin
        paths_dict = nx.single_source_dijkstra_path(self.G, origin, cutoff = None, weight = 't_a') #Compute shortest path between source and all other reachable nodes for a weighted graph. Returns dict keyed by by target with the value being a list of node ids of the shortest path
        #print(paths_dict)
        for destination in self.demand[origin].keys(): #actual A or B (one end of an actual edge in the graph)
          # print destination
          od_flow = iteration_vals[i] * self.demand[origin][destination] #to get morning flows, take 5.3% of daily driver values. 11.5/(4.5*6+11.5*10+14*4+4.5*4) from Figure S10 of http://www.nature.com/srep/2012/121220/srep01001/extref/srep01001-s1.pdf
          #get path
          path_list = paths_dict[destination] #list of nodes
        
          trip_completed = True
          total_dem += od_flow
        
          #increment flow on the paths and update t_a
          for index in range(0, len(path_list) - 1):
            u = path_list[index]
            v = path_list[index + 1]
            
            
            
            if self.G.is_multigraph():
              num_multi_edges =  len( self.G[u][v]) #if not multigraph, this just returns the number of edge attributes
              if num_multi_edges >1: #multi-edge
                #identify multi edge with lowest t_a
                best = 0
                best_t_a = float('inf')
                # print self.G[u][v].keys()
                for multi_edge in self.G[u][v].keys():
                  new_t_a = self.G[u][v][multi_edge]['t_a'] #causes problems
                  if (new_t_a < best_t_a) and (self.G[u][v][multi_edge]['capacity']>0):
                    best = multi_edge
                    best_t_a = new_t_a
              else:
                best = 0
              if (self.G[u][v][best]['capacity']>0):
                self.G[u][v][best]['flow'] += od_flow
                #total_flow += od_flow
                t = util.TravelTime(self.G[u][v][best]['t_0'], self.G[u][v][best]['capacity'])
                travel_time= t.get_new_travel_time(od_flow) #TODO #min(t.get_new_travel_time(od_flow), self.G[u][v][best]['distance_0']*1.0/3600.0) #distance in miles, t_a in seconds!! So we are saying that the minimum of the t_a and distance (in miles) * (1 hr/ 1 mile) * (1hr / 3600s)
                self.G[u][v][best]['t_a'] = travel_time
              else:
                trip_completed = False
            else:
              try:
                if (self.G[u][v]['capacity']>0):
                  self.G[u][v]['flow'] += od_flow
                  #total_flow += od_flow
                  t = util.TravelTime(self.G[u][v]['t_0'], self.G[u][v]['capacity'])
                  travel_time= t.get_new_travel_time(od_flow) #TODO #min(t.get_new_travel_time(od_flow), self.G[u][v][best]['distance_0']*1.0/3600.0) #distance in miles, t_a in seconds!! So we are saying that the minimum of the t_a and distance (in miles) * (1 hr/ 1 mile) * (1hr / 3600s)
                  self.G[u][v]['t_a'] = travel_time #in seconds
                else:
                  trip_completed = False
              except KeyError as e:
                print('found key error: ', e)
                pdb.set_trace()
          total_flow += od_flow * trip_completed
    #print 'in ita:', total_flow, total_dem
    return self.G, total_flow, total_dem

      





def test():
  #create graph info
  G = nx.MultiDiGraph()
  G.add_node(1)
  G.add_node(2)
  G.add_node(3)
  G.add_edge(1,2,capacity_0=1000,capacity=1000,t_0=15,t_a=15,flow=0, distance=10)
  G.add_edge(1,2,capacity_0=3000,capacity=3000,t_0=20,t_a=20,flow=0, distance=10)
  G.add_edge(2,3, capacity_0=0, capacity=0, t_0 = 10, t_a = 10, flow = 0, distance = 10)
  #get od info. This is in format of a dict keyed by od, like demand[sd1][sd2] = 200000.
  demand = {}
  demand[1] = {}
  demand[1][2] = 8000 #divide by 0.053 since that is what we multiply by above
  demand[2] = {}
  demand[2][3] = 4000

  #call ita
  it = ITA(G,demand)
  newG, total_flow, total_demand = it.assign()
  print(newG)
  print('total flow:', total_flow)
  print('total demand:', total_demand)
  print('proportion of demand met:', (total_flow/total_demand))
  for n,nbrsdict in newG.adjacency_iter():
    for nbr,keydict in nbrsdict.items():
      for key,eattr in keydict.items():
        print (n, nbr, eattr['flow'])
  print('should have flow of 3200 and 4800')

if __name__ == '__main__':
  test()

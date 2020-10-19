#Author: Mahalia Miller
#Date: Jan. 21, 2013

import pdb

def build_demand(trip_filename, centroid_filename):
  '''demand dict has keys of actual nodes (one per travel district)'''
  demand_dict = {}
  sd_dict = {}
  with open(centroid_filename,'rb') as f:
    read_data = f.read().splitlines()
    for row in read_data[1:]:
      tokens = row.split(',')
      sd_dict[str(int(tokens[0]))] =  str(int(tokens[2]))
      demand_dict[str(int(tokens[2]))] = {} #demand_dict[A] is important for line 19 below where we assign travel from origin to destination

  with open(trip_filename,'rb') as f:
    read_data = f.read().splitlines()
    # for row in read_data[4:]:
    for row in read_data: #modified feb 10, 2014
      tokens = row.split(',')
#      print tokens
      #grab number of trips where someone is driving only (to approximate the number of vehicles on the road)
      demand_dict[str(sd_dict[str(int(tokens[0]))])][str(sd_dict[str(int(tokens[1]))])] = int(tokens[2]) # + int(tokens[3])  #int(tokens[12]) Note: tokens[3] has passengers but we just want drivers since we want cars.
      if (int(tokens[2]) + int(tokens[3])) > int(tokens[12]):
        print 'what is going on?'
  return demand_dict


if __name__=='__main__':
  demand_dict = build_demand('input/BATS2000_34SuperD_TripTableData.csv', 'input/superdistricts_centroids_dummies.csv')
  pdb.set_trace()
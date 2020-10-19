Created by Mahalia Miller
July 2, 2014

This folder contains example code and data to illustrate the efficient transportation model using iterative traffic assignment described in Chapter 2 of:

M. Miller, “Seismic risk assessment of complex transportation networks,” PhD Thesis, Stanford University, 2014.

This chapter contains citation to the relevant prior work.

Running this code uses the Python programming language and some Python packages. If you are an academic or if your company will buy the license, you can save effort of installing all the different packages by downloading a python bundle from: https://www.enthought.com/products/canopy/academic/
If you are new to Python, you can see the official documentation and help here: https://www.python.org/docs/

****
You can run the model with the example case study data by typing into Terminal (on a Mac), into your Python command prompt, or other suitable command prompt:
python mahmodel_road_only.py
*****

The following files are included in this folder:
bd.py -- a function for building the travel demand
bridge_metadata_NBI.xlsx -- a file that has the background data about the case study road bridges
input/20140114_master_bridge_dict.pkl -- sample data for the SF Bay Area road components (bridges)
input/20140114_master_transit_dict.pkl -- sample data for the SF Bay Area BART components
input/BATS2000_34SuperD_TripTableData.csv -- average daily trips between different superdistricts. See http://analytics.mtc.ca.gov/foswiki/Main/DataDictionary for more info
input/graphMTC_CentroidsLength3int.gpickle -- the graph of the SF Bay Area highways and key local roads
input/sample_ground_motion_intensity_map_JUST_THREE.txt -- ground-motion intensity map data for just three ground-motion intensity maps. The columns refer to:  first column is simulation number, second is fault id, third is magnitude, fourth is the annual occurrence rate (SUPER USEFUL), fifth is Sa (NOT logSa) in site new ID 1, sixth is Sa in site new ID 2, ...site ID n
input/sample_ground_motion_intensity_maps_road_only_filtered.txt -- same columns as the previous file but this has a full hazard-consistent set of events
input/superdistricts_centroids_dummies.csv -- file that has a centroidal/dummy link node for each superdistrict (for traffic assignment)
input/superdistricts_clean.csv -- file that has a few nodes in each superdistrict (for traffic assignment)
ita.py -- the core function that does the iterative traffic assignment
mahmodel_road_only.py -- the main file with only road damage considered
mahmodel.py -- alternative main file that also keeps track of which transit components are damaged
make_bridge_dict.py -- a sample file for showing how to create your own master_bridge_dict.pkl
output -- a folder foroutput
README_quick_traffic_model.txt -- documentation for this folder
transit_to_damage.py -- a file that gives some helper functions for translating damaged components to nonoperational transit lines for the case study
util.py -- helper functions


Note: you may want to know what will happen to transit. In that case, you can run the mahmodel.py file. The saved output will have bridge ids. The new ids bigger than 1743 are damaged BART structures. See transit_to_damage.py for info about how the new ids correspond to BART, VTA, MUNI, and Caltrain.

Additional note: if you want to change the fragility data, you will want to change the line 83 of mahmodel_road_only.py to use a different key of the master_bridge_dict.pkl, which is a Python dictionary with the fragility data. To make that work, you'l have to add a new key and value for each bridge with your new fragility data to master_brdige_dict.pkl and save this updated version. See make_bridge_dict.py for ideas about how to add key, value pairs to Python dictionaries. You can also email me at mahaliakmiller@gmail.com for support. :)

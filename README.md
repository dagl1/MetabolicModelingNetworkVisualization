# Metabolic network visualization
![image](https://github.com/dagl1/MetabolicModelingNetworkVisualization/assets/24440380/8b9d14c1-7626-4796-9809-e6ffb53342f0)


## What?
This tool can load ESHER (https://escher.github.io/#/) JSON maps and ESHER-ready JSON Models created using COBRA (https://opencobra.github.io/cobratoolbox/stable/index.html) to allow for easy semi-automatic visiualization and editing of (small) network maps, which subsequently can be saved as ESHER-ready maps and further explored within ESCHER.

## Why? 
At the department where I contributed to a project, we have several methods of creating submodels of larger metabolic networks. These submodels need to be explored to validate and verify their usefulness, normally done in ESCHER. This method is tedious and requires substantial time investment, therefore I thought it would be an interesting project to create an initial visualization tool in which users can (approximately) create the required network with only a few clicks of the mouse. Once acceptable, this can then be further loaded into ESCHER (as it is a lot prettier).

![ezgif com-optimize](https://github.com/dagl1/MetabolicModelingNetworkVisualization/assets/24440380/dfa014f6-63ff-481b-8fe8-b0bd35c150de)
Example of using ESCHER to manually add each reaction.

## Why not use some other tools? 
One of the primary reasons for this project was to develop something which fit directly into our current workflow; most students were not using any other tools out there, possibly due to not knowing they existed, or because the usage wasn't integrated enough with their workflow. 
Other tools that look promising (and much more developed) include: 
(https://bioinformatics.mdanderson.org/public-software/sammi/) 
(https://cave.biodesign.ac.cn/)
(https://metexplore.toulouse.inrae.fr/metexploreViz/doc/index.php)

## Installation
This application was created and tested using Python 3.11 on Windows 10 Pro, but most likely will work with most Python 3 installations. In the future I will most likely package the application as an .exe and/or compile in CPython for additional performance, in those cases users on non-Windows operating systems might be required to utilize the standard python version, or compile their own binaries.
Required packages:
json, pygame, numpy, pandas, tkinter, threading, queue, copy, re, and openpyxl (does not come standard with pandas but is required for handling the excel files).

## Running
To run the application after installing all required packages, run main.py through the console using:
``` cmd
cd C:/location/of/your/main/file
python main.py 
```
A window should open where a JSON file can be uploaded, this JSON file should have been created from a COBRA model using the "save_json_model()" function:
``` python
    # Load the MATLAB model
    model = load_matlab_model(os.path.join(new_models_directory, file))
    
    # Create the JSON file in the "ModelsJSON" directory
    json_file = name + '.json'
    json_path = os.path.join(json_files_directory, json_file)
    save_json_model(model, json_path)
```
Confirming the file by uploading it to the application will open the settings screen. 
In case you already have a pre-made JSON map you can upload this by clicking the "upload previously made JSON model" button. 
(the JSON maps are not to be confused with the aforementioned JSON model; a map is the positions of the particular nodes, either created through this app or taken from ESCHER using the Map > Save map JSON option.)
![image](https://github.com/dagl1/MetabolicModelingNetworkVisualization/assets/24440380/61a832e5-69b5-41ab-b4c4-e9e57de79206)

In addition one can upload an Excel by pressing "Browse excel file" and subsequently "upload excel file". If you are running the program through an IDE, a text prompt will indicate if the excel file was accepted or not. This file should have no headers and contain the flux per reation in column 2, expression data per reaction in column 3, reaction IDs in MARxxxx notaion in column 4, and regular reaction names in column 5:
![image](https://github.com/dagl1/MetabolicModelingNetworkVisualization/assets/24440380/09d280e5-396b-4115-aa57-d372b794d561)

Running a pre-made map uploaded to the program can be ran by clicking the "Restart the latest simulation" button. 
To start a new simulation using just the JSON model and the current settings, press the "Start simulation" button.

## Settings
The settings of this application update the simulation real time, it is therefore not necessary to set particular settings before running the simulation. However it is of course the case your final outcome may differ if you change settings during a run versus starting with those settings. In those cases, click the "Stop simulation" button. To start the latest run, press the "Restart the latest simulation" button, while for running a new simulation with the current settings, click the "Start simulation" button again.

### Keyboard and mouse controls:
```
WASD or arrow keys to move around the canvas
E to zoom in
Q to zoom out
Left-mouse button (selected node + hold) to drag node. Can also be used to drag multiple nodes when they are selected by selection rectangle
Left-mouse button (no node selected + hold) to create selection rectangle
Shift + left-mouse button (no node selected + hold, pressing shift at release) to create selection rectangle that marks/unmarks nodes as "special"
Right-mouse button to unselect nodes selected through rectangle
Shift + Left-mouse button (click node) to mark/unmark node as "special"
X + Left-Mouse button to split nodes based on marked "special" nodes
F + Left-Mouse button to fix/unfix node
V + Left-mouse button OR Middle-mouse button to make/unmake nodes visible (invisible nodes can be added back through the visibility list (see "Included/not included and visible/not visible" below)
G to enter search mode
(In search mode) Escape to exit search mode
(In search mode) Tab to select suggestion
(In search mode) Up/down arrow keys to change suggestion
(In search mode) Enter, mark/unmark selected node for finding purposes
```

### Attraction and Repulsion sliders:
```
Normal repulsion indicates the amount of repulsion each node receives from each other node in the simulation.

Intercompartment repulsion indicates the amount of repulsion the special "compartment" nodes feel from each other;
    does not directly affect any metabolite or reaction nodes.

Boundary Repulsion indicates the strength the boundary nodes (invisible set of nodes surrounding the canvas) excert
    on the rest of the nodes, useful for centering the nodes slightly (might make compartment nodes "clip out", in
    which case you should fix them in place or lower the Boundary repulsion and dragging the compartment nodes back into the canvas).

Input_metabolite-Input_metabolite Repulsion is used to spread the input and output metabolites.
    These can only move on their own on the X-axis, you can still drag them yourself if you don't want to rely on automatic organization.

Reaction-Metabolite Attraction is the strength attraction that nodes connected by lines/edges will feel to each other.

Metabolite-Compartrment Attraction is the strenght the metabolites of a particular compartment will feel with their connected compartment node. Is mostly useful for separating clusters of reactions part of a single compartment.
```
![image](https://github.com/dagl1/MetabolicModelingNetworkVisualization/assets/24440380/37a26ed8-d4a3-4738-b284-b783213f681b)

### Included/not included and visible/not visible
There are 2 lists per catagory, the included/not included catagory designates which nodes are part of the simulation. Not included nodes can be shown on screen, but will excert no attraction or repulsion forces. Metabolites can be automatically excluded based on the amount of connections they have using the slider beneath the lists. An alternative is to split highly connected metabolites (see "Splitting and Merging" below). When starting a simulation, excluded nodes will be invisible, however when changing inclusion criteria, visibility doesn't automatically change. To update the visibility to match the excluded/included nodes, press the "Update visibility based on included metabolites" button.
One can manually update the lists by double clicking 

![image](https://github.com/dagl1/MetabolicModelingNetworkVisualization/assets/24440380/08602957-3ad2-4356-8e13-5581cda20f4b)


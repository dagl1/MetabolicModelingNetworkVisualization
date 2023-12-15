# Metabolic network visualization
![image](https://github.com/dagl1/MetabolicModelingNetworkVisualization/assets/24440380/8b9d14c1-7626-4796-9809-e6ffb53342f0)


## What?
This tool can load ESHER (https://escher.github.io/#/) JSON maps and ESHER-ready JSON Models created using COBRA (https://opencobra.github.io/cobratoolbox/stable/index.html), specifically originating from a task-approach such as CellFie (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8577426/), to allow for easy semi-automatic visiualization and editing of (small) network maps, which subsequently can be saved as ESHER-ready maps and further explored within ESCHER.

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

Alternatively you can get the exact versions used in this project by using the requirements.txt with the following pip install command in the terminal:
```
pip install -r requirements.txt
```
Warning! It is better to create a virtual environment for your python deployment as to not interfere with any other projects you have going on:
```
cd C:/chosen_file_location/
python -m venv myenv 
.\myenv\Scripts\activate # for Windows
# or
source myenv/bin/activate # for Unix/Mac
pip install -r requirements.txt
```
For both of these approaches, you will need to have the requirements.txt file in the folder where you are running the terminal from!

## Running
To run the application after installing all required packages, run main.py through the console using:
``` cmd
cd C:/location_of_your_file/
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
**Alternatively, you can run the file directly by providing the file's location (either absolute path, or filename if it is in the correct folder).**
``` cmd
python main.py json_file_location/json_file.json
```

In case you already have a pre-made JSON map you can upload this by clicking the "upload previously made JSON model" button. 
(the JSON maps are not to be confused with the aforementioned JSON model; a map is the positions of the particular nodes, either created through this app or taken from ESCHER using the Map > Save map JSON option.)
![image](https://github.com/dagl1/MetabolicModelingNetworkVisualization/assets/24440380/61a832e5-69b5-41ab-b4c4-e9e57de79206)

In addition one can upload an Excel by pressing "Browse excel file" and subsequently "upload excel file". If you are running the program through an IDE, a text prompt will indicate if the excel file was accepted or not. This file should have no headers and contain the flux per reation in column 2, expression data per reaction in column 3, reaction IDs in MARxxxx notaion in column 4, and regular reaction names in column 5:
![image](https://github.com/dagl1/MetabolicModelingNetworkVisualization/assets/24440380/09d280e5-396b-4115-aa57-d372b794d561)

Running a pre-made map uploaded to the program can be ran by clicking the "Restart the latest simulation" button, note that all nodes will start frozen, press the "Unfreeze all nodes" button to allow for movement of the nodes.
To start a new simulation using just the JSON model and the current settings, press the "Start simulation" button.

**As of right now, this program is only fully functional for an irreversible COBRA model, reversible models will not show directionality correctly, have too many nodes, and might not save correctly too ESCHER** 

## Settings
The settings of this application update the simulation real time, it is therefore not necessary to set particular settings before running the simulation. However it is of course the case your final outcome may differ if you change settings during a run versus starting with those settings. In those cases, click the "Stop simulation" button. To start the latest run, press the "Restart the latest simulation" button, while for running a new simulation with the current settings, click the "Start simulation" button again.

### Keyboard and mouse controls:
```
Escape to stop application when in menu
F1 to stop simulation
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
Ctrl + Z will undo the last movement action
Ctrl + Y will redo the last undone movement action
Ctrl + F OR Shift + F can unfix/fix all selected nodes
Ctrl + 1-9 Will designate currently selected nodes as part of a control-group (designated by the number)
1-9 Select nodes part of the control group of the corresponding number
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

### Manually fixing nodes, changing visibility, or assigning a "special" mark"
By clicking a node to move it, you will automatically freeze the node, these nodes are "fixed in place". You can manually fix/unfix such nodes by holding f + Left-mouse click. You can also change visibility of nodes in such a way by holding v + Left-mouse click on a node, this will change the visibility of the node. An invisible node will still be part of the simulation, but it, and connections to it, will be invisible. If you cannot find a particular node's location, you can always use the visible/not visible list (see above), to turn the node visible again. Last, to mark nodes as "special" for visualization purposes but also used for splitting nodes (see below), you can shift + Left_mouse click them. Clicking again will subsequently unmark them.

When selecting a group of nodes and subsequently moving them, these nodes will all become fixed. You can unfix them by pressing Ctrl + f or Shift + f, which will fix/unfix the selected nodes. In cases where some of the selected nodes are fixed, while others are unfixed, the program will fix/unfix all of them to one status (they will all beocme fixed/unfixed depending on the first node's status).

### Merging/Splitting nodes
Sometimes some metabolites or reactions connect two different clusters of nodes and you might want to split these so that this one node does not pull both clusters to each other (and lead to a line stretching the entire screen). In this case you can manually seperate (split) the node by first marking some of its connections as "special" by Shift + Left-mouse button clicking (or holding Shift while drag selecting), then splitting the node with x + Left-mouse button. Split nodes will not be connected to each other, but can be merged by dragging one over one of its partners. Nodes can be split as many times as necessary and will remember all their partner nodes.
![ezgif com-optimize ESCHER DEMO Merging Splitting ](https://github.com/dagl1/MetabolicModelingNetworkVisualization/assets/24440380/74c74dfd-6ec9-498b-b8cd-dfaf81d56bcb)

### Search mode
You might want to find/highlight a node without finding it in the list of visible/invisible nodes. For that you can press g and enter search mode. While in search mode, most controls (other than clicks from the mouse), will not work and instead you can type into the search box on the bottom right (you do not need to press the button). To exit search mode, press Escape. While searching you will see recommendations, to autocomplete a recommendation press tab. To select a node press enter. To change which recommnendation is selected use the up/down arrow keys (selected recommendation is highlighted). The text in the yellow boxes will show where approximately in the names of the suggestions your input exists.

### Zooming in and out
You might feel like you cannot see exactly what is going on somewhere, in those cases zooming in (e) or out (q), and moving (w, a, s, d keys or arrow keys) the canvas might help. Zooming in is only visual, so the physics calculations will remain exactly the same.

### Undo and redo
Just like most text editing software, you can use Ctrl + Z to undo movement actions, and Ctrl + Y to redo undone movement actions. If you undo, then perform a new move action, your action still stored in the "redo" actions will be lost.

### Canvas settings
Some users reported that the canvas was too large for their screen, so now you manually set the canvas of the simulation in the options menu. To lock-in your new canvas width and height, click the "Set Dimensions" button.

### Other controls
There are several buttons that can be used to change some aspects of the application.
The show/unshow lines button will change the visibility of the lines connecting the nodes, while the show/unshow (default at unshow) compartments will show the special "compartment" nodes, to which each metabolite of that compartment is connected to. Color/uncolor (default uncolor) will change the color of the lines from white, to the color of the compartment the line is connected to. 
Freeze all nodes and unfreeze all nodes will fix or unfix all nodes, if you have made a nice looking arrangement, be careful wiht the unfree button! Then there are the fontsize and rectangle size sliders, which change the fontsize of the text and the rectangle size respectively. Beneath those there are options for showing reaction and metabolite names, fixed metabolites are input/output metabolites of a specific task using the CellFie (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8577426/) approach.
Fluxes and expression data can only be shown by uploading the designated excel file (see above).
![image](https://github.com/dagl1/MetabolicModelingNetworkVisualization/assets/24440380/5bc2e0c0-6cf0-423b-b2d8-36c0a71734cd)

### Saving the simulation
By clicking the "save latest simulation as .json" button, you can save your map, both for reloading into the application, but also applicable for importing it as a map for ESCHER where it can be further edited.
![image](https://github.com/dagl1/MetabolicModelingNetworkVisualization/assets/24440380/8a916909-148c-4389-ba0f-178605972b35)

#### Author notes:
If there are any issues in the usage or you have some new ideas for killer features, please raise an issue on this Github. 
You are free to use this however and whenever you like, if you use this for any analysis for papers, an acknoledgment or link would be nice but is not necessary! 

Dagl1 / Jelle Bonthuis 

Created: 25-11-2023

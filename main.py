# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
import os, json, pygame, math
import numpy as np
import pandas as pd
import tkinter as tk
from tkinter import filedialog

def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print_hi('PyCharm')

# See PyCharm help at https://www.jetbrains.com/help/pycharm/



class MainLoop:

    def __init__(self):
        self.WIDTH, self.HEIGHT = 1200, 800
        self.FINAL_WIDTH, self.FINAL_HEIGHT = 2400, 1600
        pygame.init()
        self.canvas = pygame.Surface((self.WIDTH,self.HEIGHT))
        self.screen = pygame.display.set_mode((self.FINAL_WIDTH,self.FINAL_HEIGHT))

class TkinterSetup:
    pass
class JSONViewer:
    def __init__(self, root):
        self.root = root
        self.root.title("JSON Viewer")

        # Create a StringVar to store the selected file path
        self.entry_var = tk.StringVar()

        # Create and place widgets
        label = tk.Label(self.root, text="Selected File:")
        label.pack(pady=10)

        entry = tk.Entry(self.root, textvariable=self.entry_var, state="readonly", width=40)
        entry.pack(pady=10)

        browse_button = tk.Button(self.root, text="Browse", command=self.browse_file)
        browse_button.pack(pady=10)

        upload_button = tk.Button(self.root, text="Upload", command=self.upload_file)
        upload_button.pack(pady=10)

    def browse_file(self):
        file_path = filedialog.askopenfilename(filetypes=[("JSON files", "*.json")])
        self.entry_var.set(file_path)

    def upload_file(self):
        file_path = self.entry_var.get()
        try:
            with open(file_path, 'r') as file:
                # Read JSON data
                json_data = json.load(file)

                # Open a new window for JSON data visualization
                self.open_json_viewer(json_data)
        except FileNotFoundError:
            print("File not found.")

    def open_json_viewer(self, json_data):
        viewer_window = tk.Toplevel(self.root)
        viewer_window.title("JSON Data Viewer")

        # Extract names from JSON data
        names = [entry.get("name", "") for entry in json_data]

        # Create a variable to store selected indices
        selected_indices = tk.StringVar(value="")

        # Create a listbox with checkbuttons
        listbox = tk.Listbox(viewer_window, selectmode=tk.MULTIPLE, listvariable=selected_indices, height=len(names))
        for i, name in enumerate(names):
            listbox.insert(tk.END, name)
        listbox.pack(pady=10)

        # Create a slider
        slider = tk.Scale(viewer_window, from_=0, to=len(names) - 1, orient=tk.HORIZONTAL, label="Select a range")
        slider.pack(pady=10)

        # Create a button to perform some action based on the selected items
        action_button = tk.Button(viewer_window, text="Perform Action", command=lambda: self.perform_action(listbox.curselection(), slider.get()))
        action_button.pack(pady=10)

    def perform_action(self, selected_indices, slider_value):
        print("Selected Indices:", selected_indices)
        print("Slider Value:", slider_value)
# Create the main Tkinter window
root = tk.Tk()

# Create an instance of the JSONViewer class
# json_viewer = JSONViewer(root)

# Start the Tkinter event loop
# root.mainloop()

filepath ="E:\Git\MaCSBio-GEM\General\Functions\Metabolic_tasks\Task322Sample1Alpha0.50Beta0.00_Date2023-10-25_21-51-37\ModelsJSON/newModelMILP_OR_AdjustTask322Sample1Scored2.152average4.105alpha0.500beta0.000.json"

with open(filepath, 'r') as file:
    # Read JSON data
    json_data = json.load(file)

# print(json_data)
class MetaboliteNodes:
    def __init__(self, id, name, compartment):
        self.id = id
        self.name = name
        self.compartment = compartment
        self.x, self.y = 0,0
        self.amount_of_reactions = 0
        self.reactions = []
        self.visible = True
        self.fixed_node_input = False
        self.fixed_node_output = False
class ReactionEdges:
    def __init__(self, id, name, metabolites, lb, ub):
        self.id = id
        self.name = name
        self.metabolites = metabolites
        self.lb = lb
        self.ub = ub
        self.visible = True #for compartment



class MetabolicNetwork:
    def __init__(self, nodes, edges):
        self.metabolites_nodes = nodes
        self.reaction_edges = edges
        self.set_amount_of_reactions_per_node()
        self.set_slider()
        self.filter_input_output_reactions()
        self.set_reactions_for_nodes()
        self.commonly_excluded_names = ["H+", "H20", "NADH", "NAD+", "NADP+", "NADPH", "FAD", "FADH", "CO2"]
        self.included_metabolites = []
        self.not_included_metabolites = []
        self.create_lookup_dicts()
        self.update_included_not_included_based_on_slider()


    def filter_input_output_reactions(self):
        reaction_edges_copy = []
        for reaction in self.reaction_edges:
            if not (reaction.lb >0 and reaction.ub <1000) or not reaction.name.startswith('temporary_exchange_') :
                reaction_edges_copy.append(reaction)
            else:
                self.set_fixed_metabolite(reaction)
        self.reaction_edges = reaction_edges_copy

    def set_fixed_metabolite(self, reaction):
        mets = reaction.metabolites
        for metabolite, stoichiometry in mets.items():
            for met in self.metabolites_nodes:
                if metabolite == met.id:
                    if stoichiometry >0:
                        met.fixed_node_input = True
                    if stoichiometry <0:
                        met.fixed_node_output = True
                    break

    def find_reactions_that_should_generally_be_excluded(self):
        for commonly_excluded_name in self.commonly_excluded_names:
            for metabolite in self.metabolites_nodes:
                if metabolite.name.startswith(commonly_excluded_name):
                    self.not_included_metabolites.append(metabolite)
        # self.included_metabolites = [item for item in self.metabolites_nodes if item not in self.not_included_metabolites]

    def create_lookup_dicts(self):
        self.dict_metabolite_names = {}
        self.dict_metabolite_alternative_names = {}
        for metabolite in self.metabolites_nodes:

            self.dict_metabolite_names[metabolite.id] = metabolite
            self.dict_metabolite_alternative_names[metabolite.name] = metabolite
        self.dict_reactions = {}
        for reaction in self.reaction_edges:
            self.dict_reactions[reaction.id] = reaction


    def find_reaction(self, reaction_id):
        reaction = self.dict_reactions[reaction_id]
        return reaction

    def find_metabolite(self, metabolite_id):
        if metabolite_id.startswith("MAM0"): # should be changed if metabolite IDs ever go past 10000
            metabolite = self.dict_metabolite_names[metabolite_id]
        else:
            metabolite = self.dict_metabolite_alternative_names[metabolite_id]
        return metabolite

    def update_included_not_included_based_on_slider(self):
        self.not_included_metabolites = []
        self.included_metabolites = []
        self.find_reactions_that_should_generally_be_excluded()
        for metabolite in self.metabolites_nodes:
            if metabolite.amount_of_reactions >= self.slider_value:
                self.not_included_metabolites.append(metabolite)
        self.included_metabolites = [item for item in self.metabolites_nodes if item not in self.not_included_metabolites]
    def set_reactions_for_nodes(self):
          for reaction in self.reaction_edges:
            if not (reaction.lb >0 and reaction.ub <1000) or not reaction.name.startswith('temporary_exchange_'):
                for key, value in reaction.metabolites.items():
                    for metabolite in self.metabolites_nodes:
                        if key == metabolite.id:
                            metabolite.reactions.append(reaction.id)
                            break

    def set_amount_of_reactions_per_node(self):
        # dont count temporary reactions
        for reaction in self.reaction_edges:
            if not (reaction.lb >0 and reaction.ub <1000) or not reaction.name.startswith('temporary_exchange_'):
                for key, value in reaction.metabolites.items():
                    for metabolite in self.metabolites_nodes:
                        if key == metabolite.id:
                            metabolite.amount_of_reactions += 1
                            break

    def set_slider(self):
        self.slider_start = 1
        highest_value = 0
        total_edges = 0
        for metabolite in self.metabolites_nodes:
            total_edges += metabolite.amount_of_reactions
            if metabolite.amount_of_reactions > highest_value:
                highest_value = metabolite.amount_of_reactions
        self.slider_end = highest_value
        self.slider_value = math.floor(highest_value * 0.8)


list_edges_reactions = []
list_nodes_metabolites = []

for key, value in json_data.items():
    list_in_dic = value# Get the list of metabolites or an empty list if not present
    if key == "metabolites":
        for metabolite in list_in_dic:
            # Access information about each metabolite
            metabolite_id = metabolite['id']
            metabolite_name = metabolite['name']
            metabolite_compartment = metabolite['compartment']
            metabolite_name = f"%s[%s]"%(metabolite_name, metabolite_compartment)
            list_nodes_metabolites.append(MetaboliteNodes(metabolite_id, metabolite_name,metabolite_compartment))
    if key == "reactions":
        for reaction in list_in_dic:
            reaction_id = reaction["id"]
            reaction_name = reaction["name"]
            reaction_metabolites = reaction["metabolites"]
            reaction_lb = reaction["lower_bound"]
            reaction_ub = reaction["upper_bound"]
            list_edges_reactions.append(ReactionEdges(reaction_id, reaction_name, reaction_metabolites, reaction_lb, reaction_ub))


met_network = MetabolicNetwork(list_nodes_metabolites,list_edges_reactions)

for met in met_network.metabolites_nodes:
    print(met.amount_of_reactions,met.id, met.name, met.reactions, met.fixed_node_input, met.fixed_node_output)
for rec in met_network.reaction_edges:
    # print(rec.id)
    pass
print(len(met_network.dict_metabolite_alternative_names))
print(len(met_network.dict_metabolite_names))
print(len(met_network.included_metabolites), len(met_network.not_included_metabolites), met_network.slider_value)
print(met_network.included_metabolites)
print(met_network.dict_metabolite_alternative_names)
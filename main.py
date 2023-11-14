# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
import os, json, pygame, math, random
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
        self.flux = 0
        self.expression = 0


class MetabolicNetwork:
    def __init__(self, nodes, edges):
        self.metabolites_nodes = nodes
        self.reaction_edges = edges
        self.set_amount_of_reactions_per_node()
        self.set_slider()
        self.filter_input_output_reactions()
        self.set_reactions_for_nodes()
        self.commonly_excluded_names = ["H+", "H2O", "NADH", "NAD+", "NADP+", "NADPH", "FAD", "FADH", "CO2"]
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
        if reaction_id in self.dict_reactions:
            reaction_ = self.dict_reactions[reaction_id]
        try:
            return reaction_
        except:
            return None

    def find_metabolite(self, metabolite_id):
        if metabolite_id.startswith("MAM0"): # should be changed if metabolite IDs ever go past 10000
            if metabolite_id in self.dict_metabolite_names:
                metabolite = self.dict_metabolite_names[metabolite_id]
        else:
            if metabolite_id in self.dict_metabolite_alternative_names:
                metabolite = self.dict_metabolite_alternative_names[metabolite_id]
        try:
            return metabolite
        except:
            return None

    def update_included_not_included_based_on_slider(self):
        self.not_included_metabolites = []
        self.included_metabolites = []
        self.find_reactions_that_should_generally_be_excluded()
        for metabolite in self.metabolites_nodes:
            if metabolite.amount_of_reactions >= self.slider_value:
                self.not_included_metabolites.append(metabolite)
        self.included_metabolites = [item for item in self.metabolites_nodes if item not in self.not_included_metabolites]
        for metabolite in self.metabolites_nodes:
            if metabolite.fixed_node_input or metabolite.fixed_node_output:
                self.included_metabolites.append(metabolite)
        self.not_included_metabolites = [item for item in self.metabolites_nodes if item not in self.included_metabolites]
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

    def calculate_network(self):
        self.create_nodes_and_edges_for_network_calculation()
        self.run_network_simulation()
        nodes, edges = self.extract_non_pseudo_nodes_and_edges()
        return nodes, edges

    def create_nodes_and_edges_for_network_calculation(self):
        self.network_nodes = []
        self.network_edges = []
        #create initial nodes and connect edges

        for metabolite in self.included_metabolites:
            self.network_nodes.append(metabolite.id)
            for reaction in metabolite.reactions:
                reaction_ = self.find_reaction(reaction)
                for met, stoichiometry in reaction_.metabolites.items():
                    if metabolite.id == met and stoichiometry > 0:
                        for met, stoichiometry in reaction_.metabolites.items():
                            if stoichiometry > 0 and met != metabolite.id and met in [metabolit.id for metabolit in self.included_metabolites]:
                                met_ = self.find_metabolite(met)
                                self.network_edges.append((metabolite.id, met_.id,reaction_))
        self.create_compartment_pseudo_nodes_and_edges()
    def create_compartment_pseudo_nodes_and_edges(self):
        compartments = []
        for metabolite in self.included_metabolites:
            if metabolite.compartment not in compartments:
                compartments.append(metabolite.compartment)
        self.network_nodes.extend(compartments)
        for compartment in compartments:
            for metabolite in self.included_metabolites:
                if metabolite.compartment == compartment:
                    self.network_edges.append((compartment,metabolite.id, None))
        self.compartments = compartments


    def run_network_simulation(self):
        self.network_node_positions = {}

        self.set_fixded_and_non_fixed_nodes()
        self.provide_random_positions_for_non_fixed_nodes()

        # translate problem to instances of classes of nodes+connected edges to other nodes
        # way to calculate distance between all nodes
        # calculate forces per node
        ## repulsion 1/(distance cubed + 0.0001)
        ## attraction per edge towards edge partner, attraction, distance squared
        # move based on force

        # add in cyclic repulsion for creation of circles



    def set_fixded_and_non_fixed_nodes(self):
        # add code for setting height and width depending on amount of fixed reactions
        # and on max path lenght between metabolites
        output_counter = 1
        input_counter = 1
        random.seed(666)
        for metabolite_id in self.network_nodes:
            metabolite = self.find_metabolite(metabolite_id)
            if metabolite is not None and metabolite.fixed_node_output == True:
                self.network_node_positions[metabolite_id] = (output_counter*100, 800)
                output_counter += 1
            elif  metabolite is not None and metabolite.fixed_node_input == True:
                self.network_node_positions[metabolite_id] = (input_counter*100, 100)
                input_counter += 1
            else:
                self.network_node_positions[metabolite_id] = (random.randint(200, 700 ), random.randint(200, 700))

    def provide_random_positions_for_non_fixed_nodes(self):
        pass

    def extract_non_pseudo_nodes_and_edges(self):
        #remove
        nodes, edges = 0, 0
        #remove
        return nodes, edges
        pass


class Application:
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

        self.browse_button = tk.Button(self.root, text="Browse", command=self.browse_file)
        self.browse_button.pack(pady=10)

        self.upload_button = tk.Button(self.root, text="Upload", command=self.upload_file)
        self.upload_button.pack(pady=10)

        self.root.bind("<KeyPress>", self.key_interactions)
        self.network_window = None

    def browse_file(self):
        file_path = filedialog.askopenfilename(filetypes=[("JSON files", "*.json")])
        self.entry_var.set(file_path)

    def key_interactions(self,event):
        if event.keycode == 27:
            self.root.destroy()
    def upload_file(self):
        file_path = self.entry_var.get()
        try:
            with open(file_path, 'r') as file:
                # Read JSON data
                json_data = json.load(file)

                # Open a new window for JSON data visualization
                self.parse_json_metabolic_network_data(json_data)
                self.met_network = MetabolicNetwork(self.list_nodes_metabolites, self.list_edges_reactions)
                # self.upload_button.destroy()
                # self.browse_button.destroy()
                self.root.destroy()
                self.root = tk.Tk()
                self.create_listboxes()
                self.create_slider()
                self.create_load_network_button()
                self.root.bind("<KeyPress>", self.key_interactions)
        except FileNotFoundError:
            print("File not found.")

    def parse_json_metabolic_network_data(self, json_data):
        self.list_edges_reactions = []
        self.list_nodes_metabolites = []

        for key, value in json_data.items():
            list_in_dic = value  # get the list of metabolites or an empty list if not present
            if key == "metabolites":
                for metabolite in list_in_dic:
                    # access information about each metabolite
                    metabolite_id = metabolite['id']
                    metabolite_name = metabolite['name']
                    metabolite_compartment = metabolite['compartment']
                    metabolite_name = f"%s[%s]" % (metabolite_name, metabolite_compartment)
                    self.list_nodes_metabolites.append(MetaboliteNodes(metabolite_id, metabolite_name, metabolite_compartment))
            if key == "reactions":
                for reaction in list_in_dic:
                    reaction_id = reaction["id"]
                    reaction_name = reaction["name"]
                    reaction_metabolites = reaction["metabolites"]
                    reaction_lb = reaction["lower_bound"]
                    reaction_ub = reaction["upper_bound"]
                    self.list_edges_reactions.append(ReactionEdges(reaction_id, reaction_name, reaction_metabolites, reaction_lb, reaction_ub))


    def create_load_network_button(self):
        self.load_network_button = tk.Button(self.root, text = "(re)load network", command = self.load_network )
        self.load_network_button.grid(row=1, column= 2, padx =10, pady=5)

    def load_network(self):
        self.met_network.calculate_network()
        self.create_and_draw_network_canvas()

    def create_and_draw_network_canvas(self):
        # Replace this function with your network drawing logic
        if self.network_window is None:
            self.network_window = tk.Toplevel(self.root)
            self.network_window.geometry("800x900+900+0")
        try:
            self.canvas.place_forget()
            self.canvas.delete("all")
        except:
            pass

        self.canvas = tk.Canvas(self.network_window, width=1200, height=900)
        self.draw_network(self.canvas)
        self.canvas.place(in_=self.network_window, x=0, y=0)


        #nodes =
        #edges
    def create_listboxes(self):

        label1 = tk.Label(self.root, text="included")
        label1.grid(row= 0, column = 0 , padx = 10, pady=5)
        self.left_list = tk.Listbox(self.root, width = 20, height = 20)
        self.left_list.grid(row= 1, column = 0 , padx = 10, pady=5)
        label2 = tk.Label(self.root, text="not included")
        label2.grid(row= 0, column = 1 , padx = 10, pady=5)
        self.right_list = tk.Listbox(self.root, width = 20, height = 20)
        self.right_list.grid(row= 1, column = 1 , padx = 10, pady=5)

        self.names_left = [name.name for name in self.met_network.included_metabolites]
        self.names_right = [name.name for name in self.met_network.not_included_metabolites]

        for item in self.names_left:
            self.left_list.insert(tk.END, item)
        for item in self.names_right:
            self.right_list.insert(tk.END, item)

        self.left_list.bind('<Double-Button-1>', lambda event: self.on_double_click(event, self.left_list, self.right_list))
        self.right_list.bind('<Double-Button-1>', lambda event: self.on_double_click(event,self.right_list, self.left_list))

    def move_item(self,event, from_list, to_list):
        selected_item = from_list.get(from_list.curselection())
        to_list.insert(tk.END, selected_item)
        from_list.delete(from_list.curselection())

    def on_double_click(self,event, from_list, to_list):
        self.move_item(event, from_list, to_list)

    def update_slider(self, value):
        value = int(value)
        self.met_network.slider_value = value
        self.right_list.destroy()
        self.left_list.destroy()
        self.met_network.update_included_not_included_based_on_slider()
        self.create_listboxes()

    def create_slider(self):
        scale = tk.Scale(self.root, from_ = self.met_network.slider_start, to = self.met_network.slider_end, orient=tk.HORIZONTAL, command = self.update_slider)
        scale.grid(row=2, column = 0, columnspan = 2, padx = 10, pady = 5, sticky = "ew")
        scale.set(self.met_network.slider_value)

    def draw_network(self,canvas):
        print(len(self.met_network.network_node_positions))
        # Draw nodes
        node_radius = 5
        for node, (x, y) in self.met_network.network_node_positions.items():
            canvas.create_oval(x - node_radius, y - node_radius, x + node_radius, y + node_radius, fill='blue')
            canvas.create_text(x, y - 15, text=node)
        # Draw edges
        for start, end, bla in self.met_network.network_edges:
            x1, y1 = self.met_network.network_node_positions[start]
            x2, y2= self.met_network.network_node_positions[end]
            canvas.create_line(x1, y1, x2, y2, fill='black')


# Create the main Tkinter window
root = tk.Tk()

# Create an instance of the JSONViewer class
TK_application = Application(root)

# Start the Tkinter event loop
root.mainloop()


### CURRENT ISSUE:
# THE code does not correctly keep the position of the different fixed nodes, even when nothing changes, potentially the order in which they are read is different?
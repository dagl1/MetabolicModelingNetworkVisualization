# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
import os, json, pygame, math, random, markdown
import numpy as np
import pandas as pd
import tkinter as tk
from tkinter import filedialog
import threading
import queue
import openpyxl
import pickle
import copy
import re
import sys
from datetime import datetime

def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print_hi('PyCharm')


# See PyCharm help at https://www.jetbrains.com/help/pycharm/
#
#
# class MainLoop:
#
#     def __init__(self):
#         self.WIDTH, self.HEIGHT = 1200, 800
#         self.FINAL_WIDTH, self.FINAL_HEIGHT = 2400, 1600
#         pygame.init()
#         self.canvas = pygame.Surface((self.WIDTH, self.HEIGHT))
#         self.screen = pygame.display.set_mode((self.FINAL_WIDTH, self.FINAL_HEIGHT))
#
#
# class TkinterSetup:
#     pass
#
#
# filepath = "E:\Git\MaCSBio-GEM\General\Functions\Metabolic_tasks\Task322Sample1Alpha0.50Beta0.00_Date2023-10-25_21-51-37\ModelsJSON/newModelMILP_OR_AdjustTask322Sample1Scored2.152average4.105alpha0.500beta0.000.json"
#
# with open(filepath, 'r') as file:
#     # Read JSON data
#     json_data = json.load(file)
#

class MetaboliteNodes:
    def __init__(self, id, name, compartment):
        self.id = id
        self.name = name
        self.compartment = compartment
        self.x, self.y = 0, 0
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
        self.visible = True  # for compartment
        self.flux = 0
        self.expression = 0


class CompartmentNode:
    def __init__(self, id, name):
        self.id = id
        self.name = name
        self.metabolites = []
        self.visible = False
        self.flux = 0
        self.expression = 0


class BoundaryNode:
    def __init__(self):
        pass


class NetworkNodesForPhysicsSimulation:
    def __init__(self, x, y, instance, id_number, met_network):
        self.x, self.y = x, y
        self.instance = instance
        self.x_vector, self.y_vector = 0, 0
        self.id_number = id_number
        self.met_network = met_network
        if isinstance(self.instance, ReactionEdges):
            self.node_type = "reaction"
            self.fixed = False
        elif isinstance(self.instance, MetaboliteNodes):
            self.node_type = "metabolite"
            if self.instance.fixed_node_input or self.instance.fixed_node_output:
                self.fixed = True
            else:
                self.fixed = False
        elif isinstance(self.instance, CompartmentNode):
            self.node_type = "compartment"
            self.fixed = False

        elif isinstance(self.instance, BoundaryNode):
            self.node_type = "boundary"
            self.fixed = True

        self.input_instances = []
        self.output_instances = []
        self.find_input_and_output_instances()

        self.rect = pygame.Rect(self.x, self.y, 15, 15)
        self.fixed_by_mouse = False
        self.visible_by_mouse = True
        # Unfortunately necessary to make sure all lines are drawn (something goes wrong with indicing
        if isinstance(self.instance, MetaboliteNodes) and self.instance.id == "000fake":
            self.visible_by_mouse = False
        self.highlighted_by_mouse = False
        self.split_of_nodes_partners =  []
        self.selected_by_drag_selection = False
    def find_input_and_output_instances(self):
        if isinstance(self.instance, ReactionEdges):
            for metabolite, stoichiometry in self.instance.metabolites.items():
                if metabolite in [metabolit.id for metabolit in self.met_network.included_metabolites] and stoichiometry > 0:
                    self.input_instances.append(self.met_network.find_metabolite(metabolite))
                elif metabolite in [metabolit.id for metabolit in self.met_network.included_metabolites] and stoichiometry < 0:
                    self.output_instances.append(self.met_network.find_metabolite(metabolite))
        elif isinstance(self.instance, MetaboliteNodes):
            for reaction in self.instance.reactions:
                reaction_ = self.met_network.find_reaction(reaction)
                for metabolite, stoichiometry in reaction_.metabolites.items():
                    if metabolite == self.instance.id and stoichiometry < 0:
                        self.input_instances.append(reaction_)
                    elif metabolite == self.instance.id and stoichiometry > 0:
                        self.output_instances.append(reaction_)
        elif isinstance(self.instance, CompartmentNode):
            for metabolite in self.instance.metabolites:
                if metabolite.id in [metabolit.id for metabolit in self.met_network.included_metabolites]:
                    self.input_instances.append(metabolite)


class MetabolicNetwork:
    def __init__(self, nodes, edges, tkinter_app):
        self.metabolites_nodes = nodes
        self.reaction_edges = edges
        self.tk_application = tkinter_app
        self.set_amount_of_reactions_per_node()
        self.set_slider()
        self.filter_input_output_reactions()
        self.set_reactions_for_nodes()
        self.commonly_excluded_names = []  # ["H+", "H2O", "NADH", "NAD+", "NADP+", "NADPH", "FAD", "FADH", "CO2"]
        self.included_metabolites = []
        self.not_included_metabolites = []
        self.part_of_physics_metabolites = []
        self.not_part_of_physics_metabolites = []
        self.create_lookup_dicts()
        self.update_included_not_included_based_on_slider()
        self.edge_logic_dict = {}
        self.all_network_nodes = []
        self.fixed_metabolite_repulsion_value = 1000
        self.metabolite_compartment_attraction_value = 3000
        self.rxn_met_attraction_value = 3000
        self.boundary_repulsion_value = 3000
        self.comp_comp_repulsion_value = 3000
        self.normal_repulsion_value = 1000
        self.rect_size = 15
        self.font_size = 18
        self.show_lines = True
        self.show_fixed_metabolite_names = "No Names"
        self.show_metabolite_names = "No Names"
        self.show_reaction_names = "Reaction Names"
        self.show_compartments = False
        self.show_flux_expression = "None"
        self.visible_nodes =[]
        self.not_visible_nodes = []
        self.queue = queue.Queue()
        self.df_fluxes_reactions = pd.DataFrame()
        self.fluxes_dict = {}
        self.expression_dict = {}
        self.full_json_dict = {}
        self.zoom_factor = 1
        self.show_color_lines = False
        self.drawing_selection_rect_on_screen = False
        self.begin_selection_rect_pos = (0,0)
        self.end_selection_rect_pos = (0,0)
        self.search_mode_enabled = False
        self.search_input_text = ""
        self.search_suggestion_list = []
        self.full_search_list =[]
        self.current_search_index = 0
        self.selected_search_item = None
        self.show_expressionless_reaction_highlight = False
        self.canvas_width = 900
        self.canvas_height = 800
        self.current_command = None
        self.previous_commands = []
        self.next_commands_undone =[]
        self.control_groups_for_selection = {1:[],2:[],3:[],4:[],5:[],6:[],7:[],8:[],9:[]}
        self.start_command_recording = True

    def save_as_json(self):
        node_xy_save_dict = {}
        for node in self.all_network_nodes:
            if not isinstance(node.instance, BoundaryNode):
                node_xy_save_dict[node.instance.id] = (node.x,node.y)

        # self.fluxes_dict
        # self.expression_dict
        # freeze all mets at loading
        # xy coordinates of each node
        self.full_json_dict = {}
        self.full_json_dict["fluxes_dict"] = self.fluxes_dict
        self.full_json_dict["expression_dict"] = self.expression_dict
        self.full_json_dict["node_xy_save_dict"] = node_xy_save_dict

        self.formulate_metadata_dict()
        self.formulate_main_ESHER_dict()
        self.remove_duplicate_main_nodes()
        self.list_with_dicts_to_save = [ self.meta_data_dict, self.main_ESHER_dict]
        #self.formulate_labels_dict()
        #self.formulate_canvas_dict()

        pass


    def formulate_metadata_dict(self):
        self.meta_data_dict = {}
        date = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
        self.meta_data_dict = {
            "map_name": "new_map",
            "map_id": "000000000",
            "map_description": f"\nLast Modified {date}",
            "homepage": "https://escher.github.io", "schema": "https://escher.github.io/escher/jsonschema/1-0-0#",
            "not_escher_dict": self.full_json_dict
        }

    def formulate_main_ESHER_dict(self):
        self.main_ESHER_dict ={}
        node_counter =0
        segment_counter = 0
        nodes_dict = {}
        reaction_counter  = 0
        expansion_constant = 5
        self.main_ESHER_dict["reactions"]={}
        for node in self.all_network_nodes:
            segments_dict = {}
            if node.node_type == "reaction" and node.visible_by_mouse:
                # calculate closests input and output nodes and designate them as main, all else will be build separately
                distance = 100000000000
                closest_input_node = None
                for input_instance in node.input_instances:
                    for temp_node in self.all_network_nodes:
                        if temp_node.node_type == "metabolite"and input_instance==temp_node.instance and temp_node.visible_by_mouse:
                            input_node = temp_node
                            break
                    node_distance =  math.sqrt((input_node.x - node.x)**2 + (input_node.y - node.y)**2)
                    if closest_input_node is None or node_distance < distance:
                        distance = node_distance
                        closest_input_node = input_node

                distance = 10000000000000
                closest_output_node = None
                for output_instance in node.output_instances:
                    for temp_node in self.all_network_nodes:
                        if temp_node.node_type == "metabolite"and output_instance == temp_node.instance and temp_node.visible_by_mouse:
                            output_node = temp_node
                    node_distance =  math.sqrt((output_node.x - node.x)**2 + (output_node.y - node.y)**2)
                    if closest_output_node is None or node_distance<distance:
                        distance = node_distance
                        closest_output_node = output_node


                ## create node and reaction entry

                start_node_counter = node_counter
                main_x, main_y = node.x, node.y

                # static segments:
                segments_dict[segment_counter] = {
                    "from_node_id": str(start_node_counter),
                    "to_node_id": str(start_node_counter+1),
                    "b1": None,
                    "b2": None
                }
                segment_counter += 1

                segments_dict[segment_counter] = {
                    "from_node_id": str(start_node_counter+1),
                    "to_node_id": str(start_node_counter+2),
                    "b1": None,
                    "b2": None
                }
                segment_counter += 1

                nodes_dict[node_counter] = {
                    "node_type":"multimarker",
                    "x": expansion_constant*main_x,
                    "y": expansion_constant*(main_y) -40,
                }
                node_counter +=1

                nodes_dict[node_counter] = {
                    "node_type":"midmarker",
                    "x": expansion_constant*main_x,
                    "y": expansion_constant*(main_y) -10,
                }
                node_counter +=1

                nodes_dict[node_counter] = {
                    "node_type":"multimarker",
                    "x": expansion_constant*main_x,
                    "y": expansion_constant*(main_y) +20,
                }
                node_counter +=1

                before_input_counter = node_counter
                if closest_output_node is not None:
                    nodes_dict[node_counter] = {
                        "node_type": "metabolite",
                        "x": expansion_constant * closest_output_node.x,
                        "y": expansion_constant * closest_output_node.y,
                        "bigg_id": closest_output_node.instance.id,
                        "name": closest_output_node.instance.name,
                        "label_x": expansion_constant * closest_output_node.x - 20,
                    "label_y": expansion_constant * closest_output_node.y - 20,
                    "node_is_primary": True
                    }
                    node_counter +=1

                    segments_dict[segment_counter] = {
                        "from_node_id": str(node_counter-1),
                        "to_node_id": str(start_node_counter),
                        "b1": { "x": expansion_constant * closest_output_node.x+50,
                                "y": expansion_constant * closest_output_node.y+50},
                        "b2": {"x": expansion_constant * main_x,
                               "y": expansion_constant * (main_y)-40 }
                    }
                    segment_counter += 1


                for idx, input_instance in enumerate(node.output_instances):
                    if input_instance != closest_output_node.instance:
                        nodes_dict[node_counter] = {
                            "node_type": "metabolite",
                            "x": expansion_constant* (main_x) + (idx*40)-30,
                            "y": expansion_constant* (main_y)-120,
                            "bigg_id": input_instance.id,
                            "name": input_instance.name,
                            "label_x": expansion_constant *(main_x) + (idx*40) -50,
                            "label_y": expansion_constant * (main_y)-120 - 20,
                            "node_is_primary": False
                        }
                        node_counter += 1

                        segments_dict[segment_counter] = {
                            "from_node_id": str(node_counter - 1),
                            "to_node_id": str(start_node_counter),
                            "b1": {"x": expansion_constant * (main_x) + (idx*40)-30,
                                   "y":  expansion_constant * (main_y)-120},
                            "b2": {"x": expansion_constant*main_x ,
                                   "y": expansion_constant* (main_y) - 40}
                        }
                        segment_counter += 1

                before_output_node_counter = node_counter

                if closest_input_node is not None:
                    nodes_dict[node_counter] = {
                        "node_type": "metabolite",
                        "x": expansion_constant* closest_input_node.x,
                        "y": expansion_constant* closest_input_node.y,
                        "bigg_id": closest_input_node.instance.id,
                        "name": closest_input_node.instance.name,
                        "label_x": expansion_constant * closest_input_node.x - 20,
                        "label_y": expansion_constant * closest_input_node.y + 20,
                        "node_is_primary": True
                    }
                    node_counter += 1
                    segments_dict[segment_counter] = {
                        "from_node_id": str(start_node_counter + 2),
                        "to_node_id": str(node_counter - 1),
                        "b1": {"x":expansion_constant*main_x,
                               "y":expansion_constant*main_y + 20},
                        "b2": { "x": expansion_constant * closest_input_node.x-50,
                                "y": expansion_constant * closest_input_node.y-50}
                    }
                    segment_counter += 1

                for idx, input_instance in enumerate(node.input_instances):
                    if input_instance != closest_input_node.instance:
                        nodes_dict[node_counter] = {
                            "node_type": "metabolite",
                            "x": expansion_constant*main_x + (idx * 40)-30,
                            "y": expansion_constant*main_y + 120,
                            "bigg_id": input_instance.id,
                            "name": input_instance.name,
                            "label_x": expansion_constant*main_x + (idx * 40) - 50,
                            "label_y": expansion_constant*main_y +120 - 20,
                            "node_is_primary": False
                        }
                        node_counter += 1

                        segments_dict[segment_counter] = {
                            "from_node_id": str(start_node_counter + 2),
                            "to_node_id": str(node_counter-1),
                            "b1": {"x":expansion_constant*main_x,
                                   "y":expansion_constant*main_y +20},
                            "b2": {"x":expansion_constant*main_x + (idx * 40)-30-50,
                                   "y": expansion_constant*main_y+120 -50}
                        }
                        segment_counter += 1
                metabolites_list_with_dicts = []
                for met, stoichiometry in node.instance.metabolites.items():
                    dict_to_add = {"bigg_id": met, "coefficient": stoichiometry}
                    metabolites_list_with_dicts.append(dict_to_add)


                self.main_ESHER_dict["reactions"][reaction_counter] = {
                    "name": node.instance.name,
                    "bigg_id": node.instance.id,
                    "reversibility": False,
                    "label_x": expansion_constant*main_x + 50,
                    "label_y": expansion_constant*main_y,
                    "gene_reaction_rule": "",
                    "genes": [{}],
                    "metabolites": metabolites_list_with_dicts,
                    "segments" : segments_dict
                }
                reaction_counter+=1
        self.main_ESHER_dict["nodes"]=nodes_dict
        self.main_ESHER_dict["text_labels"] = {}
        self.main_ESHER_dict["canvas"] = {
            "x": -1000,
            "y": -500,
            "width": 10000,
            "height": 7000
        }

    def remove_duplicate_main_nodes(self):
        list_to_remove = []
        for id, node_dict in self.main_ESHER_dict["nodes"].items():
            if node_dict["node_type"] == "metabolite" and node_dict["node_is_primary"]:
                list_with_ids = []
                for id2, node_dict2 in self.main_ESHER_dict["nodes"].items():
                    if id != id2 and node_dict2["node_type"] == "metabolite" and node_dict2["node_is_primary"] and node_dict2["bigg_id"] == node_dict["bigg_id"]:
                        list_with_ids.append(id2)
                        node_dict2["bigg_id"] = f"removed {id2}"
                for reaction_id, reaction_dict in self.main_ESHER_dict["reactions"].items():
                    segments = reaction_dict["segments"]
                    for value in segments.values():
                        #print(value,list_with_ids)
                        if int(value["from_node_id"]) in list_with_ids:
                            value["from_node_id"] = id
                        if int(value["to_node_id"]) in list_with_ids:
                            value["to_node_id"] = id
                list_to_remove.extend(list_with_ids)
        for key in list_to_remove:
            if key in self.main_ESHER_dict["nodes"]:
                del self.main_ESHER_dict["nodes"][key]



    def find_main_metabolite_nodes(self):
         self.big_metabolite_nodes = []
         for node in self.all_network_nodes:
             if node.node_type == "reaction" and node.visible_by_mouse:
                 reaction_instance = node.instance
                 mets = reaction_instance.metabolites

                 lowest_connection_input_reaction_number = None
                 lowest_connection_output_reaction_number = None
                 lowest_connection_input_metabolite = None
                 lowest_connection_output_metabolite = None
                 for metabolite, stoichiometry in mets.items():
                     metabolite_instance = self.find_metabolite(metabolite)
                     if stoichiometry < 0: #input reaction
                         amount_of_reactions_where_metabolite_is_input = self.determine_true_amount_of_reactions(metabolite_instance, "input")
                         if lowest_connection_input_reaction_number == None or lowest_connection_input_reaction_number > amount_of_reactions_where_metabolite_is_input:
                             lowest_connection_input_reaction_number = amount_of_reactions_where_metabolite_is_input
                             lowest_connection_input_metabolite = metabolite_instance
                     else:
                         amount_of_reactions_where_metabolite_is_output = self.determine_true_amount_of_reactions(metabolite_instance, "output")
                         if lowest_connection_output_reaction_number == None or lowest_connection_output_reaction_number > amount_of_reactions_where_metabolite_is_output:
                             lowest_connection_output_reaction_number = amount_of_reactions_where_metabolite_is_output
                             lowest_connection_output_metabolite = metabolite_instance
                 node.main_metabolite_connections = [lowest_connection_input_metabolite, lowest_connection_output_metabolite]


    def determine_true_amount_of_reactions(self, metabolite_instance, input_or_output):
        name_to_check = re.sub(r'\[.\]$', '', metabolite_instance.name)
        other_instances_with_similiar_name = [metabolite_instance]
        for metabolite in self.metabolites_nodes:
            if re.sub(r'\[.\]$', '', metabolite.name) == name_to_check:
                other_instances_with_similiar_name.append(metabolite)
        amount_of_reactions = 0
        if input_or_output == "input":
            for instance in other_instances_with_similiar_name:
                for reaction in self.reaction_edges:
                    for met, stoichiometry in reaction.metabolites.items():
                        if stoichiometry < 0 and instance.id == met:
                            amount_of_reactions +=1
            return amount_of_reactions

        else:
            for instance in other_instances_with_similiar_name:
                for reaction in self.reaction_edges:
                    for met, stoichiometry in reaction.metabolites.items():
                        if stoichiometry > 0 and instance.id == met:
                            amount_of_reactions +=1
            return amount_of_reactions

    def formulate_labels_dict(self):
        self.main_ESHER_dict["text_labels"] = {}
    def formulate_canvas_dict(self):
        self.main_ESHER_dict["canvas"] = {
            "x": -2000, "y": -1000, "width": 10000, "height": 10000
        }

    def read_from_json(self,data):
        escher_exists = False
        visualization_network_exists = False
        for dicts in data:
            if "not_escher_dict" in dicts:
                visualization_network_exists = True
                for key_, element in dicts["not_escher_dict"].items():
                    self.full_json_dict[key_] = element
                break

        if visualization_network_exists:
            self.calculate_network(visualization_network_exists)

    def filter_input_output_reactions(self):
        reaction_edges_copy = []
        for reaction in self.reaction_edges:
            if not (reaction.lb > 0 and reaction.ub < 1000) or not reaction.name.startswith('temporary_exchange_'):
                reaction_edges_copy.append(reaction)
            else:
                self.set_fixed_metabolite(reaction)
        self.reaction_edges = reaction_edges_copy

    def set_fixed_metabolite(self, reaction):
        mets = reaction.metabolites
        for metabolite, stoichiometry in mets.items():
            for met in self.metabolites_nodes:
                if metabolite == met.id:
                    if stoichiometry > 0:
                        met.fixed_node_input = True
                    if stoichiometry < 0:
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
        if metabolite_id.startswith("MAM0"):  # should be changed if metabolite IDs ever go past 10000
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
        self.part_of_physics_metabolites=[]
        self.not_part_of_physics_metabolites = []
        self.find_reactions_that_should_generally_be_excluded()
        for metabolite in self.metabolites_nodes:
            if metabolite.amount_of_reactions > self.slider_value and metabolite not in self.not_part_of_physics_metabolites:
                self.not_part_of_physics_metabolites.append(metabolite)
        for metabolite in self.metabolites_nodes:
            if metabolite not in self.included_metabolites:
                self.included_metabolites.append(metabolite)
        for metabolite in self.metabolites_nodes:
             if (metabolite.fixed_node_input or metabolite.fixed_node_output) and metabolite not in self.part_of_physics_metabolites:
                self.part_of_physics_metabolites.append(metabolite)
                try:
                    self.not_part_of_physics_metabolites.remove(metabolite)
                except:
                    pass
        self.part_of_physics_metabolites = [item for item in self.metabolites_nodes if item not in self.not_part_of_physics_metabolites]
        self.not_part_of_physics_metabolites= [item for item in self.metabolites_nodes if item not in self.part_of_physics_metabolites]


    def set_reactions_for_nodes(self):
        for reaction in self.reaction_edges:
            # if not (reaction.lb > 0 and reaction.ub < 1000) or not reaction.name.startswith('temporary_exchange_'):
            for key, value in reaction.metabolites.items():
                for metabolite in self.metabolites_nodes:
                    if key == metabolite.id:
                        metabolite.reactions.append(reaction.id)
                        break

    def set_amount_of_reactions_per_node(self):
        # dont count temporary reactions
        for reaction in self.reaction_edges:
            # if not (reaction.lb > 0 and reaction.ub < 1000) or not reaction.name.startswith('temporary_exchange_'):
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

    def restart_simulation(self):
        simulation_thread = threading.Thread(target=self.run_network_simulation)
        simulation_thread.start()

    def calculate_network(self, load_from_json=False):
        # Create compartment_pseudo_reaction_nodes:
        self.create_compartment_pseudo_reaction_nodes()
        self.create_network_nodes_for_simulation()

        # create edge matrices
        self.create_edges_matrices_dicts()
        self.update_edge_network()
        # set initial x, y positions for non-boundary nodes
        self.set_x_y_positions_non_boundary_nodes()
        if load_from_json:
            for node in self.all_network_nodes:
                if node.node_type != "boundary":
                    key = node.instance.id
                    x,y = self.full_json_dict["node_xy_save_dict"][key]
                    node.x = x
                    node.y = y
            for node in self.all_network_nodes:
                node.fixed_by_mouse = True


        # run network
        simulation_thread = threading.Thread(target=self.run_network_simulation)
        simulation_thread.start()

        nodes, edges = self.extract_non_pseudo_nodes_and_edges()
        return nodes, edges

    def create_compartment_pseudo_reaction_nodes(self):
        compartments = []
        self.compartments = []
        for metabolite in self.included_metabolites:
            if metabolite.compartment not in compartments:
                compartments.append(metabolite.compartment)
        for idx, compartment in enumerate(compartments):
            self.compartments.append(CompartmentNode(compartment, compartment))
            for metabolite in self.included_metabolites:
                if metabolite.compartment == compartment:
                    self.compartments[idx].metabolites.append(metabolite)

    def create_network_nodes_for_simulation(self):
        # Add a fake metabolite, unfortunately necessary for making sure all lines are drawn correttly (weirdly enough)
        if self.included_metabolites[0].id != "000fake":
            self.included_metabolites.insert(0, MetaboliteNodes("000fake", "000fake", self.compartments[0].id))

        self.all_network_nodes = []
        # Create all types of nodes
        id_number = 0
        for metabolite in self.included_metabolites:
            self.all_network_nodes.append(NetworkNodesForPhysicsSimulation(0, 0, metabolite, id_number, self))
            if metabolite in [met for met in self.not_part_of_physics_metabolites]:
                self.all_network_nodes[-1].visible_by_mouse = False
            id_number += 1

        for reaction in self.reaction_edges:
            self.all_network_nodes.append(NetworkNodesForPhysicsSimulation(0, 0, reaction, id_number, self))
            id_number += 1

        for compartment in self.compartments:
            self.all_network_nodes.append(NetworkNodesForPhysicsSimulation(0, 0, compartment, id_number, self))
            id_number += 1

        # Create an edge with boundary repulsors (will have a very high power so their value will fall off quickly, just to make sure network stays within the
        # bounds of the canvas, takes global height and width inputted as canvas as I didn't want to do some weird directing
        global width
        global height
        amount_of_boundary_rows = (height // 50) + 1
        amount_of_boundary_columns = (width // 50) + 1
        for boundary_idx in range(amount_of_boundary_columns):
            self.all_network_nodes.append(NetworkNodesForPhysicsSimulation(boundary_idx * 50, 0, BoundaryNode(), id_number, self))
            id_number += 1
            self.all_network_nodes.append(NetworkNodesForPhysicsSimulation(boundary_idx * 50, height, BoundaryNode(), id_number, self))
            id_number += 1

        for boundary_idx in range(amount_of_boundary_rows):
            self.all_network_nodes.append(NetworkNodesForPhysicsSimulation(0, boundary_idx * 50, BoundaryNode(), id_number, self))
            id_number += 1
            self.all_network_nodes.append(NetworkNodesForPhysicsSimulation(width, boundary_idx * 50, BoundaryNode(), id_number, self))
            id_number += 1

    def create_edges_matrices_dicts(self):
        node_length = len(self.all_network_nodes)

        temp_connection_matrix = np.zeros((node_length, node_length))
        edges = []
        for node in self.all_network_nodes:
            edge = [node.id_number]
            if node.node_type == "reaction" or node.node_type == "metabolite":
                for input_node in node.input_instances:
                    for temp_node in self.all_network_nodes:
                        if input_node == temp_node.instance:
                            connected_node_id_number = temp_node.id_number
                            edge.append(connected_node_id_number)
                            break
                for output_node in node.output_instances:
                    for temp_node in self.all_network_nodes:
                        if output_node == temp_node.instance:
                            connected_node_id_number = temp_node.id_number
                            edge.append(connected_node_id_number)
                            break
            edges.append(edge)
        for edge in edges:
            node_id = edge[0]
            connected_nodes = edge[1:]
            temp_connection_matrix[node_id, connected_nodes] = 1
        self.edge_logic_dict["reaction_metabolite_attraction"] = temp_connection_matrix

        temp_connection_matrix = np.zeros((node_length, node_length))
        edges = []
        for node in self.all_network_nodes:
            edge = [node.id_number]
            if node.node_type == "compartment":
                for other_node in self.all_network_nodes:
                    if other_node != node and other_node.node_type == "compartment":
                        edge.append(other_node.id_number)
            edges.append(edge)
        for edge in edges:
            node_id = edge[0]
            connected_nodes = edge[1:]
            temp_connection_matrix[node_id, connected_nodes] = 1
        self.edge_logic_dict["compartment_compartment_repulsion"] = temp_connection_matrix

        temp_connection_matrix = np.zeros((node_length, node_length))
        edges = []
        for node in self.all_network_nodes:
            edge = [node.id_number]
            if node.node_type == "metabolite" and node.instance.fixed_node_input:
                for other_node in self.all_network_nodes:
                    if node != other_node and other_node.node_type == "metabolite" and other_node.instance.fixed_node_input:
                        edge.append(other_node.id_number)

            if node.node_type == "metabolite" and node.instance.fixed_node_output:
                for other_node in self.all_network_nodes:
                    if other_node.node_type == "metabolite" and other_node.instance.fixed_node_output:
                        edge.append(other_node.id_number)
            edges.append(edge)
        for edge in edges:
            node_id = edge[0]
            connected_nodes = edge[1:]

            temp_connection_matrix[node_id, connected_nodes] = 1
        self.edge_logic_dict["fixed_metabolite_repulsion"] = temp_connection_matrix

        temp_connection_matrix = np.zeros((node_length, node_length))
        edges = []
        for node in self.all_network_nodes:
            edge = [node.id_number]
            for boundary_node in self.all_network_nodes:
                if boundary_node.node_type == "boundary":
                    edge.append(boundary_node.id_number)
            edges.append(edge)
        for edge in edges:
            node_id = edge[0]
            connected_nodes = edge[1:]
            temp_connection_matrix[node_id, connected_nodes] = 1
        self.edge_logic_dict["boundary_repulsion"] = temp_connection_matrix

        temp_connection_matrix = np.zeros((node_length, node_length))
        edges = []
        for node in self.all_network_nodes:
            edge = [node.id_number]
            if node.node_type != "boundary":
                for other_node in self.all_network_nodes:
                    if other_node != node and other_node.node_type != "boundary":
                        edge.append(other_node.id_number)
            edges.append(edge)
        for edge in edges:
            node_id = edge[0]
            connected_nodes = edge[1:]

            temp_connection_matrix[node_id, connected_nodes] = 1
        self.edge_logic_dict["normal_repulsion"] = temp_connection_matrix

        temp_connection_matrix = np.zeros((node_length, node_length))
        edges = []
        for node in self.all_network_nodes:
            edge = [node.id_number]
            if node.node_type == "compartment":
                for connected_instance in node.input_instances:
                    for other_node in self.all_network_nodes:
                        if connected_instance == other_node.instance and not (other_node.instance.fixed_node_input or other_node.instance.fixed_node_output):
                            edge.append(other_node.id_number)
            edges.append(edge)
        for edge in edges:
            node_id = edge[0]
            connected_nodes = edge[1:]
            temp_connection_matrix[node_id, connected_nodes] = 1
            for extra_edge in edge[1:]:
                temp_connection_matrix[extra_edge, node_id] = 1
        self.edge_logic_dict["metabolite_compartment_attraction"] = temp_connection_matrix

    def update_edge_network(self):
        self.create_edges_matrices_dicts()
        idx_list = []

        for node in self.all_network_nodes:
            if node.node_type == "metabolite" and node.instance in [met for met in self.not_part_of_physics_metabolites]:
                idx_list.append(node.id_number)
        for key,matrix in self.edge_logic_dict.items():
            if key != "draw_lines_connections":
                for i in idx_list:
                    matrix[i, :] = 0  # Zero out rows
                    matrix[:, i] = 0  # Zero out columns
                self.edge_logic_dict[key] = matrix
        self.queue.put(self.tk_application.create_listboxes)
        self.queue.put(self.tk_application.create_listboxes_visibility)
        node_length = len(self.all_network_nodes)

        temp_connection_matrix = np.zeros((node_length, node_length))
        edges = []
        for node in self.all_network_nodes:
            edge = [node.id_number]
            if node.node_type == "reaction" or node.node_type == "metabolite":
                # for input_node in node.input_instances:
                #     for temp_node in self.all_network_nodes:
                #         if input_node == temp_node.instance:
                #             connected_node_id_number = temp_node.id_number
                #             edge.append(connected_node_id_number)
                #             break
                for output_node in node.output_instances:
                    for temp_node in self.all_network_nodes:
                        if output_node == temp_node.instance:
                            connected_node_id_number = temp_node.id_number
                            edge.append(connected_node_id_number)
                            break
            edges.append(edge)

        for edge in edges:
            node_id = edge[0]
            connected_nodes = edge[1:]
            temp_connection_matrix[node_id, connected_nodes] = 1
        self.edge_logic_dict["draw_lines_connections"] = temp_connection_matrix

    def set_x_y_positions_non_boundary_nodes(self):
        fixed_input_counter = 0
        fixed_output_counter = 0
        global width
        global height
        random.seed(666)
        for node in self.all_network_nodes:
            if node.node_type == "boundary":
                pass
            elif node.node_type == "metabolite" and node.instance.fixed_node_output:
                node.x = 200 + (fixed_output_counter * 25)
                fixed_output_counter += 1
                node.y = height - 50
            elif node.node_type == "metabolite" and node.instance.fixed_node_input:
                node.x = 200 + (fixed_input_counter * 25)
                fixed_input_counter += 1
                node.y = 0 + 50
            else:
                node.x = random.randint(200, width - 200)
                node.y = random.randint(200, height - 200)

    def run_network_simulation(self):
        self.network_physics_simulation_loop()

    def network_physics_simulation_loop(self):
        # loop
        try:
            pygame.quit()
        except:
            pass
        global running
        running = True
        pygame_on = True
        if pygame_on:
            os.environ['SDL_VIDEO_WINDOW_POS'] = '100,30' # TODO ADD IN DEFAULT SETTING
            pygame.init()
            display = pygame.display.set_mode([self.canvas_width, self.canvas_height])
            text_input_rect = pygame.Rect(50, 50, 300, 40)
            self.viewport = pygame.Rect(0, 0, 1500, 1200)
            self.font = pygame.font.SysFont(None, self.font_size)
            self.font2 = pygame.font.SysFont(None, 28)
            self.font3 = pygame.font.Font(None,32)
            key_actions = {"left_mouse_clicked": False, "right_mouse_clicked": False, "middle_mouse_clicked": False,
                           "shift_clicked": False, "v_clicked": False, "split_clicked": False, "up_clicked": False, "down_clicked": False,
                           "right_clicked": False, "left_clicked": False, "f_clicked": False, "g_clicked": False, "tab_clicked": False,
                           "enter_clicked": False, "escape clicked": False, "control_clicked": False}
            selected_node = None
        while running:
            if pygame_on:
                display.fill((0, 0, 0))
                for event in pygame.event.get():
                    if event.type == pygame.QUIT:
                        running = False
                    elif event.type == pygame.KEYDOWN:
                        if self.search_mode_enabled:
                            if event.key == pygame.K_ESCAPE:
                                self.search_mode_enabled = False
                            elif event.key == pygame.K_BACKSPACE:
                                self.search_input_text = self.search_input_text[:-1]  # Remove the last character
                            elif event.key == pygame.K_RETURN:
                                if len(self.search_suggestion_list) == 0:
                                    pass
                                elif len(self.search_suggestion_list) <2:
                                    self.current_search_index = 0
                                    self.selected_search_item = self.search_suggestion_list[self.current_search_index]["node"]
                                    self.selected_search_item.selected_by_drag_selection = not self.selected_search_item.selected_by_drag_selection
                                else:
                                    self.selected_search_item = self.search_suggestion_list[self.current_search_index]["node"]
                                    self.selected_search_item.selected_by_drag_selection = not self.selected_search_item.selected_by_drag_selection
                                self.search_input_text = ''  # Clear the input text after processing
                            elif event.key == pygame.K_TAB:
                                #finish search
                                if len(self.search_suggestion_list)>0:
                                    self.search_input_text = self.search_suggestion_list[self.current_search_index]["text"]
                                    self.current_search_index = 0
                            elif event.key == pygame.K_UP:
                                if len(self.search_suggestion_list)==0:
                                    pass
                                elif len(self.search_suggestion_list) ==1:
                                    self.current_search_index =0
                                else:
                                    self.current_search_index -=1
                                    self.current_search_index = self.current_search_index % (len(self.search_suggestion_list))
                            elif event.key == pygame.K_DOWN:

                                if len(self.search_suggestion_list)==0:
                                    pass
                                elif len(self.search_suggestion_list) ==1:
                                    self.current_search_index =0
                                else:
                                    self.current_search_index +=1
                                    self.current_search_index = self.current_search_index % (len(self.search_suggestion_list))

                            else:
                                self.search_input_text += event.unicode  # Add the typed character to the input text
                            self.update_search_recommendations()

                        elif event.key == pygame.K_F1:
                            running = False
                        elif event.key == pygame.K_1 and key_actions["control_clicked"]:
                            self.add_currently_selected_to_control_group(1)
                        elif event.key == pygame.K_1:
                            self.select_control_group(1)

                        elif event.key == pygame.K_2 and key_actions["control_clicked"]:
                            self.add_currently_selected_to_control_group(2)
                        elif event.key == pygame.K_2:
                            self.select_control_group(2)

                        elif event.key == pygame.K_3 and key_actions["control_clicked"]:
                            self.add_currently_selected_to_control_group(3)
                        elif event.key == pygame.K_3:
                            self.select_control_group(3)

                        elif event.key == pygame.K_4 and key_actions["control_clicked"]:
                            self.add_currently_selected_to_control_group(4)
                        elif event.key == pygame.K_4:
                            self.select_control_group(4)

                        elif event.key == pygame.K_5 and key_actions["control_clicked"]:
                            self.add_currently_selected_to_control_group(5)
                        elif event.key == pygame.K_5:
                            self.select_control_group(5)

                        elif event.key == pygame.K_6 and key_actions["control_clicked"]:
                            self.add_currently_selected_to_control_group(6)
                        elif event.key == pygame.K_6:
                            self.select_control_group(6)

                        elif event.key == pygame.K_7 and key_actions["control_clicked"]:
                            self.add_currently_selected_to_control_group(7)
                        elif event.key == pygame.K_7:
                            self.select_control_group(7)

                        elif event.key == pygame.K_8 and key_actions["control_clicked"]:
                            self.add_currently_selected_to_control_group(8)
                        elif event.key == pygame.K_8:
                            self.select_control_group(8)

                        elif event.key == pygame.K_9 and key_actions["control_clicked"]:
                            self.add_currently_selected_to_control_group(9)
                        elif event.key == pygame.K_9:
                            self.select_control_group(9)

                        elif event.key == pygame.K_LCTRL:
                            key_actions["control_clicked"] = True

                        elif event.key == pygame.K_LSHIFT:
                            key_actions["shift_clicked"] = True
                        elif event.key == pygame.K_e:
                            previous_zoom_factor = self.zoom_factor
                            self.zoom_factor += 1
                            if self.zoom_factor > 10:
                                self.zoom_factor = 10
                            current_zoom_factor = self.zoom_factor
                            if previous_zoom_factor != current_zoom_factor:
                                zoom_ratio = current_zoom_factor / previous_zoom_factor
                                self.viewport.x = int((self.viewport.x * zoom_ratio) + (display.get_width()/previous_zoom_factor/2))
                                self.viewport.y = int((self.viewport.y * zoom_ratio) + (display.get_height()/previous_zoom_factor/2))
                        elif event.key == pygame.K_q:
                            previous_zoom_factor = self.zoom_factor
                            self.zoom_factor -= 1
                            if self.zoom_factor < 0.5:
                                self.zoom_factor = 0.5
                            current_zoom_factor = self.zoom_factor
                            if previous_zoom_factor != current_zoom_factor:
                                zoom_ratio = current_zoom_factor / previous_zoom_factor
                                self.viewport.x = int((self.viewport.x * zoom_ratio) - (display.get_width()/previous_zoom_factor/2) )
                                self.viewport.y = int((self.viewport.y * zoom_ratio) - (display.get_height()/previous_zoom_factor/2))
                        elif event.key == pygame.K_v:
                            key_actions["v_clicked"] = True
                        elif event.key == pygame.K_x:
                            key_actions["split_clicked"] = True
                        elif event.key == pygame.K_UP or event.key == pygame.K_w:
                            key_actions["up_clicked"] = True
                        elif event.key == pygame.K_DOWN or event.key == pygame.K_s:
                            key_actions["down_clicked"] = True
                        elif event.key == pygame.K_LEFT or event.key == pygame.K_a:
                            key_actions["left_clicked"] = True
                        elif event.key == pygame.K_RIGHT or event.key == pygame.K_d:
                            key_actions["right_clicked"] = True
                        elif event.key == pygame.K_f:
                            key_actions["f_clicked"] = True
                            if key_actions["shift_clicked"] or key_actions["control_clicked"]:
                                set_to_freeze = False
                                for node in self.all_network_nodes:
                                    if node.selected_by_drag_selection:
                                        set_to_freeze = node.fixed_by_mouse
                                        break
                                for node in self.all_network_nodes:
                                    if node.selected_by_drag_selection:
                                        node.fixed_by_mouse = not set_to_freeze
                        elif event.key == pygame.K_g:
                            self.search_input_text = ""
                            self.search_mode_enabled = True
                        elif event.key == pygame.K_z and key_actions["control_clicked"]:
                            self.undo_last_action()

                        elif event.key == pygame.K_y and key_actions["control_clicked"]:
                            self.redo_last_action()
                    elif event.type == pygame.KEYUP:
                        if event.key == pygame.K_LCTRL:
                            key_actions["control_clicked"] = False
                        elif self.search_mode_enabled:
                            pass
                        elif event.key == pygame.K_LSHIFT:
                            key_actions["shift_clicked"] = False
                        elif event.key == pygame.K_v:
                            key_actions["v_clicked"] = False
                        elif event.key == pygame.K_x:
                            key_actions["split_clicked"] = False
                        elif event.key == pygame.K_UP or event.key == pygame.K_w:
                            key_actions["up_clicked"] = False
                        elif event.key == pygame.K_DOWN or event.key == pygame.K_s:
                            key_actions["down_clicked"] = False
                        elif event.key == pygame.K_LEFT or event.key == pygame.K_a:
                            key_actions["left_clicked"] = False
                        elif event.key == pygame.K_RIGHT or event.key == pygame.K_d:
                            key_actions["right_clicked"] = False
                        elif event.key == pygame.K_f:
                            key_actions["f_clicked"] = False


                    elif event.type == pygame.MOUSEBUTTONDOWN:
                        if event.button == 1 and key_actions["split_clicked"]:
                            for node in self.all_network_nodes:
                                mouse_pos = pygame.mouse.get_pos()
                                scaled_mouse_x = (mouse_pos[0] + self.viewport.x) / self.zoom_factor
                                scaled_mouse_y = (mouse_pos[1] + self.viewport.y) / self.zoom_factor
                                mouse_pos = (scaled_mouse_x,scaled_mouse_y)
                                if node.rect.collidepoint(mouse_pos):
                                    self.split_node_based_on_highlighted_connections(node)
                                    self.update_edge_network()
                                    break

                        elif event.button == 1 and key_actions["shift_clicked"]:
                            mouse_pos = pygame.mouse.get_pos()
                            scaled_mouse_x = (mouse_pos[0] + self.viewport.x) / self.zoom_factor
                            scaled_mouse_y = (mouse_pos[1] + self.viewport.y) / self.zoom_factor
                            mouse_pos = (scaled_mouse_x, scaled_mouse_y)
                            for node in self.all_network_nodes:
                                if node.rect.collidepoint(mouse_pos):
                                    node.highlighted_by_mouse = not node.highlighted_by_mouse
                                    break

                        elif event.button == 1 and key_actions["v_clicked"]:
                            mouse_pos = pygame.mouse.get_pos()
                            scaled_mouse_x = (mouse_pos[0] + self.viewport.x) / self.zoom_factor
                            scaled_mouse_y = (mouse_pos[1] + self.viewport.y) / self.zoom_factor
                            mouse_pos = (scaled_mouse_x, scaled_mouse_y)
                            for node in self.all_network_nodes:
                                if node.rect.collidepoint(mouse_pos):
                                    node.visible_by_mouse = not node.visible_by_mouse
                                    self.queue.put(self.tk_application.create_listboxes)
                                    self.queue.put(self.tk_application.create_listboxes_visibility)
                                    break

                        elif event.button == 1 and key_actions["f_clicked"]:
                            key_actions["right_mouse_clicked"] = True
                            mouse_pos = pygame.mouse.get_pos()
                            scaled_mouse_x = (mouse_pos[0] + self.viewport.x) / self.zoom_factor
                            scaled_mouse_y = (mouse_pos[1] + self.viewport.y) / self.zoom_factor
                            mouse_pos = (scaled_mouse_x, scaled_mouse_y)
                            for node in self.all_network_nodes:
                                if node.rect.collidepoint(mouse_pos):
                                    node.fixed_by_mouse = not node.fixed_by_mouse
                                    break
                        elif event.button == 1 and not self.drawing_selection_rect_on_screen:
                            key_actions["left_mouse_clicked"] = True
                            mouse_pos = pygame.mouse.get_pos()
                            scaled_mouse_x = (mouse_pos[0] + self.viewport.x) / self.zoom_factor
                            scaled_mouse_y = (mouse_pos[1] + self.viewport.y) / self.zoom_factor
                            mouse_pos = (scaled_mouse_x, scaled_mouse_y)
                            for node in self.all_network_nodes:
                                if node.rect.collidepoint(mouse_pos):
                                    node.fixed_by_mouse = True
                                    selected_node = node
                                    break

                        elif event.button == 2:
                            key_actions["middle_mouse_clicked"] = True
                            mouse_pos = pygame.mouse.get_pos()
                            scaled_mouse_x = (mouse_pos[0] + self.viewport.x) / self.zoom_factor
                            scaled_mouse_y = (mouse_pos[1] + self.viewport.y) / self.zoom_factor
                            mouse_pos = (scaled_mouse_x, scaled_mouse_y)
                            for node in self.all_network_nodes:
                                if node.rect.collidepoint(mouse_pos):
                                    node.visible_by_mouse = not node.visible_by_mouse
                                    self.queue.put(self.tk_application.create_listboxes)
                                    self.queue.put(self.tk_application.create_listboxes_visibility)
                                    break
                        elif event.button == 3:
                            key_actions["right_mouse_clicked"] = True
                            self.drawing_selection_rect_on_screen = False
                            for node in self.all_network_nodes:
                                node.selected_by_drag_selection = False
                    elif event.type == pygame.MOUSEBUTTONUP:
                        if event.button == 1 and not key_actions["shift_clicked"]:
                            key_actions["left_mouse_clicked"] = False
                            selected_node = None
                            self.start_command_recording = True
                            if self.drawing_selection_rect_on_screen:
                                x_ = min(self.begin_selection_rect_pos[0], self.end_selection_rect_pos[0])
                                y_ = min(self.begin_selection_rect_pos[1], self.end_selection_rect_pos[1])
                                width_ = abs(self.end_selection_rect_pos[0] - self.begin_selection_rect_pos[0])
                                height_ = abs(self.end_selection_rect_pos[1] - self.begin_selection_rect_pos[1])
                                rect_to_draw = pygame.Rect(x_, y_, width_, height_)
                                for node in self.all_network_nodes:
                                    if rect_to_draw.colliderect(node.rect):
                                        node.selected_by_drag_selection = True
                            self.drawing_selection_rect_on_screen = False
                        elif event.button == 1 and key_actions["shift_clicked"]:
                             key_actions["left_mouse_clicked"] = False
                             selected_node = None
                             self.start_command_recording = True
                             if self.drawing_selection_rect_on_screen:
                                 x_ = min(self.begin_selection_rect_pos[0], self.end_selection_rect_pos[0])
                                 y_ = min(self.begin_selection_rect_pos[1], self.end_selection_rect_pos[1])
                                 width_ = abs(self.end_selection_rect_pos[0] - self.begin_selection_rect_pos[0])
                                 height_ = abs(self.end_selection_rect_pos[1] - self.begin_selection_rect_pos[1])
                                 rect_to_draw = pygame.Rect(x_, y_, width_, height_)
                                 for node in self.all_network_nodes:
                                     if rect_to_draw.colliderect(node.rect):
                                         node.highlighted_by_mouse = not node.highlighted_by_mouse
                             self.drawing_selection_rect_on_screen = False

                        elif event.button == 3:
                            key_actions["right_mouse_clicked"] = False
                if key_actions["up_clicked"]:
                    self.viewport.y -= 8
                elif key_actions["down_clicked"]:
                    self.viewport.y += 8
                if key_actions["right_clicked"]:
                    self.viewport.x += 8
                elif key_actions["left_clicked"]:
                    self.viewport.x -= 8
                if key_actions["left_mouse_clicked"]:

                    mouse_pos = pygame.mouse.get_pos()
                    scaled_mouse_x = (mouse_pos[0] + self.viewport.x) / self.zoom_factor
                    scaled_mouse_y = (mouse_pos[1] + self.viewport.y) / self.zoom_factor
                    scaled_mouse_pos = (scaled_mouse_x, scaled_mouse_y)
                    if selected_node is not None:
                        if not selected_node.selected_by_drag_selection:
                            if self.start_command_recording == True:
                                self.start_command_recording = False
                                self.previous_commands.append([(selected_node, scaled_mouse_pos[0], scaled_mouse_pos[1])])
                                self.next_commands_undone = []
                            selected_node.x = scaled_mouse_pos[0]
                            selected_node.y = scaled_mouse_pos[1]
                            if len(selected_node.split_of_nodes_partners)>0:
                                for node in self.all_network_nodes:
                                    if node in selected_node.split_of_nodes_partners:
                                        if node.rect.colliderect(selected_node.rect):
                                            selected_node= self.merge_split_nodes_back_together(selected_node, node)
                        else:
                            previous_x, previous_y = selected_node.x, selected_node.y
                            selected_node.x = scaled_mouse_pos[0]
                            selected_node.y = scaled_mouse_pos[1]
                            distance_x, distance_y = selected_node.x - previous_x, selected_node.y - previous_y
                            temp_commands = []
                            for node in self.all_network_nodes:
                                if node.selected_by_drag_selection and node != selected_node:
                                    if self.show_compartments:
                                        node.x += distance_x
                                        node.y += distance_y
                                        node.fixed_by_mouse = True
                                        temp_commands.append((node, node.x, node.y))
                                    elif node.node_type != "compartment":
                                        node.x += distance_x
                                        node.y += distance_y
                                        node.fixed_by_mouse = True
                                        temp_commands.append((node, node.x, node.y))


                            if self.start_command_recording == True:
                                temp_commands.append((selected_node, scaled_mouse_pos[0], scaled_mouse_pos[1]))
                                self.start_command_recording = False
                                self.previous_commands.append(temp_commands)
                                self.next_commands_undone = []

                    else:
                        if not self.drawing_selection_rect_on_screen:
                            mouse_pos = pygame.mouse.get_pos()
                            scaled_mouse_x = (mouse_pos[0] + self.viewport.x) / self.zoom_factor
                            scaled_mouse_y = (mouse_pos[1] + self.viewport.y) / self.zoom_factor
                            self.begin_selection_rect_pos = (scaled_mouse_x, scaled_mouse_y)
                            self.end_selection_rect_pos= self.begin_selection_rect_pos
                            self.drawing_selection_rect_on_screen = True
                            for node in self.all_network_nodes:
                                node.drawing_selection_rect = False
                        else:
                            mouse_pos = pygame.mouse.get_pos()
                            scaled_mouse_x = (mouse_pos[0] + self.viewport.x) / self.zoom_factor
                            scaled_mouse_y = (mouse_pos[1] + self.viewport.y) / self.zoom_factor
                            self.end_selection_rect_pos = (scaled_mouse_x, scaled_mouse_y)

                if self.drawing_selection_rect_on_screen:
                    x_ = min(self.begin_selection_rect_pos[0], self.end_selection_rect_pos[0])* self.zoom_factor - self.viewport.x
                    y_ = min(self.begin_selection_rect_pos[1], self.end_selection_rect_pos[1])* self.zoom_factor - self.viewport.y
                    width_ = abs(self.end_selection_rect_pos[0] - self.begin_selection_rect_pos[0]) * self.zoom_factor
                    height_ = abs(self.end_selection_rect_pos[1] - self.begin_selection_rect_pos[1]) *self.zoom_factor
                    rect_to_draw = pygame.Rect(x_, y_, width_, height_)
                    pygame.draw.rect(display, (255, 0, 0), rect_to_draw, 2)

                ## # Centre Dot
                # unscaled_middle_x = display.get_width() / 2
                # unscaled_middle_y = display.get_height() / 2
                # dot_size = 5  # Define the dot size
                # dot_color = (255, 165, 0)  # Orange color
                # pygame.draw.circle(display, dot_color, (int(unscaled_middle_x), int(unscaled_middle_y)), dot_size)

                # for idx, node in enumerate(self.all_network_nodes):
                #     if node.visible_by_mouse:
                #         if node.node_type == "reaction":
                #             pygame.draw.rect(display, (0, 200, 200), node.rect)
                #             if node.highlighted_by_mouse:
                #                 outer_rect = node.rect.inflate(10, 10)
                #                 outer_color = (255, 0, 0)
                #                 pygame.draw.rect(display, outer_color, outer_rect,2)  # Inn
                #             if self.show_reaction_names == "Reaction Names":
                #                 text_surface = self.font.render(node.instance.name, True, (255, 255, 255))
                #                 text_rect = text_surface.get_rect()
                #                 text_rect.center = (node.rect.center[0] , node.rect.center[1] - self.rect_size)
                #                 display.blit(text_surface, text_rect)
                #             elif self.show_reaction_names == "Reaction IDs":
                #                 text_surface = self.font.render(node.instance.id, True, (255, 255, 255))
                #                 text_rect = text_surface.get_rect()
                #                 text_rect.center = (node.rect.center[0] , node.rect.center[1] - self.rect_size)
                #                 display.blit(text_surface, text_rect)
                #             elif self.show_reaction_names == "No names":
                #                 pass
                #         elif node.node_type == "metabolite":
                #             if node.highlighted_by_mouse:
                #                 outer_rect = node.rect.inflate(10, 10)
                #                 outer_color = (255, 0, 0)
                #                 pygame.draw.rect(display, outer_color, outer_rect,2)  # Inn
                #             for idx, compartment in enumerate(self.compartments):
                #                 if node.instance.compartment == compartment.id:
                #                     pygame.draw.rect(display, (40 + (idx * 60) % 255, 100, 0), node.rect)
                #             if node.instance.fixed_node_input or node.instance.fixed_node_output:
                #                 if self.show_fixed_metabolite_names == "Fixed Metabolite Names":
                #                     text_surface = self.font.render(node.instance.name, True, (255, 255, 255))
                #                     text_rect = text_surface.get_rect()
                #                     text_rect.center = (node.rect.center[0] , node.rect.center[1] - self.rect_size)
                #                     display.blit(text_surface, text_rect)
                #                 elif self.show_fixed_metabolite_names == "Fixed Metabolite IDs":
                #                     text_surface = self.font.render(node.instance.id, True, (255, 255, 255))
                #                     text_rect = text_surface.get_rect()
                #                     text_rect.center = (node.rect.center[0] , node.rect.center[1] - self.rect_size)
                #                     display.blit(text_surface, text_rect)
                #                 elif self.show_fixed_metabolite_names == "No Fixed Names":
                #                     pass
                #             else:
                #                 if self.show_metabolite_names == "Metabolite Names":
                #                     text_surface = self.font.render(node.instance.name, True, (255, 255, 255))
                #                     text_rect = text_surface.get_rect()
                #                     text_rect.center = (node.rect.center[0] , node.rect.center[1] - self.rect_size)
                #                     display.blit(text_surface, text_rect)
                #                 elif self.show_metabolite_names == "Metabolite IDs":
                #                     text_surface = self.font.render(node.instance.id, True, (255, 255, 255))
                #                     text_rect = text_surface.get_rect()
                #                     text_rect.center = (node.rect.center[0] , node.rect.center[1] - self.rect_size)
                #                     display.blit(text_surface, text_rect)
                #                 elif self.show_metabolite_names == "No Names":
                #                     pass
                #         elif node.node_type == "compartment":
                #             for idx, compartment in enumerate(self.compartments):
                #                 if node.instance == compartment and self.show_compartments:
                #                     rect = node.rect
                #                     # Drawing lines for the cross inside the rectangle
                #                     pygame.draw.rect(display, (40 + (idx * 60) % 255, 100, 0), node.rect)
                #                     pygame.draw.line(display, (255, 255, 255), (rect.x, rect.y), (rect.x + rect.width, rect.y + rect.height), 3)
                #                     pygame.draw.line(display, (255, 255, 255), (rect.x + rect.width, rect.y), (rect.x, rect.y + rect.height), 3)

                #highlighted rects
                for node in self.all_network_nodes:
                    if node.visible_by_mouse and node.node_type != "boundary" and node.selected_by_drag_selection:
                        if self.show_compartments:
                            outer_rect = node.rect.inflate(15, 15)
                            scaled_outer_rect = pygame.Rect(outer_rect.x * self.zoom_factor - self.viewport.x,
                                                            outer_rect.y * self.zoom_factor - self.viewport.y,
                                                            outer_rect.width * self.zoom_factor,
                                                            outer_rect.height * self.zoom_factor)
                            outer_color = (0, 200, 0)
                            pygame.draw.rect(display, outer_color, scaled_outer_rect, 2)  # Inn

                        elif node.node_type != "compartment":
                            outer_rect = node.rect.inflate(15, 15)
                            scaled_outer_rect = pygame.Rect(outer_rect.x * self.zoom_factor - self.viewport.x,
                                                            outer_rect.y * self.zoom_factor - self.viewport.y,
                                                            outer_rect.width * self.zoom_factor,
                                                            outer_rect.height * self.zoom_factor)
                            outer_color = (0, 200, 0)
                            pygame.draw.rect(display, outer_color, scaled_outer_rect, 2)  # Inn

                for idx, node in enumerate(self.all_network_nodes):
                    if node.visible_by_mouse:
                        if node.node_type == "reaction":
                            scaled_rect = pygame.Rect(node.rect.x * self.zoom_factor - self.viewport.x,
                                                      node.rect.y * self.zoom_factor - self.viewport.y,
                                                      node.rect.width * self.zoom_factor,
                                                      node.rect.height * self.zoom_factor)
                            pygame.draw.rect(display, (0, 200, 200), scaled_rect)
                            # Other operations based on zoom...

                            if node.highlighted_by_mouse:
                                outer_rect = node.rect.inflate(10, 10)
                                scaled_outer_rect = pygame.Rect(outer_rect.x * self.zoom_factor - self.viewport.x,
                                                                outer_rect.y * self.zoom_factor - self.viewport.y,
                                                                outer_rect.width * self.zoom_factor,
                                                                outer_rect.height * self.zoom_factor)
                                outer_color = (255, 0, 0)
                                pygame.draw.rect(display, outer_color, scaled_outer_rect, 2)  # Inn

                            # Apply text rendering based on conditions and adjust positions with viewport
                            if self.show_reaction_names == "Reaction Names":
                                text_surface = self.font.render(node.instance.name, True, (255, 255, 255))
                                text_rect = text_surface.get_rect()
                                text_rect.center = (node.rect.center[0]  * self.zoom_factor- self.viewport.x,
                                                    (node.rect.center[1] * self.zoom_factor - self.rect_size*self.zoom_factor - self.viewport.y))
                                display.blit(text_surface, text_rect)
                            elif self.show_reaction_names == "Reaction IDs":
                                text_surface = self.font.render(node.instance.id, True, (255, 255, 255))
                                text_rect = text_surface.get_rect()
                                text_rect.center = (node.rect.center[0]  * self.zoom_factor- self.viewport.x,
                                                    (node.rect.center[1] * self.zoom_factor - self.rect_size*self.zoom_factor - self.viewport.y))
                                display.blit(text_surface, text_rect)
                            elif self.show_reaction_names == "No names":
                                pass

                        elif node.node_type == "metabolite":
                            # Similar adjustments for metabolites based on conditions and viewport
                            if node.highlighted_by_mouse:
                                outer_rect = node.rect.inflate(10, 10)
                                scaled_outer_rect = pygame.Rect(outer_rect.x * self.zoom_factor - self.viewport.x,
                                                                outer_rect.y * self.zoom_factor - self.viewport.y,
                                                                outer_rect.width * self.zoom_factor,
                                                                outer_rect.height * self.zoom_factor)
                                outer_color = (255, 0, 0)
                                pygame.draw.rect(display, outer_color, scaled_outer_rect, 2)  # Inn

                            for idx, compartment in enumerate(self.compartments):
                                if node.instance.compartment == compartment.id:
                                    scaled_compartment_rect = pygame.Rect(node.rect.x * self.zoom_factor - self.viewport.x,
                                                                          node.rect.y * self.zoom_factor - self.viewport.y,
                                                                          node.rect.width * self.zoom_factor,
                                                                          node.rect.height * self.zoom_factor)
                                    pygame.draw.rect(display, (40 + (idx * 60) % 255, 100, 0), scaled_compartment_rect)

                            # Apply text rendering based on conditions and adjust positions with viewport
                            if node.instance.fixed_node_input or node.instance.fixed_node_output:
                                if self.show_fixed_metabolite_names == "Fixed Metabolite Names":
                                    text_surface = self.font.render(node.instance.name, True, (255, 255, 255))
                                    text_rect = text_surface.get_rect()
                                    text_rect.center = (node.rect.center[0]  * self.zoom_factor- self.viewport.x,
                                                    (node.rect.center[1] * self.zoom_factor - self.rect_size*self.zoom_factor - self.viewport.y))
                                    display.blit(text_surface, text_rect)
                                elif self.show_fixed_metabolite_names == "Fixed Metabolite IDs":
                                    text_surface = self.font.render(node.instance.id, True, (255, 255, 255))
                                    text_rect = text_surface.get_rect()
                                    text_rect.center = (node.rect.center[0]  * self.zoom_factor- self.viewport.x,
                                                    (node.rect.center[1] * self.zoom_factor - self.rect_size*self.zoom_factor - self.viewport.y))
                                    display.blit(text_surface, text_rect)
                                elif self.show_fixed_metabolite_names == "No Fixed Names":
                                    pass
                            else:
                                # Apply text rendering based on conditions and adjust positions with viewport
                                if self.show_metabolite_names == "Metabolite Names":
                                    text_surface = self.font.render(node.instance.name, True, (255, 255, 255))
                                    text_rect = text_surface.get_rect()
                                    text_rect.center = (node.rect.center[0]  * self.zoom_factor- self.viewport.x,
                                                    (node.rect.center[1] * self.zoom_factor - self.rect_size*self.zoom_factor - self.viewport.y))
                                    display.blit(text_surface, text_rect)
                                elif self.show_metabolite_names == "Metabolite IDs":
                                    text_surface = self.font.render(node.instance.id, True, (255, 255, 255))
                                    text_rect = text_surface.get_rect()
                                    text_rect.center = (node.rect.center[0]  * self.zoom_factor- self.viewport.x,
                                                        (node.rect.center[1] * self.zoom_factor - self.rect_size*self.zoom_factor - self.viewport.y))
                                    display.blit(text_surface, text_rect)
                                elif self.show_metabolite_names == "No Names":
                                    pass
                        elif node.node_type == "compartment" and self.show_compartments:
                            for idx, compartment in enumerate(self.compartments):
                                if node.instance == compartment:
                                    rect = node.rect.copy()  # Make a copy of the node's rectangle

                                    # Adjust the rectangle's position based on the viewport
                                    rect.x -= self.viewport.x
                                    rect.y -= self.viewport.y

                                    # Adjust the rectangle's size based on the zoom factor
                                    rect.width *= self.zoom_factor
                                    rect.height *= self.zoom_factor

                                    # Drawing lines for the cross inside the rectangle with adjusted coordinates
                                    pygame.draw.rect(display, (40 + (idx * 60) % 255, 100, 0), rect)
                                    pygame.draw.line(display, (255, 255, 255), (rect.x, rect.y), (rect.x + rect.width, rect.y + rect.height), 3)
                                    pygame.draw.line(display, (255, 255, 255), (rect.x + rect.width, rect.y), (rect.x, rect.y + rect.height), 3)
            # compartment legends
                for idx, compartment in enumerate(self.compartments):
                    scaled_rect = pygame.Rect((width - 50) * self.zoom_factor - self.viewport.x,
                                              (50 + idx * 50) * self.zoom_factor - self.viewport.y,
                                              15 * self.zoom_factor,
                                              15 * self.zoom_factor)
                    pygame.draw.rect(display, (40 + (idx * 60) % 255, 100, 0), scaled_rect)

                    text_surface = self.font2.render(str(compartment.name), True, (50 + (idx * 50) % 255, 100, 0))
                    text_rect = text_surface.get_rect(center=((width - 100) * self.zoom_factor - self.viewport.x,
                                                              (58 + (idx * 50)) * self.zoom_factor - self.viewport.y))
                    display.blit(text_surface, text_rect)

                # Code for drawing lines, fluxes, and expressions with adjustments for the viewport
                if self.show_lines:
                    for idx, node in enumerate(self.all_network_nodes):
                        if node.visible_by_mouse:
                            for idx2 in range(1, len(self.edge_logic_dict["draw_lines_connections"][idx])):
                                if self.edge_logic_dict["draw_lines_connections"][idx2][idx] == 1:
                                    if self.all_network_nodes[idx2].visible_by_mouse:
                                        x_start = self.all_network_nodes[idx2].x
                                        y_start = self.all_network_nodes[idx2].y
                                        x_end = self.all_network_nodes[idx].x
                                        y_end = self.all_network_nodes[idx].y

                                        scaled_x_start = x_start * self.zoom_factor - self.viewport.x
                                        scaled_y_start = y_start * self.zoom_factor - self.viewport.y
                                        scaled_x_end = x_end * self.zoom_factor - self.viewport.x
                                        scaled_y_end = y_end * self.zoom_factor - self.viewport.y
                                        if self.show_color_lines:
                                            if self.all_network_nodes[idx].node_type == "metabolite":
                                                for idx3, compartment in enumerate(self.compartments):
                                                    if self.all_network_nodes[idx].instance.compartment == compartment.id:
                                                        color_ = (40 + (idx3 * 60) % 255, 100, 0)
                                                        if node.highlighted_by_mouse or self.all_network_nodes[idx2].highlighted_by_mouse :
                                                            color_ = (255,0,0)
                                                        self.draw_arrow_head((scaled_x_start, scaled_y_start), (scaled_x_end, scaled_y_end), display, color_)
                                                        break
                                            if self.all_network_nodes[idx2].node_type == "metabolite":
                                                for idx3, compartment in enumerate(self.compartments):
                                                    if self.all_network_nodes[idx2].instance.compartment == compartment.id:
                                                        color_ = (40 + (idx3 * 60) % 255, 100, 0)
                                                        if node.highlighted_by_mouse or self.all_network_nodes[idx2].highlighted_by_mouse:
                                                            color_ = (255,0,0)
                                                        self.draw_arrow_head((scaled_x_start, scaled_y_start), (scaled_x_end, scaled_y_end), display, color_)
                                                        break
                                        else:
                                            color_= (255,255,255)
                                            if node.highlighted_by_mouse or self.all_network_nodes[idx2].highlighted_by_mouse:
                                                color_ = (255, 0, 0)
                                            self.draw_arrow_head((scaled_x_start, scaled_y_start), (scaled_x_end, scaled_y_end), display, color_)

                if self.show_flux_expression == "Fluxes" or self.show_flux_expression == "Both flux and expression":
                    for node in self.all_network_nodes:
                        if node.node_type == "reaction":
                            try:
                                text_surface = self.font.render(str(self.fluxes_dict[node.instance.id]), True, (255, 0, 0))
                                text_rect = text_surface.get_rect()
                                text_rect.topright = ((node.rect.topright[0] + self.rect_size) * self.zoom_factor - self.viewport.x,
                                                      (node.rect.topright[1] + self.rect_size) * self.zoom_factor - self.viewport.y)
                                display.blit(text_surface, text_rect)
                            except Exception as e:
                                pass
                if self.show_flux_expression == "Expression" or self.show_flux_expression == "Both flux and expression":
                    for node in self.all_network_nodes:
                        if node.node_type == "reaction":
                            try:
                                text_surface = self.font.render(str(self.expression_dict[node.instance.id]), True, (0, 255, 0))
                                text_rect = text_surface.get_rect()
                                text_rect.topright = ((node.rect.topright[0] + self.rect_size) * self.zoom_factor - self.viewport.x,
                                                      (node.rect.topright[1] - self.rect_size-5) * self.zoom_factor - self.viewport.y)
                                display.blit(text_surface, text_rect)
                                if self.show_expressionless_reaction_highlight and self.expression_dict[node.instance.id] ==-1:
                                    inflated_rect = text_rect.inflate(5,5)
                                    pygame.draw.rect(display, (255,0,0), inflated_rect, width = 1)
                            except Exception as e:
                                pass
                if self.search_mode_enabled:
                    rendered_text = self.font3.render(self.search_input_text, True, (0, 0, 0))
                    text_input_rect.width = max(300, rendered_text.get_width() + 10)  # Adjust width based on rendered text width
                    text_input_rect.top = display.get_height() -400
                    text_input_rect.left = display.get_width() - 50-text_input_rect.width
                    pygame.draw.rect(display, (200, 200, 200), text_input_rect)  # Draw the input box
                    display.blit(rendered_text, (text_input_rect.x + 5, text_input_rect.y + 10))

                    for idx1, result in enumerate(self.search_suggestion_list):
                        rendered_recommendation = self.font3.render(result["text"], True, (255, 255, 255))
                        recommendation_rect = rendered_recommendation.get_rect()
                        recommendation_rect.topleft = (text_input_rect.left, text_input_rect.bottom + (idx1 * recommendation_rect.height+20))
                        if idx1 == self.current_search_index:
                            pygame.draw.rect(display, (150,150,150), recommendation_rect)
                        display.blit(rendered_recommendation, recommendation_rect.topleft)
                        start_pos = result["start_pos"]
                        end_pos = result["end_pos"]
                        if start_pos != -1 and end_pos != -1:
                            nudging_value = 13
                            highlight_rect = pygame.Rect(recommendation_rect.left + (start_pos*nudging_value),
                                                         recommendation_rect.top,
                                                         (end_pos*nudging_value) - (start_pos*nudging_value),
                                                         recommendation_rect.height)
                            pygame.draw.rect(display, (255,255,0), highlight_rect, 2)



                # for idx, compartment in enumerate(self.compartments):
                #     rect = pygame.Rect(width - 50, 50 + idx * 50, 15, 15)
                #     pygame.draw.rect(display, (40 + (idx * 60) % 255, 100, 0), rect)
                #     text_surface = self.font2.render(str(compartment.name), True, (50 + (idx * 50) % 255, 100, 0))
                #     text_rect = text_surface.get_rect(center=(width - 100, 58 + (idx * 50)))
                #     display.blit(text_surface, text_rect)
                #
                # if self.show_lines:
                #     for idx, node in enumerate(self.all_network_nodes):
                #         if node.visible_by_mouse:
                #             for idx2 in range(1, len(self.edge_logic_dict["draw_lines_connections"][idx])):
                #                 if self.edge_logic_dict["draw_lines_connections"][idx2][idx] == 1 :
                #                     if self.all_network_nodes[idx2].visible_by_mouse:
                #                         x_start = self.all_network_nodes[idx2].x
                #                         y_start = self.all_network_nodes[idx2].y
                #                         x_end = self.all_network_nodes[idx].x
                #                         y_end = self.all_network_nodes[idx].y
                #                         color_ = (255,255,255)
                #
                #                         # pygame.draw.line(display, (255, 255, 255), (x_start, y_start), (x_end, y_end))
                #                         self.draw_arrow_head((x_start, y_start), (x_end,y_end), display, color_)
                # if self.show_flux_expression == "Fluxes":
                #     for node in self.all_network_nodes:
                #         if node.node_type == "reaction":
                #             try:
                #                 text_surface = self.font.render(str(self.fluxes_dict[node.instance.id]), True, (255, 0, 0))
                #                 text_rect = text_surface.get_rect()
                #                 text_rect.topright = (node.rect.topright[0] +  self.rect_size, node.rect.topright[1] +self.rect_size)
                #                 display.blit(text_surface, text_rect)
                #             except:
                #                 pass
                # elif self.show_flux_expression == "Expression":
                #     for node in self.all_network_nodes:
                #         if node.node_type == "reaction":
                #             try:
                #                 text_surface = self.font.render(str(self.expression_dict[node.instance.id]), True, (0, 255, 0))
                #                 text_rect = text_surface.get_rect()
                #                 text_rect.topright = (node.rect.topright[0] +  self.rect_size, node.rect.topright[1] -self.rect_size)
                #                 display.blit(text_surface, text_rect)
                #             except:
                #                 pass
                # elif self.show_flux_expression == "Both flux and expression":
                #     for node in self.all_network_nodes:
                #         if node.node_type == "reaction":
                #             try:
                #                 text_surface = self.font.render(str(self.fluxes_dict[node.instance.id]), True, (255, 0,0))
                #                 text_rect = text_surface.get_rect()
                #                 text_rect.topright = (node.rect.topright[0] +  self.rect_size, node.rect.topright[1] +self.rect_size)
                #                 display.blit(text_surface, text_rect)
                #                 text_surface = self.font.render(str(self.expression_dict[node.instance.id]), True, (0, 255, 0))
                #                 text_rect = text_surface.get_rect()
                #                 text_rect.topright = (node.rect.topright[0] +  self.rect_size, node.rect.topright[1] -self.rect_size)
                #                 display.blit(text_surface, text_rect)
                #             except Exception as e:
                #                 pass
                pygame.display.flip()

            self.calculate_forces_through_network()
            self.move_nodes_based_on_forces()
            # counter += 1
            # if counter >= iterations:
            #     running = False
        if pygame_on:
            pygame.quit()

    def draw_arrow_head(self,end,start, display, color_):
        dx = end[0] - start[0]
        dy = end[1] - start[1]
        angle = math.atan2(dy, dx)
        arrow_len = max(3, self.rect_size-6)
        arrow_angle = math.pi / 6
        arrow1 = (end[0] - arrow_len * math.cos(angle - arrow_angle),
                  end[1] - arrow_len * math.sin(angle - arrow_angle))
        arrow2 = (end[0] - arrow_len * math.cos(angle + arrow_angle),
                  end[1] - arrow_len * math.sin(angle + arrow_angle))

        pygame.draw.line(display, color_, (start[0], start[1]), (end[0], end[1]))
        pygame.draw.polygon(display,color_ , (end, arrow1, arrow2))

    def pygame_draw_network_code(self, display):
        pass

    def pygame_event_handling_code(self):
        pass

    def undo_last_action(self):
        if len(self.previous_commands) >0:
            command = self.previous_commands.pop()
            temp_commands = []
            for command_tuple in command:
                temp_commands.append((command_tuple[0], command_tuple[0].x, command_tuple[0].y))
                command_tuple[0].x = command_tuple[1]
                command_tuple[0].y = command_tuple[2]

            self.next_commands_undone.append(temp_commands)

    def redo_last_action(self):
        if len(self.next_commands_undone)>0:
            temp_commands = []
            command = self.next_commands_undone.pop()
            for command_tuple in command:
                temp_commands.append((command_tuple[0], command_tuple[0].x, command_tuple[0].y))
                command_tuple[0].x = command_tuple[1]
                command_tuple[0].y = command_tuple[2]
            self.previous_commands.append(temp_commands)

    def add_currently_selected_to_control_group(self, ctrl_number):
        self.control_groups_for_selection[ctrl_number] = []
        for node in self.all_network_nodes:
            if node.selected_by_drag_selection:
                self.control_groups_for_selection[ctrl_number].append(node)

    def select_control_group(self, ctrl_number):
        for node in self.all_network_nodes:
            node.selected_by_drag_selection = False

        for node in self.control_groups_for_selection[ctrl_number]:
            node.selected_by_drag_selection = True


    def update_search_recommendations(self):
        self.search_suggestion_list = []
        self.full_search_list =[]
        for node in self.all_network_nodes:
            if node.node_type == "compartment" and self.show_compartments:
                names_string = node.instance.name + "  " + node.instance.id
                self.full_search_list.append((names_string, node))
            elif node.node_type != "boundary" and  node.visible_by_mouse:
                names_string = node.instance.name + "  " + node.instance.id
                self.full_search_list.append((names_string, node))
        for recommendation in self.full_search_list:
            start_pos = recommendation[0].lower().find(self.search_input_text.lower())
            if start_pos != -1:
                end_pos = start_pos + len(self.search_input_text)
                self.search_suggestion_list.append({
                    "text": recommendation[0],
                    "start_pos": start_pos,
                    "end_pos": end_pos,
                    "node": recommendation[1]
                })


    def calculate_distance_matrix(self):
        coordinates = np.array([(obj.x, obj.y) for obj in self.all_network_nodes])

        diff = coordinates[:, np.newaxis, :] - coordinates[np.newaxis, :, :]
        distance_matrix = np.linalg.norm(diff, axis=2)
        angle_matrix = np.arctan2(diff[:, :, 1], diff[:, :, 0])
        return distance_matrix, angle_matrix

    def move_nodes_based_on_forces(self):
        SMALL_CONSTANT = 0.000015
        for node in self.all_network_nodes:
            if node.fixed_by_mouse:
                pass
            elif node.fixed:
                node.x += node.x_vector * SMALL_CONSTANT
            else:
                node.x += node.x_vector * SMALL_CONSTANT
                node.x = node.x % (width-50)
                node.y += node.y_vector * SMALL_CONSTANT
                node.y = node.y % (height-50)
            if node.x > width:
                node.x = width
            elif node.x < 0:
                node.x = 0
            if node.y > width:
                node.y = width
            elif node.y < 0:
                node.y = 0

            node.rect.update(node.x-(self.rect_size)/2, node.y-(self.rect_size)/2, self.rect_size, self.rect_size)

    def calculate_forces_through_network(self):
        distance_matrix, angle_matrix = self.calculate_distance_matrix()
        distance_matrix = distance_matrix + 1
        cos_angles = np.cos(angle_matrix)
        sin_angles = np.sin(angle_matrix)

        # # General repulsion
        np.fill_diagonal(distance_matrix, 1000000)
        force_magnitudes = (200 / np.power(np.sqrt(distance_matrix), 2))
        force_magnitudes = np.where(force_magnitudes < 0.0000000001, 0, force_magnitudes)
        x_forces = np.sum(force_magnitudes * cos_angles * self.normal_repulsion_value * self.edge_logic_dict["normal_repulsion"], axis=1)
        y_forces = np.sum(force_magnitudes * sin_angles * self.normal_repulsion_value * self.edge_logic_dict["normal_repulsion"], axis=1)

        # Compartment repulsion
        force_magnitudes = 300 / np.power(np.sqrt(distance_matrix), 1)
        force_magnitudes = np.where(force_magnitudes < 0.0000001, 0, force_magnitudes)
        x_forces = x_forces + np.sum(force_magnitudes * cos_angles * self.comp_comp_repulsion_value * self.edge_logic_dict["compartment_compartment_repulsion"],
                                     axis=1)
        y_forces = y_forces + np.sum(force_magnitudes * sin_angles * self.comp_comp_repulsion_value * self.edge_logic_dict["compartment_compartment_repulsion"],
                                     axis=1)

        # Fixed metabolite repulsion
        force_magnitudes = 3000 / np.power(np.sqrt(distance_matrix), 2)
        force_magnitudes = np.where(force_magnitudes < 0.0000001, 0, force_magnitudes)
        x_forces = x_forces + np.sum(force_magnitudes * cos_angles * self.fixed_metabolite_repulsion_value * self.edge_logic_dict["fixed_metabolite_repulsion"],
                                     axis=1)

        # Boundary Repulsion
        force_magnitudes = 500000 / np.power(np.sqrt(distance_matrix), 5)
        force_magnitudes = np.where(force_magnitudes < 0.0000001, 0, force_magnitudes)
        x_forces = x_forces + np.sum(force_magnitudes * cos_angles * self.boundary_repulsion_value * self.edge_logic_dict["boundary_repulsion"], axis=1)
        y_forces = y_forces + np.sum(force_magnitudes * sin_angles * self.boundary_repulsion_value * self.edge_logic_dict["boundary_repulsion"], axis=1)

        # rxn_met attraction
        np.fill_diagonal(distance_matrix, 0)
        force_magnitudes = -np.power(distance_matrix, 1) / 10
        force_magnitudes = np.where(force_magnitudes > -0.0001, 0, force_magnitudes)
        x_forces = x_forces + np.sum(force_magnitudes * cos_angles * self.rxn_met_attraction_value * self.edge_logic_dict["reaction_metabolite_attraction"],
                                     axis=1)
        y_forces = y_forces + np.sum(force_magnitudes * sin_angles * self.rxn_met_attraction_value * self.edge_logic_dict["reaction_metabolite_attraction"],
                                     axis=1)

        # met_comp attraction
        force_magnitudes = -np.power(distance_matrix, 1) / 20
        force_magnitudes = np.where(force_magnitudes > -0.0001, 0, force_magnitudes)
        x_forces = x_forces + np.sum(
            force_magnitudes * cos_angles * self.metabolite_compartment_attraction_value * self.edge_logic_dict["metabolite_compartment_attraction"], axis=1)
        y_forces = y_forces + np.sum(
            force_magnitudes * sin_angles * self.metabolite_compartment_attraction_value * self.edge_logic_dict["metabolite_compartment_attraction"], axis=1)

        for node, x_force, y_force in zip(self.all_network_nodes, x_forces, y_forces):
            node.x_vector = x_force
            node.y_vector = y_force

    def split_node_based_on_highlighted_connections(self, node):
        connected_output_instances = []
        connected_input_instances = []
        for instance in node.output_instances:
            connected_output_instances.append(instance)
        for instance in node.input_instances:
            connected_input_instances.append(instance)
        to_split_off_input_instances =[]
        to_split_off_output_instances =[]
        to_split_off_output_nodes = []
        to_split_off_input_nodes = []
        for other_node in self.all_network_nodes:
            if other_node.instance in connected_output_instances and other_node.highlighted_by_mouse:
                to_split_off_output_instances.append(other_node.instance)

                to_split_off_output_nodes.append(other_node)
                try:
                    node.output_instances.remove(other_node.instance)
                except:
                    pass
            elif other_node.instance in connected_input_instances and other_node.highlighted_by_mouse:
                to_split_off_input_instances.append(other_node.instance)
                to_split_off_input_nodes.append(other_node)
                try:
                    node.input_instances.remove(other_node.instance)
                except:
                    pass
        #remove other nodes from node
        #remove node from alternative nodes that connect to current node
        current_node_instance = node.instance

        # create new node with connected nodes and similiar id ?
        if len(to_split_off_output_instances)>0 or len(to_split_off_input_instances)>0:
            id_number = len(self.all_network_nodes)
            new_current_node_instance = copy.deepcopy(current_node_instance)

            self.all_network_nodes.append(NetworkNodesForPhysicsSimulation(node.x + 20, node.y+20, node.instance,id_number,self ))
            self.all_network_nodes[-1].split_of_nodes_partners = node.split_of_nodes_partners.copy()
            self.all_network_nodes[-1].split_of_nodes_partners.append(node)

            for alt_node in node.split_of_nodes_partners:
                alt_node.split_of_nodes_partners.append(self.all_network_nodes[-1])
            node.split_of_nodes_partners.append(self.all_network_nodes[-1])
            current_node_instance.id = current_node_instance.id + "_" + str(len(self.all_network_nodes[-1].split_of_nodes_partners))
            current_node_instance.name = current_node_instance.name + "_" + str(len(self.all_network_nodes[-1].split_of_nodes_partners))
            new_current_node_instance.id = new_current_node_instance.id + "_" + str(len(self.all_network_nodes[-1].split_of_nodes_partners)+1)
            new_current_node_instance.name = new_current_node_instance.name + "_" + str(len(self.all_network_nodes[-1].split_of_nodes_partners)+1)
            self.all_network_nodes[-1].instance = new_current_node_instance
            self.all_network_nodes[-1].input_instances = to_split_off_input_instances
            self.all_network_nodes[-1].output_instances = to_split_off_output_instances
            self.all_network_nodes[-1].fixed_by_mouse = node.fixed_by_mouse

            for other_instance, other_node in zip(to_split_off_output_instances, to_split_off_output_nodes):
                try:
                    other_node.input_instances.append(self.all_network_nodes[-1].instance)
                    other_node.input_instances.remove(current_node_instance)
                except:
                    pass
            for other_instance, other_node in zip(to_split_off_input_instances, to_split_off_input_nodes):
                try:
                    other_node.output_instances.append(self.all_network_nodes[-1].instance)
                    other_node.output_instances.remove(current_node_instance)
                except:
                    pass

    def merge_split_nodes_back_together(self, node1, node2):
        match1 = re.search(r'_(\d+)', node1.instance.id)
        result1 = int(match1.group(1)) if match1 else 0
        match2 = re.search(r'_(\d+)', node2.instance.id)
        result2 = int(match2.group(1)) if match2 else 0
        highest_node, lowest_node = (node1, node2) if result1 > result2 else (node2, node1)
        to_change_input_instances =highest_node.input_instances
        to_change_output_instances =highest_node.output_instances
        lowest_node.input_instances.extend(to_change_input_instances)
        lowest_node.output_instances.extend(to_change_output_instances)

        for instance in to_change_output_instances:
            for node in self.all_network_nodes:
                if node.instance == instance:
                    node.input_instances.remove(highest_node.instance)
                    node.input_instances.append(lowest_node.instance)
                    break
        for instance in to_change_input_instances:
            for node in self.all_network_nodes:
                if node.instance ==instance:
                    node.output_instances.remove(highest_node.instance)
                    node.output_instances.append(lowest_node.instance)

        lowest_node.split_of_nodes_partners.remove(highest_node)
        if len(lowest_node.split_of_nodes_partners)==0:
            match = re.match(r'^(.*?[a-zA-Z])(?:_\d+)?$', lowest_node.instance.id)
            if match:
                lowest_node.instance.id = match.group(1)
            match = re.match(r'^(.*?[a-zA-Z\s()-]+?)(?:_\d+)?$|^(.*?\))_\d+$', lowest_node.instance.name)
            if match:
                lowest_node.instance.name = match.group(1)
        else:
            for partners in lowest_node.split_of_nodes_partners:
                partners.split_of_nodes_partners.remove(highest_node)

        self.all_network_nodes.remove(highest_node)
        for idx,node in enumerate(self.all_network_nodes):
            node.id_number = idx
        remaining_node = lowest_node
        self.update_edge_network()
        return remaining_node
    def extract_non_pseudo_nodes_and_edges(self):
        # remove
        nodes, edges = 0, 0
        # remove
        return nodes, edges
        pass



class Application:
    def __init__(self, root, width, height):
        self.root = root
        self.root.title("JSON Viewer")
        self.WIDTH = width
        self.HEIGHT = height

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

        self.selected_option = tk.StringVar()
        self.get_user_input()

    def get_user_input(self):
        if len(sys.argv) > 1:
            filename = ' '.join(sys.argv[1:])
            full_path = os.path.abspath(filename)
            self.entry_var.set(full_path)
            self.upload_file()

            return False

        else:
            # Not running in an interactive terminal
             return True

    def check_queue(self):
        while True:
            try:
                action = self.met_network.queue.get(block=False)
                self.root.after(0, action)
            except queue.Empty:
                break
        self.root.after(100, self.check_queue)  # Repeat periodically

    def browse_file(self):
        file_path = filedialog.askopenfilename(filetypes=[("JSON files", "*.json")])
        self.entry_var.set(file_path)

    def key_interactions(self, event):
        if event.keycode == 27:
            self.root.destroy()
            try:
                pygame.quit()
            except:
                pass

    def upload_file(self):
        file_path = self.entry_var.get()
        try:
            with open(file_path, 'r') as file:
                # Read JSON data
                json_data = json.load(file)

                # Open a new window for JSON data visualization
                self.parse_json_metabolic_network_data(json_data)
                self.met_network = MetabolicNetwork(self.list_nodes_metabolites, self.list_edges_reactions, self)
                # self.upload_button.destroy()
                # self.browse_button.destroy()
                self.create_GUI_after_loading()
        except FileNotFoundError:
            print("File not found.")

    def create_GUI_after_loading(self):
        self.root.destroy()
        self.root = tk.Tk()
        self.create_listboxes()
        self.create_slider()
        self.create_slider_normal_repulsion()
        self.create_slider_comp_comp_repulsion()
        self.create_slider_boundary_repulsion()
        self.create_slider_fixed_metabolite_slider()
        self.create_slider_rxn_met_attraction_slider()
        self.create_slider_met_comp_attraction_slider()
        self.create_load_network_button()
        self.create_stop_simulation_button()
        self.create_freeze_all_nodes_button()
        self.create_unfreeze_all_nodes_button()
        self.create_font_size_slider()
        self.create_rect_size_slider()
        self.root.bind("<KeyPress>", self.key_interactions)
        self.create_reaction_name_checkbutton()
        self.create_metabolite_name_checkbutton()
        self.create_fixed_metabolite_name_checkbutton()
        self.root.update_idletasks()
        self.create_show_unshow_lines_button()
        self.create_show_unshow_compartments_button()
        self.create_browse_and_upload_excel_file()
        self.create_show_flux_expression_checkbutton()
        self.create_update_visible_button()
        self.create_restart_latest_simulation()
        self.create_save_latest_simulation()
        self.create_load_previously_made_model()
        self.create_color_uncolor_lines_button()
        self.create_highlight_expressionless_button()
        self.create_canvas_settings()
        width_ = self.root.winfo_width()
        height_ = self.root.winfo_height()
        self.root.geometry(f"{width_ + 270}x{height_+50}+0+0")
        self.check_queue()

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
        self.load_network_button = tk.Button(self.root, text="Start Simulation", command=self.load_network)
        self.load_network_button.grid(row=3, column=2, padx=10, pady=5, sticky = "ew")

    def load_network(self):
        self.met_network.calculate_network()
        # self.create_and_draw_network_canvas()
        self.create_listboxes_visibility()

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

        self.canvas = tk.Canvas(self.network_window, width=self.WIDTH, height=self.HEIGHT)
        self.draw_network(self.canvas)
        self.canvas.place(in_=self.network_window, x=0, y=0)

    def create_listboxes(self):
        label1 = tk.Label(self.root, text="included")
        label1.grid(row=0, column=0, padx=10, pady=5)
        self.left_list = tk.Listbox(self.root, width=20, height=20)
        self.left_list.grid(row=1, column=0, padx=10, pady=5)
        label2 = tk.Label(self.root, text="not included")
        label2.grid(row=0, column=1, padx=10, pady=5)
        self.right_list = tk.Listbox(self.root, width=20, height=20)
        self.right_list.grid(row=1, column=1, padx=10, pady=5)

        self.names_left = [name.name for name in self.met_network.part_of_physics_metabolites]
        self.names_right = [name.name for name in self.met_network.not_part_of_physics_metabolites]

        for item in self.names_left:
            self.left_list.insert(tk.END, item)
        for item in self.names_right:
            self.right_list.insert(tk.END, item)

        self.left_list.bind('<Double-Button-1>', lambda event: self.on_double_click(event, self.left_list, self.right_list))
        self.right_list.bind('<Double-Button-1>', lambda event: self.on_double_click(event, self.right_list, self.left_list))

    def move_item(self, event, from_list, to_list):
        try:
            selected_item = from_list.get(from_list.curselection())
            to_list.insert(tk.END, selected_item)
            from_list.delete(from_list.curselection())
        except:
            pass
    def on_double_click(self, event, from_list, to_list):
        self.move_item(event, from_list, to_list)
        self.met_network.part_of_physics_metabolites = [self.met_network.find_metabolite(met) for met in self.left_list.get("0", tk.END)]
        self.met_network.not_part_of_physics_metabolites = [self.met_network.find_metabolite(met) for met in self.right_list.get("0", tk.END)]
        try:
            self.met_network.update_edge_network()
        except Exception as e:
            print(e)

    def update_slider(self, value):
        value = int(value)
        self.met_network.slider_value = value
        self.right_list.destroy()
        self.left_list.destroy()
        self.met_network.update_included_not_included_based_on_slider()
        self.create_listboxes()
        try:
            self.met_network.update_edge_network()
        except:
            pass

    def create_slider(self):
        label = tk.Label(self.root, text="Metabolites with equal or more connections are excluded")
        label.grid(row=2, column=0, columnspan=2, padx=10, pady=5, sticky="ew")

        scale = tk.Scale(self.root, from_=self.met_network.slider_start, to=self.met_network.slider_end, orient=tk.HORIZONTAL, command=self.update_slider)
        scale.grid(row=3, column=0, columnspan=2, padx=10, pady=5, sticky="ew")
        scale.set(self.met_network.slider_value)

    def create_slider_normal_repulsion(self):
        # self.root.grid_rowconfigure(1, weight=1, maxsize=50)  # Adjust minsize as needed
        label = tk.Label(self.root, text="Normal Repulsion")
        label.grid(row=2, column=3, columnspan=2, padx=10, pady=5, sticky="ew")

        scale = tk.Scale(self.root, from_=1, to=5000, orient=tk.HORIZONTAL, command=self.update_normal_repulsion_slider)
        scale.grid(row=3, column=3, columnspan=2, padx=10, pady=5, sticky="ew")
        scale.set(self.met_network.normal_repulsion_value)

    def update_normal_repulsion_slider(self, value):
        self.met_network.normal_repulsion_value = int(value)

    def create_slider_comp_comp_repulsion(self):
        # self.root.grid_rowconfigure(1, weight=1, maxsize=50)  # Adjust minsize as needed
        label = tk.Label(self.root, text="Intercompartment Repulsion")
        label.grid(row=4, column=3, columnspan=2, padx=10, pady=5, sticky="ew")

        scale = tk.Scale(self.root, from_=0, to=50000, orient=tk.HORIZONTAL, command=self.update_comp_comp_repulsion_slider)
        scale.grid(row=5, column=3, columnspan=2, padx=10, pady=5, sticky="ew")
        scale.set(self.met_network.comp_comp_repulsion_value)

    def update_comp_comp_repulsion_slider(self, value):
        self.met_network.comp_comp_repulsion_value = int(value)

    def create_slider_boundary_repulsion(self):
        # self.root.grid_rowconfigure(1, weight=1, maxsize=50)  # Adjust minsize as needed
        label = tk.Label(self.root, text="Boundary Repulsion")
        label.grid(row=6, column=3, columnspan=2, padx=10, pady=5, sticky="ew")

        scale = tk.Scale(self.root, from_=0, to=50000, orient=tk.HORIZONTAL, command=self.update_boundary_repulsion_slider)
        scale.grid(row=7, column=3, columnspan=2, padx=10, pady=5, sticky="ew")
        scale.set(self.met_network.boundary_repulsion_value)

    def update_boundary_repulsion_slider(self, value):
        self.met_network.boundary_repulsion_value = int(value)

    def create_slider_rxn_met_attraction_slider(self):
        # self.root.grid_rowconfigure(1, weight=1, maxsize=50)  # Adjust minsize as needed
        label = tk.Label(self.root, text="Reaction-Metabolite Attraction")
        label.grid(row=10, column=3, columnspan=2, padx=10, pady=5, sticky="ew")

        scale = tk.Scale(self.root, from_=0, to=50000, orient=tk.HORIZONTAL, command=self.update_rxn_met_attraction_slider)
        scale.grid(row=11, column=3, columnspan=2, padx=10, pady=5, sticky="ew")
        scale.set(self.met_network.rxn_met_attraction_value)

    def update_rxn_met_attraction_slider(self, value):
        self.met_network.rxn_met_attraction_value = int(value)
        # self.met_network.stop_simulation()
        # self.met_network.start_simulation()

    def create_slider_met_comp_attraction_slider(self):
        # self.root.grid_rowconfigure(1, weight=1, maxsize=50)  # Adjust minsize as needed
        label = tk.Label(self.root, text="Metabolite-Compartment Attraction")
        label.grid(row=12, column=3, columnspan=2, padx=10, pady=5, sticky="ew")

        scale = tk.Scale(self.root, from_=0, to=50000, orient=tk.HORIZONTAL, command=self.update_met_comp_slider)
        scale.grid(row=13, column=3, columnspan=2, padx=10, pady=5, sticky="ew")
        scale.set(self.met_network.metabolite_compartment_attraction_value)

    def update_met_comp_slider(self, value):
        self.met_network.metabolite_compartment_attraction_value = int(value)

    def create_slider_fixed_metabolite_slider(self):
        # self.root.grid_rowconfigure(1, weight=1, maxsize=50)  # Adjust minsize as needed
        label = tk.Label(self.root, text="Input_metabolite-Output_metabolite Repulsion")
        label.grid(row=8, column=3, columnspan=2, padx=10, pady=5, sticky="ew")

        scale = tk.Scale(self.root, from_=0, to=50000, orient=tk.HORIZONTAL, command=self.update_fixed_metabolite_repulsion_slider)
        scale.grid(row=9, column=3, columnspan=2, padx=10, pady=5, sticky="ew")
        scale.set(self.met_network.fixed_metabolite_repulsion_value)

    def update_fixed_metabolite_repulsion_slider(self, value):
        self.met_network.fixed_metabolite_repulsion_value = int(value)

    def create_stop_simulation_button(self):
        button = tk.Button(self.root, text="Stop simulation", command=self.on_stop_simulation_button_click)
        button.grid(row=4, column=2, padx=5, pady=5, sticky="ew")

    def on_stop_simulation_button_click(self):
        self.stop_simulation()

    def stop_simulation(self):
        global running
        running = False

    def create_freeze_all_nodes_button(self):
        button = tk.Button(self.root, text="Freeze all nodes", command=self.on_freeze_all_nodes_button_click)
        button.grid(row=5, column=1, padx=5, pady=5, sticky="ew")

    def on_freeze_all_nodes_button_click(self):
        for node in self.met_network.all_network_nodes:
            node.fixed_by_mouse = True

    def create_unfreeze_all_nodes_button(self):
        button = tk.Button(self.root, text="Unfreeze all nodes", command=self.on_unfreeze_all_nodes_button_click)
        button.grid(row=6, column=1, padx=5, pady=5, sticky="ew")

    def on_unfreeze_all_nodes_button_click(self):
        for node in self.met_network.all_network_nodes:
            node.fixed_by_mouse = False

    def create_font_size_slider(self):
        label = tk.Label(self.root, text="Fontsize")
        label.grid(row=8, column=0, columnspan=2, padx=10, pady=5, sticky="ew")
        slider = tk.Scale(self.root, from_=1, orient=tk.HORIZONTAL, to=50, command=self.on_font_size_slider)
        slider.grid(row=9, column=0, columnspan=2, padx=10, pady=5, sticky="ew")
        slider.set(self.met_network.font_size)

    def on_font_size_slider(self, value):
        try:
            self.met_network.font = pygame.font.SysFont(None, int(value))
        except:
            pass

    def create_rect_size_slider(self):
        label = tk.Label(self.root, text="Rectangle size")
        label.grid(row=10, column=0, columnspan=2, padx=10, pady=5, sticky="ew")
        slider = tk.Scale(self.root, from_=1, to=50, orient=tk.HORIZONTAL, command=self.on_rect_size_slider)
        slider.grid(row=11, column=0, columnspan=2, padx=10, pady=5, sticky="ew")
        slider.set(self.met_network.font_size)

    def on_rect_size_slider(self, value):
        self.met_network.rect_size = int(value)

    def create_line_width_slider(self):
        label = tk.Label(self.root, text="Line Width")
        label.grid(row=12, column=0, columnspan=2, padx=10, pady=5, sticky="ew")
        slider = tk.Scale(self.root, from_=1, to=50, orient=tk.HORIZONTAL, command=self.on_line_width_slider)
        slider.grid(row=13, column=0, columnspan=2, padx=10, pady=5, sticky="ew")
        slider.set(self.met_network.font_size)

    def on_line_width_slider(self, value):
        self.met_network.line_width = int(value)

    def create_reaction_name_checkbutton(self):
        self.selected_option = tk.StringVar(self.root, "Reaction Names")
        option_frame = tk.LabelFrame(self.root, text="Show rxn names options")
        option_frame.grid(row=15, column=0, padx=10, pady=5, sticky="ew")
        options = ["No Names", "Reaction Names", "Reaction IDs"]
        for idx, option in enumerate(options):
            tk.Radiobutton(option_frame, text=option, variable=self.selected_option, value=option, command=self.on_reaction_name_checkbutton).grid(row=idx,
                                                                                                                                                   column=0,
                                                                                                                                                   sticky=tk.W)
    def on_reaction_name_checkbutton(self):
        option = self.selected_option.get()
        self.met_network.show_reaction_names = option

    def create_metabolite_name_checkbutton(self):
        self.selected_option_metabolite = tk.StringVar(self.root, "No Names")
        option_frame = tk.LabelFrame(self.root, text="Show met names options")
        option_frame.grid(row=15, column=1, padx=10, pady=5, sticky="ew")
        options = ["No Names", "Metabolite Names", "Metabolite IDs"]
        for idx, option in enumerate(options):
            tk.Radiobutton(option_frame, text=option, variable=self.selected_option_metabolite, value=option, command=self.on_metabolite_name_checkbutton).grid(
                row=idx, column=1, sticky=tk.W)

    def on_metabolite_name_checkbutton(self):
        option = self.selected_option_metabolite.get()
        self.met_network.show_metabolite_names = option

    def create_fixed_metabolite_name_checkbutton(self):
        self.selected_option_fixed_metabolite = tk.StringVar(self.root, "No Fixed Names")
        option_frame = tk.LabelFrame(self.root, text="Show fixed met names options")
        option_frame.grid(row=16, column=0, padx=10, pady=5, sticky="ew")
        options = ["No Fixed Names", "Fixed Metabolite Names", "Fixed Metabolite IDs"]
        for idx, option in enumerate(options):
            tk.Radiobutton(option_frame, text=option, variable=self.selected_option_fixed_metabolite, value=option,
                           command=self.on_fixed_metabolite_name_checkbutton).grid(row=idx, column=0, sticky=tk.W)

    def on_fixed_metabolite_name_checkbutton(self):
        option = self.selected_option_fixed_metabolite.get()
        self.met_network.show_fixed_metabolite_names = option

    def create_show_flux_expression_checkbutton(self):
        self.selected_option_flux_expression = tk.StringVar(self.root, "None")
        option_frame = tk.LabelFrame(self.root, text="Show Excel flux/expression")
        option_frame.grid(row=16, column=1, padx=10, pady=5, sticky="ew")
        options = ["None", "Fluxes", "Expression", "Both flux and expression"]
        for idx, option in enumerate(options):
            tk.Radiobutton(option_frame, text=option, variable=self.selected_option_flux_expression, value=option,
                           command=self.on_show_flux_expression_checkbutton).grid(row=idx, column=0, sticky=tk.W)

    def on_show_flux_expression_checkbutton(self):
        option = self.selected_option_flux_expression.get()
        self.met_network.show_flux_expression = option

    def create_show_unshow_lines_button(self):
        button = tk.Button(self.root, text="Show/unshow lines toggle", command=self.on_show_unshow_lines_button_click)
        button.grid(row=4, column=0, padx=5, pady=5, sticky="ew")

    def on_show_unshow_lines_button_click(self):
        self.met_network.show_lines = not self.met_network.show_lines

    def create_color_uncolor_lines_button(self):
        button = tk.Button(self.root, text="Color/uncolor lines toggle", command=self.on_color_uncolor_lines_button_click)
        button.grid(row=6, column=0, padx=5, pady=5, sticky="ew")

    def on_color_uncolor_lines_button_click(self):
        self.met_network.show_color_lines = not self.met_network.show_color_lines

    def create_highlight_expressionless_button(self):
        button = tk.Button(self.root, text="Highlight expressionless reactions expression", command=self.on_highlight_expressionless_button_click)
        button.grid(row=7, column=0, padx=5, pady=5, sticky="ew")

    def on_highlight_expressionless_button_click(self):
        self.met_network.show_expressionless_reaction_highlight = not self.met_network.show_expressionless_reaction_highlight
    def create_restart_latest_simulation(self):
        button = tk.Button(self.root, text="Restart the latest simulation", command=self.on_restart_latest_simulation_click)
        button.grid(row=6, column=2, padx=5, pady=5, sticky="ew")

    def on_restart_latest_simulation_click(self):
        self.met_network.restart_simulation()

    def create_save_latest_simulation(self):
        button = tk.Button(self.root, text="Save latest simulation as .json", command=self.on_save_latest_simulation_click)
        button.grid(row=7, column=2, padx=5, pady=5, sticky="ew")

    def on_save_latest_simulation_click(self):
        self.met_network.save_as_json()
        file_path = filedialog.asksaveasfilename(defaultextension=".json" , filetypes=[("JSON files", "*.json")])

        with open(file_path, "w") as file:
            json.dump(self.met_network.list_with_dicts_to_save, file, indent =4 )

    def create_load_previously_made_model(self):
        button = tk.Button(self.root, text="Load previously made JSON model", command=self.on_load_previously_made_model_click)
        button.grid(row=8, column=2, padx=5, pady=5, sticky="ew")

    def on_load_previously_made_model_click(self):
        file_path = filedialog.askopenfilename(filetypes=[("JSON files", "*.json")])
        with open(file_path, 'r') as file:
            data = json.load(file)
            self.met_network.read_from_json(data)


    def create_update_visible_button(self):
        button = tk.Button(self.root, text="Update visibility based on included nodes", command=self.on_update_visible_button_click)
        button.grid(row=4, column=1, padx=5, pady=5, sticky="ew")

    def on_update_visible_button_click(self):
        try:
            for node in self.met_network.all_network_nodes:
                if node.instance in [met for met in self.met_network.part_of_physics_metabolites]:
                    node.visible_by_mouse = True
                elif node.instance in [met for met in self.met_network.not_part_of_physics_metabolites]:
                    node.visible_by_mouse = False
        except:
            pass
        self.create_listboxes_visibility()
    def create_show_unshow_compartments_button(self):
        button = tk.Button(self.root, text="Show/unshow compartments toggle", command=self.on_show_unshow_compartments_button_click)
        button.grid(row=5, column=0, padx=5, pady=5, sticky="ew")

    def on_show_unshow_compartments_button_click(self):
        self.met_network.show_compartments = not self.met_network.show_compartments

    def create_browse_and_upload_excel_file(self):
        self.entry_var_excel_file = tk.StringVar(self.root)
        self.var_excel_file_path = tk.StringVar(self.root)

        self.entry2 = tk.Entry(self.root, textvariable=self.entry_var_excel_file, state="readonly", width=40)
        self.entry2.grid(row=10, column=2, padx=10, pady=0, sticky="ew")

        button = tk.Button(self.root, text="Browse excel File", command=self.on_browse_excel_file)
        button.grid(row=11, column=2, padx=10, pady=15, sticky="ew")

        button2 = tk.Button(self.root, text="Upload excel File", command=self.on_upload_excel_file)
        button2.grid(row=12, column=2, padx=10, pady=0, sticky="ew")

    def on_browse_excel_file(self):
        file_path = filedialog.askopenfilename(filetypes=[("Excel files", "*.xls;*.xlsx")])
        if file_path:
            self.entry_var_excel_file.set(os.path.basename(file_path))
            self.var_excel_file_path.set(file_path)

    def on_upload_excel_file(self):
        file_path = self.var_excel_file_path.get()
        try:
            with open(file_path, 'r') as file:
                # Read excel data
                self.met_network.df_fluxes_reactions = pd.read_excel(file_path, header = None)
                self.check_if_correct_fluxes_reactions_file()
        except FileNotFoundError:
            print("File not found.")

    def check_if_correct_fluxes_reactions_file(self):
        columns_to_convert = [1, 2]
        for col_idx in columns_to_convert:
            self.met_network.df_fluxes_reactions.iloc[:, col_idx] = pd.to_numeric(self.met_network.df_fluxes_reactions.iloc[:, col_idx], errors='coerce')
        selected_columns = self.met_network.df_fluxes_reactions.iloc[:, [1, 2 ]]  # Adjust column indices as per your DataFrame
        are_numeric = selected_columns.iloc[:, [0, 1]].apply(pd.to_numeric, errors='coerce').notnull().all()
        has_temp_exchange = self.met_network.df_fluxes_reactions.iloc[:, 3].str.startswith(('temporary_exchange', 'MAR')).any()
        if are_numeric.all() and has_temp_exchange:
            self.create_fluxes_and_expression_dicts()
        else:
            print("bad excel file uploaded")
            self.met_network.df_fluxes_reactions = pd.DataFrame()

    def create_fluxes_and_expression_dicts(self):
        self.met_network.df_fluxes_reactions = self.met_network.df_fluxes_reactions.iloc[:, [1, 2, 3 ]]
        self.met_network.df_fluxes_reactions.columns = ["fluxes", "expression", "reaction"]
        for index, row in self.met_network.df_fluxes_reactions.iterrows():
            value_2 = row["expression"]
            value_1 = row["fluxes"]
            value_3 = row["reaction"]
            if '[' in value_3:
                value_3 = value_3.split('[')[0]

            self.met_network.fluxes_dict[str(value_3)] = round(value_1,3)
            self.met_network.expression_dict[str(value_3)] = round(value_2,3)

    def create_canvas_settings(self):
        self.canvas_width_entry = tk.Entry(self.root, width=10)
        self.canvas_height_entry = tk.Entry(self.root, width=10)
        label_width = tk.Label(self.root, text="Canvas Width")
        label_width.grid(row=13, column=0)

        default_width = tk.StringVar(value="900")  # Default width
        self.canvas_width_entry = tk.Entry(self.root, textvariable=default_width)
        self.canvas_width_entry.grid(row=14, column=0)

        label_height = tk.Label(self.root, text="Canvas Height")
        label_height.grid(row=13, column=1)

        default_height = tk.StringVar(value="800")  # Default height
        self.canvas_height_entry = tk.Entry(self.root, textvariable=default_height)
        self.canvas_height_entry.grid(row=14, column=1)

        button = tk.Button(self.root, text="Set Dimensions", command=self.on_canvas_settings_button)
        button.grid(row=14, column=2)

    def on_canvas_settings_button(self):
        self.met_network.canvas_width = int(self.canvas_width_entry.get())
        self.met_network.canvas_height = int(self.canvas_height_entry.get())
    def create_listboxes_visibility(self):
        label1 = tk.Label(self.root, text="visible")
        label1.grid(row=0, column=2, padx=10, pady=5)
        self.left_list_visibility = tk.Listbox(self.root, width=20, height=20)
        self.left_list_visibility.grid(row=1, column=2, padx=10, pady=5)
        label2 = tk.Label(self.root, text="not visible")
        label2.grid(row=0, column=3, padx=10, pady=5)
        self.right_list_visibility= tk.Listbox(self.root, width=20, height=20)
        self.right_list_visibility.grid(row=1, column=3, padx=10, pady=5)

        self.names_left_visibility = [node.instance.name if isinstance(node.instance, MetaboliteNodes) else node.instance.id for node in
                                      self.met_network.all_network_nodes if not isinstance(node.instance, BoundaryNode) and node.visible_by_mouse and node.instance.id !="000fake"]
        self.names_right_visibility =[node.instance.name if isinstance(node.instance, MetaboliteNodes) else node.instance.id for node in
                                      self.met_network.all_network_nodes if not isinstance(node.instance, BoundaryNode) and not node.visible_by_mouse and node.instance.id != "000fake"]

        for item in self.names_left_visibility:
            self.left_list_visibility.insert(tk.END, item)
        for item in self.names_right_visibility:
            self.right_list_visibility.insert(tk.END, item)

        self.left_list_visibility.bind('<Double-Button-1>', lambda event: self.on_double_click_visibility(event, self.left_list_visibility, self.right_list_visibility))
        self.right_list_visibility.bind('<Double-Button-1>', lambda event: self.on_double_click_visibility(event, self.right_list_visibility, self.left_list_visibility))

    def move_item_visibility(self, event, from_list, to_list):
        try:
            selected_item = from_list.get(from_list.curselection())
            for node in self.met_network.all_network_nodes:
                if not isinstance(node.instance, BoundaryNode) and (node.instance.id == selected_item or node.instance.name == selected_item):
                    node.visible_by_mouse = not node.visible_by_mouse
            to_list.insert(tk.END, selected_item)
            from_list.delete(from_list.curselection())
        except:
            pass

    def on_double_click_visibility(self, event, from_list, to_list):
        self.move_item_visibility(event, from_list, to_list)
    def draw_network(self, canvas):
        # Draw nodes
        node_radius = 5
        for node in self.met_network.all_network_nodes:
            if node.node_type != "boundary" and node.node_type != "compartment":
                x = node.x
                y = node.y
                canvas.create_oval(x - node_radius, y - node_radius, x + node_radius, y + node_radius, fill='blue')
                canvas.create_text(x, y - 15, text=node)
        # Draw edges
        for idx, edges in enumerate(self.met_network.edge_logic_dict["reaction_metabolite_attraction"]):
            for edge_idx in range(len(edges)):
                x1 = self.met_network.all_network_nodes[idx].x
                y1 = self.met_network.all_network_nodes[idx].y
                if edges[edge_idx] == 1:
                    x2 = self.met_network.all_network_nodes[edge_idx].x
                    y2 = self.met_network.all_network_nodes[edge_idx].y
                    canvas.create_line(x1, y1, x2, y2, fill='black')


# Create the main Tkinter window
root = tk.Tk()
width = 1500
height = 1200
running = False
# Create an instance of the JSONViewer class
TK_application = Application(root, width, height)

# Start the Tkinter event loop
root.mainloop()

### CURRENT ISSUE:
# Right now we get to Nan values in the force calculation, which probably means the repulsion/attraction isn't working correctly
# also did not add in a thing for not moving fixed notes
# create test code to check repulsion attraction over iterations, instead of doing it like right now

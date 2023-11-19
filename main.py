# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
import os, json, pygame, math, random, markdown
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
        self.canvas = pygame.Surface((self.WIDTH, self.HEIGHT))
        self.screen = pygame.display.set_mode((self.FINAL_WIDTH, self.FINAL_HEIGHT))


class TkinterSetup:
    pass


filepath = "E:\Git\MaCSBio-GEM\General\Functions\Metabolic_tasks\Task322Sample1Alpha0.50Beta0.00_Date2023-10-25_21-51-37\ModelsJSON/newModelMILP_OR_AdjustTask322Sample1Scored2.152average4.105alpha0.500beta0.000.json"

with open(filepath, 'r') as file:
    # Read JSON data
    json_data = json.load(file)


# print(json_data)
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

    def find_input_and_output_instances(self):
        if isinstance(self.instance, ReactionEdges):
            for metabolite, stoichiometry in self.instance.metabolites.items():
                if metabolite in [metabolit.id for metabolit in self.met_network.included_metabolites] and stoichiometry <0:
                    self.input_instances.append(self.met_network.find_metabolite(metabolite))
                elif metabolite in [metabolit.id for metabolit in self.met_network.included_metabolites] and stoichiometry >0:
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
    def __init__(self, nodes, edges):
        self.metabolites_nodes = nodes
        self.reaction_edges = edges
        self.set_amount_of_reactions_per_node()
        self.set_slider()
        self.filter_input_output_reactions()
        self.set_reactions_for_nodes()
        self.commonly_excluded_names = []#["H+", "H2O", "NADH", "NAD+", "NADP+", "NADPH", "FAD", "FADH", "CO2"]
        self.included_metabolites = []
        self.not_included_metabolites = []
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
        self.find_reactions_that_should_generally_be_excluded()
        for metabolite in self.metabolites_nodes:
            if metabolite.amount_of_reactions > self.slider_value and metabolite not in self.not_included_metabolites:
                self.not_included_metabolites.append(metabolite)

        for metabolite in self.metabolites_nodes:
            if (metabolite.fixed_node_input or metabolite.fixed_node_output) and metabolite not in self.included_metabolites:
                self.included_metabolites.append(metabolite)
                try:
                    self.not_included_metabolites.remove(metabolite)
                except:
                    pass
        self.included_metabolites = [item for item in self.metabolites_nodes if item not in self.not_included_metabolites]
        #self.not_included_metabolites = [item for item in self.metabolites_nodes if item not in self.included_metabolites]

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
            #if not (reaction.lb > 0 and reaction.ub < 1000) or not reaction.name.startswith('temporary_exchange_'):
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
        # Create compartment_pseudo_reaction_nodes:
        self.create_compartment_pseudo_reaction_nodes()
        self.create_network_nodes_for_simulation()

        # create edge matrices
        self.create_edges_matrices_dicts()

        # set initial x, y positions for non-boundary nodes
        self.set_x_y_positions_non_boundary_nodes()

        # run network
        self.run_network_simulation()

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
        self.all_network_nodes = []
        # Create all types of nodes
        id_number = 0
        for metabolite in self.included_metabolites:
            self.all_network_nodes.append(NetworkNodesForPhysicsSimulation(0, 0, metabolite, id_number, self))
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

        temp_connection_matrix = np.zeros((node_length,node_length))
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
        self.edge_logic_dict["compartment_compartment_repulsion"]= temp_connection_matrix

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
                        if connected_instance == other_node.instance:
                            edge.append(other_node.id_number)
            edges.append(edge)
        for edge in edges:
            node_id = edge[0]
            connected_nodes = edge[1:]
            temp_connection_matrix[node_id, connected_nodes] = 1
        self.edge_logic_dict["metabolite_compartment_attraction"] = temp_connection_matrix



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
                node.x =  200 + (fixed_output_counter * 25)
                fixed_output_counter += 1
                node.y = height - 50
            elif node.node_type == "metabolite" and node.instance.fixed_node_input:
                node.x = 200+ (fixed_input_counter *25)
                fixed_input_counter += 1
                node.y = 0 + 50
            else:
                node.x = random.randint(200, width-200)
                node.y = random.randint(200, height-200)


    def run_network_simulation(self):
        self.network_physics_simulation_loop()
        # translate problem to instances of classes of nodes+connected edges to other nodes
        # way to calculate distance between all nodes
        # calculate forces per node
        ## repulsion 1/(distance cubed + 0.0001)
        ## special repulsion for compartments and cyclic sets of reactions
        ## attraction per edge towards edge partner, attraction, distance squared
        # move based on force
        # add in cyclic repulsion for creation of circles



    def network_physics_simulation_loop(self):
        # loop
        iterations = 100000
        counter = 0
        try:
            pygame.quit()
        except:
            pass

        running = True
        pygame_on = True
        if pygame_on:
            os.environ['SDL_VIDEO_WINDOW_POS'] = '700,50'
            pygame.init()
            display = pygame.display.set_mode([1800, 1200])
            font = pygame.font.SysFont(None, 18)
        while running:
            if pygame_on:
                display.fill((0, 0, 0))
                for event in pygame.event.get():
                    if event.type == pygame.QUIT:
                        running = False
                    elif event.type == pygame.KEYDOWN:
                        if event.key == pygame.K_F1:
                            running = False
                for idx, node in enumerate(self.all_network_nodes):
                    if node.node_type == "reaction":
                        rect = pygame.Rect(node.x, node.y, 15, 15)
                        pygame.draw.rect(display, (0, 200, 200), rect)
                        text_surface = font.render(node.instance.name, True, (255,255,255))
                        text_rect = text_surface.get_rect(center =rect.center )
                        display.blit(text_surface, text_rect)

                    elif node.node_type == "metabolite" :
                        for idx, compartment in enumerate(self.compartments):
                            if node.instance.compartment == compartment.id:
                                pygame.draw.rect(display, (50+ (idx*50) % 255, 50, 50), pygame.Rect(node.x, node.y, 15, 15))

                for idx, node in enumerate(self.all_network_nodes):

                    x_start = self.all_network_nodes[idx].x
                    y_start = self.all_network_nodes[idx].y

                    for idx2 in range(1, len(self.edge_logic_dict["reaction_metabolite_attraction"][idx])):
                        x_end = 0
                        y_end = 0
                        if self.edge_logic_dict["reaction_metabolite_attraction"][idx][idx2] == 1:
                            x_end = self.all_network_nodes[idx2].x
                            y_end = self.all_network_nodes[idx2].y

                            pygame.draw.line(display, (255, 255, 255), (x_start, y_start), (x_end, y_end))
                for idx, node in enumerate(self.all_network_nodes):

                    for idx2 in range(1, len(self.edge_logic_dict["reaction_metabolite_attraction"][idx])):
                        x_end = 0
                        y_end = 0
                        if self.edge_logic_dict["reaction_metabolite_attraction"][idx2][idx] == 1:
                            x_start = self.all_network_nodes[idx2].x
                            y_start = self.all_network_nodes[idx2].y
                            x_end = self.all_network_nodes[idx].x
                            y_end = self.all_network_nodes[idx].y

                            pygame.draw.line(display, (255, 255, 255), (x_start, y_start), (x_end, y_end))
                pygame.display.flip()

            self.calculate_forces_through_network()
            self.move_nodes_based_on_forces()
            counter += 1
            if counter >= iterations:
                running = False
        if pygame_on:
            pygame.quit()

    def calculate_distance_matrix(self):
        coordinates = np.array([(obj.x, obj.y) for obj in self.all_network_nodes])

        diff = coordinates[:, np.newaxis, :] - coordinates[np.newaxis, :, :]
        distance_matrix = np.linalg.norm(diff, axis=2)
        angle_matrix = np.arctan2(diff[:, :, 1], diff[:, :, 0])
        # print(distance_matrix)
        # print(self.edge_logic_dict["normal_repulsion"])
        return distance_matrix, angle_matrix

    def move_nodes_based_on_forces(self):
        SMALL_CONSTANT = 0.000015
        for node in self.all_network_nodes:
            if node.fixed:
                node.x += node.x_vector * SMALL_CONSTANT

            else:
                node.x += node.x_vector * SMALL_CONSTANT
                node.x = node.x % width
                node.y += node.y_vector * SMALL_CONSTANT
                node.y = node.y %height

            if node.x > width:
                node.x = width
            elif node.x < 0:
                node.x =0
            if node.y > width:
                node.y = width
            elif node.y <0:
                node.y = 0
    def calculate_forces_through_network(self):
        distance_matrix, angle_matrix = self.calculate_distance_matrix()
        distance_matrix = distance_matrix+1
        cos_angles = np.cos(angle_matrix)
        sin_angles = np.sin(angle_matrix)

        # # General repulsion
        np.fill_diagonal(distance_matrix, 1000000)
        force_magnitudes =  (200/ np.power(np.sqrt(distance_matrix), 2))
        force_magnitudes = np.where(force_magnitudes < 0.0000000001, 0, force_magnitudes)
        x_forces = np.sum(force_magnitudes * cos_angles * self.normal_repulsion_value * self.edge_logic_dict["normal_repulsion"], axis=1)
        y_forces = np.sum(force_magnitudes * sin_angles * self.normal_repulsion_value * self.edge_logic_dict["normal_repulsion"], axis=1)


        # Compartment repulsion
        force_magnitudes = 30000 / np.power(np.sqrt(distance_matrix), 2)
        force_magnitudes = np.where(force_magnitudes < 0.0000001, 0, force_magnitudes)
        x_forces = x_forces + np.sum(force_magnitudes * cos_angles * self.comp_comp_repulsion_value * self.edge_logic_dict["compartment_compartment_repulsion"], axis = 1)
        y_forces = y_forces + np.sum(force_magnitudes * sin_angles * self.comp_comp_repulsion_value * self.edge_logic_dict["compartment_compartment_repulsion"], axis=1)

        # Fixed metabolite repulsion
        force_magnitudes = 9000 / np.power(np.sqrt(distance_matrix), 2)
        force_magnitudes = np.where(force_magnitudes < 0.0000001, 0, force_magnitudes)
        x_forces = x_forces + np.sum(force_magnitudes * cos_angles * self.fixed_metabolite_repulsion_value * self.edge_logic_dict["fixed_metabolite_repulsion"], axis = 1)


        # Boundary Repulsion
        force_magnitudes = 500000 / np.power(np.sqrt(distance_matrix), 4)
        force_magnitudes = np.where(force_magnitudes < 0.0000001, 0, force_magnitudes)
        x_forces = x_forces + np.sum(force_magnitudes * cos_angles * self.boundary_repulsion_value * self.edge_logic_dict["boundary_repulsion"], axis = 1)
        y_forces = y_forces + np.sum(force_magnitudes * sin_angles * self.boundary_repulsion_value * self.edge_logic_dict["boundary_repulsion"], axis=1)

        # rxn_met attraction
        np.fill_diagonal(distance_matrix, 0)
        force_magnitudes = -np.power(distance_matrix, 1)/10
        force_magnitudes = np.where(force_magnitudes > -0.0001, 0, force_magnitudes)
        x_forces = x_forces + np.sum(force_magnitudes * cos_angles * self.rxn_met_attraction_value * self.edge_logic_dict["reaction_metabolite_attraction"], axis = 1)
        y_forces = y_forces + np.sum(force_magnitudes * sin_angles * self.rxn_met_attraction_value * self.edge_logic_dict["reaction_metabolite_attraction"], axis=1)

        # met_comp attraction
        force_magnitudes = -np.power(distance_matrix, 2)/100
        force_magnitudes = np.where(force_magnitudes < 0.0001, 0, force_magnitudes)
        x_forces = x_forces + np.sum(force_magnitudes * cos_angles * self.metabolite_compartment_attraction_value * self.edge_logic_dict["metabolite_compartment_attraction"], axis = 1)
        y_forces = y_forces + np.sum(force_magnitudes * sin_angles * self.metabolite_compartment_attraction_value  * self.edge_logic_dict["metabolite_compartment_attraction"], axis=1)

        for node, x_force, y_force in zip(self.all_network_nodes, x_forces, y_forces):
            node.x_vector = x_force
            node.y_vector = y_force



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

    def browse_file(self):
        file_path = filedialog.askopenfilename(filetypes=[("JSON files", "*.json")])
        self.entry_var.set(file_path)

    def key_interactions(self, event):
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
                self.create_slider_normal_repulsion()
                self.create_slider_comp_comp_repulsion()
                self.create_slider_boundary_repulsion()
                self.create_slider_fixed_metabolite_slider()
                self.create_slider_rxn_met_attraction_slider()
                self.create_slider_met_comp_attraction_slider()
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
        self.load_network_button = tk.Button(self.root, text="(re)load network", command=self.load_network)
        self.load_network_button.grid(row=1, column=2, padx=10, pady=5)

    def load_network(self):
        self.met_network.calculate_network()
        #self.create_and_draw_network_canvas()

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

        self.canvas = tk.Canvas(self.network_window, width= self.WIDTH, height=self.HEIGHT)
        self.draw_network(self.canvas)
        self.canvas.place(in_=self.network_window, x=0, y=0)

        # nodes =
        # edges

    def create_listboxes(self):

        label1 = tk.Label(self.root, text="included")
        label1.grid(row=0, column=0, padx=10, pady=5)
        self.left_list = tk.Listbox(self.root, width=20, height=20)
        self.left_list.grid(row=1, column=0, padx=10, pady=5)
        label2 = tk.Label(self.root, text="not included")
        label2.grid(row=0, column=1, padx=10, pady=5)
        self.right_list = tk.Listbox(self.root, width=20, height=20)
        self.right_list.grid(row=1, column=1, padx=10, pady=5)

        self.names_left = [name.name for name in self.met_network.included_metabolites]
        self.names_right = [name.name for name in self.met_network.not_included_metabolites]

        for item in self.names_left:
            self.left_list.insert(tk.END, item)
        for item in self.names_right:
            self.right_list.insert(tk.END, item)

        self.left_list.bind('<Double-Button-1>', lambda event: self.on_double_click(event, self.left_list, self.right_list))
        self.right_list.bind('<Double-Button-1>', lambda event: self.on_double_click(event, self.right_list, self.left_list))


    def move_item(self, event, from_list, to_list):
        selected_item = from_list.get(from_list.curselection())
        to_list.insert(tk.END, selected_item)
        from_list.delete(from_list.curselection())

    def on_double_click(self, event, from_list, to_list):
        self.move_item(event, from_list, to_list)
        self.met_network.included_metabolites = [self.met_network.find_metabolite(met) for met in self.left_list.get("0", tk.END)]
        self.met_network.not_included_metabolites = [self.met_network.find_metabolite(met) for met in self.right_list.get("0", tk.END)]

    def update_slider(self, value):
        value = int(value)
        self.met_network.slider_value = value
        self.right_list.destroy()
        self.left_list.destroy()
        self.met_network.update_included_not_included_based_on_slider()
        self.create_listboxes()

    def create_slider(self):
        label = tk.Label(self.root, text="Metabolites with more connections are excluded")
        label.grid(row=2, column=0, columnspan=2, padx=10, pady=5, sticky="ew")

        scale = tk.Scale(self.root, from_=self.met_network.slider_start, to=self.met_network.slider_end, orient=tk.HORIZONTAL, command=self.update_slider)
        scale.grid(row=3, column=0, columnspan=2, padx=10, pady=5, sticky="ew")
        scale.set(self.met_network.slider_value)

    def create_slider_normal_repulsion(self):
        #self.root.grid_rowconfigure(1, weight=1, maxsize=50)  # Adjust minsize as needed
        label = tk.Label(self.root, text="Normal Repulsion")
        label.grid(row=2, column=3, columnspan=2, padx=10, pady=5, sticky="ew")

        scale = tk.Scale(self.root, from_=1, to=5000, orient=tk.HORIZONTAL, command=self.update_normal_repulsion_slider)
        scale.grid(row=3, column=3, columnspan=2, padx=10, pady=5, sticky="ew")
        scale.set(self.met_network.normal_repulsion_value)


    def update_normal_repulsion_slider(self, value):
        self.met_network.normal_repulsion_value = int(value)
        # self.met_network.stop_simulation()
        # self.met_network.start_simulation()

    def create_slider_comp_comp_repulsion(self):
        # self.root.grid_rowconfigure(1, weight=1, maxsize=50)  # Adjust minsize as needed
        label = tk.Label(self.root, text="Intercompartment Repulsion")
        label.grid(row=4, column=3, columnspan=2, padx=10, pady=5, sticky="ew")

        scale = tk.Scale(self.root, from_=1, to=5000, orient=tk.HORIZONTAL, command=self.update_comp_comp_repulsion_slider)
        scale.grid(row=5, column=3, columnspan=2, padx=10, pady=5, sticky="ew")
        scale.set(self.met_network.comp_comp_repulsion_value)
    def update_comp_comp_repulsion_slider(self, value):
        self.met_network.comp_comp_repulsion_value = int(value)
        # self.met_network.stop_simulation()
        # self.met_network.start_simulation()

    def create_slider_boundary_repulsion(self):
        # self.root.grid_rowconfigure(1, weight=1, maxsize=50)  # Adjust minsize as needed
        label = tk.Label(self.root, text="Boundary Repulsion")
        label.grid(row=6, column=3, columnspan=2, padx=10, pady=5, sticky="ew")

        scale = tk.Scale(self.root, from_=1, to=5000, orient=tk.HORIZONTAL, command=self.update_boundary_repulsion_slider)
        scale.grid(row=7, column=3, columnspan=2, padx=10, pady=5, sticky="ew")
        scale.set(self.met_network.boundary_repulsion_value)
    def update_boundary_repulsion_slider(self, value):
        self.met_network.boundary_repulsion_value = int(value)
        # self.met_network.stop_simulation()
        # self.met_network.start_simulation()

    def create_slider_rxn_met_attraction_slider(self):
        # self.root.grid_rowconfigure(1, weight=1, maxsize=50)  # Adjust minsize as needed
        label = tk.Label(self.root, text="Reaction-Metabolite Attraction")
        label.grid(row=10, column=3, columnspan=2, padx=10, pady=5, sticky="ew")

        scale = tk.Scale(self.root, from_=1, to=5000, orient=tk.HORIZONTAL, command=self.update_rxn_met_attraction_slider)
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

        scale = tk.Scale(self.root, from_=1, to=5000, orient=tk.HORIZONTAL, command=self.update_met_comp_slider)
        scale.grid(row= 13, column=3, columnspan=2, padx=10, pady=5, sticky="ew")
        scale.set(self.met_network.metabolite_compartment_attraction_value)
    def update_met_comp_slider(self, value):
        self.met_network.metabolite_compartment_attraction_value = int(value)
        # self.met_network.stop_simulation()
        # self.met_network.start_simulation()

    def create_slider_fixed_metabolite_slider(self):
        # self.root.grid_rowconfigure(1, weight=1, maxsize=50)  # Adjust minsize as needed
        label = tk.Label(self.root, text="Input_metabolite-Input_metabolite Repulsion")
        label.grid(row=8, column=3, columnspan=2, padx=10, pady=5, sticky="ew")

        scale = tk.Scale(self.root, from_=1, to=5000, orient=tk.HORIZONTAL, command=self.update_fixed_metabolite_repulsion_slider)
        scale.grid(row=9, column=3, columnspan=2, padx=10, pady=5, sticky="ew")
        scale.set(self.met_network.fixed_metabolite_repulsion_value)

    def update_fixed_metabolite_repulsion_slider(self, value):
        self.met_network.fixed_metabolite_repulsion_value = int(value)
        # self.met_network.stop_simulation()
        # self.met_network.start_simulation()


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
                if edges[edge_idx] ==1:
                    x2 = self.met_network.all_network_nodes[edge_idx].x
                    y2 = self.met_network.all_network_nodes[edge_idx].y
                    canvas.create_line(x1, y1, x2, y2, fill='black')


# Create the main Tkinter window
root = tk.Tk()
width = 1200
height = 900
# Create an instance of the JSONViewer class
TK_application = Application(root, width, height)

# Start the Tkinter event loop
root.mainloop()

### CURRENT ISSUE:
# Right now we get to Nan values in the force calculation, which probably means the repulsion/attraction isn't working correctly
# also did not add in a thing for not moving fixed notes
# create test code to check repulsion attraction over iterations, instead of doing it like right now

#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
import networkx as nx
import networkx.drawing
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics
import hashlib
import pickle

__author__ = "KERMARREC Maxime"
__copyright__ = "Universite de Paris"
__credits__ = ["KERMARREC Maxime"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "KERMARREC Maxime"
__email__ = "maximekermarrec14@gmail.com"
__status__ = "Developpement"

def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=21, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    return parser.parse_args()


def read_fastq(fichier):
    with open(fichier, "r") as filin:
        for ligne in filin:
            yield next(filin).strip()
            next(filin)
            next(filin)


def cut_kmer(seq, k):
    j=0
    for i in range(k, len(seq)+1):
        yield seq[j:i]
        j+=1


def build_kmer_dict(fichier, k):
    dico = {}         
    for seq in read_fastq(fichier):
        seq = str(seq)
        for kmer in cut_kmer(seq, k):
            if kmer in dico:
                dico[kmer] = dico[kmer] + 1
            else:
                dico[kmer] = 1
    return dico


def build_graph(dico):
    graphe = nx.DiGraph()
    for key in dico:
        prefixe = str(key)[:-1]
        suffixe = str(key)[1:]
        graphe.add_edge(prefixe, suffixe, weight=int(dico[key]))
    
    return(graphe)

def get_starting_nodes (graphe):
    entree = []    
    for node in graphe.nodes:
        if list(graphe.predecessors(node)) == []:
            entree.append(node)
    return entree

def get_sink_nodes (graphe):
    sortie = []    
    for node in graphe.nodes:
        if list(graphe.successors(node)) == []:
            sortie.append(node)
    return sortie

def get_contigs (graphe,entree,sortie):
    list_contig = []
    for i in entree:
        for j in sortie:
            for path in nx.all_simple_paths(graphe, i, j):
                contig = path[0][:1]
                for k in range(1, len(path)):
                    if k == len(path)-1:
                        contig = contig + path[k]
                    else :
                        contig = contig + path[k][:1]
                a = (str(contig), int(len(contig)))
                list_contig.append(a)
    return list_contig

def save_contigs (list_contig, text):
    with open(text, "w") as filout:
        i = 0
        for contig, longueur in list_contig:
            filout.write(">contig_{} len={}\n".format(i,longueur))
            i = i + 1
            ligne = fill(contig)
            filout.write("{}\n".format(ligne))

def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def std (liste):
    return statistics.stdev(liste)

def path_average_weight (graphe, chemin):
    poids = []
    for i in range(0, len(chemin)-1):
        poids.append(graphe[chemin[i]][chemin[i+1]]["weight"])
    return statistics.mean(poids)

def path_weight(graphe, chemin):
    poids = 0
    for i in range(0, len(chemin)-1):
        poids += graphe[chemin[i]][chemin[i+1]]["weight"]
    return poids


def remove_paths (graphe, chemin, delete_entry_node, delete_sink_node):
    for i in chemin:
        for j in range(1, len(i)-1):
            graphe.remove_node(i[j])
        if delete_entry_node == True:
            graphe.remove_node(i[0])
        if delete_sink_node == True:
            graphe.remove_node(i[len(i)-1])
    return graphe

def select_best_path (graphe, chemin, long_chemin, poids_chemin, delete_entry_node = False, delete_sink_node = False):
    best_poids = max(poids_chemin)
    best_long = max(long_chemin)
    index_poids = []
    index_taille = []
    for i in range(0, len(poids_chemin)):
        if poids_chemin[i] == best_poids:
            index_poids.append(i)
    for i in index_poids:
        if long_chemin[i] == best_long:
            index_taille.append(i)
    if len(index_poids) == 1: 
        bon_chemin = chemin[index_poids[0]]
        del chemin[index_poids[0]]
        chemin = conservation_nodes(bon_chemin, chemin)
        graphe = remove_paths (graphe, chemin, delete_entry_node, delete_sink_node)
        return graphe
    if len(index_taille) == 1:
        bon_chemin = index_taille[0]
        del chemin[index_taille[0]]
        chemin = conservation_nodes(bon_chemin, chemin)
        graphe = remove_paths (graphe, chemin, delete_entry_node, delete_sink_node)
        return graphe
    else : 
        hasard = random.randint(0,len(index_taille)-1)
        bon_chemin = chemin[hasard]
        del chemin[hasard]
        chemin = conservation_nodes(bon_chemin, chemin)
        graphe = remove_paths (graphe, chemin, delete_entry_node, delete_sink_node)
        return graphe

def conservation_nodes (conserve, suppr):
    definitive_list = []
    a_suppr = []
    for path in suppr:
        a_suppr.append(list(set(path) - set(conserve)))
    nouveau_suppr = []
    nouveau_suppr.append(conserve[0])
    for path_1 in a_suppr:
        for path_2 in path_1:
            nouveau_suppr.append(path_2)
    nouveau_suppr.append(conserve[-1])
    
    nouveau_suppr = unique(nouveau_suppr)
    
    definitive_list.append(nouveau_suppr)
    return definitive_list

def unique(list1): 
  
    # intilize a null list 
    unique_list = [] 
      
    for x in list1: 
        # check if exists in unique_list or not 
        if x not in unique_list: 
            unique_list.append(x) 
    return unique_list

def solve_bubble (graphe, ancestor, descendant):
    pathways = nx.all_simple_paths(graphe, ancestor, descendant, cutoff=None)
    poids_chemin = []
    long_chemin = []
    chemin = []
    for path in pathways:
        chemin.append(path)
        poids_chemin.append(path_average_weight(graphe, path))
        long_chemin.append(len(path))
    select_best_path (graphe, chemin, long_chemin, poids_chemin, delete_entry_node = False, delete_sink_node = False)
    return graphe


def simplify_bubbles (graphe):
    ensemble_chemin = []
    ensemble_ancestre = []
    ensemble_descendant = []
    start = get_starting_nodes (graphe)
    end = get_sink_nodes(graphe)
    for node_A in start:
        for node_B in end:
            if node_A != node_B:
                chemin = nx.all_simple_paths(graphe, node_A, node_B)
                if len(list(chemin)) > 1:
                    chemin = nx.all_simple_paths(graphe, node_A, node_B)
                    ensemble_ancestre.append(node_A)
                    ensemble_descendant.append(node_B)
                    for path in chemin:
                        ensemble_chemin.append(list(chemin))
    i=0
    save_path = 0
    save_long_path = 0
    for path_1 in ensemble_chemin:
        for path_2 in path_1:
            if save_long_path < len(path_2):
                save_path = path_1
                save_long_path = len(path_2)
                ancestre = path_2[0]
                descendant = path_2[-1]
    graphe = solve_bubble(graphe, ancestre, descendant)
    return graphe

def deter_min(ensemble_chemin):
    min = [ensemble_chemin[0]]
    taille = len(ensemble_chemin[0])
    for chemin in ensemble_chemin:
        if len(chemin) < taille:
            taille = len(chemin)
            min = chemin
    return min

def solve_entry_tips (graphe, entree):
    ensemble_chemin = []
    ensemble_poids = []
    ensemble_taille = []
    starting_nodes = get_sink_nodes (graphe)
    for start in starting_nodes:
        for node in entree:
            chemin = list(nx.all_simple_paths(graphe, node, start))
            for path in chemin:
                ensemble_chemin.append(path)
    for chemins in ensemble_chemin:
        ensemble_poids.append(path_weight(graphe, chemins))
        ensemble_taille.append(len(chemins))
    graphe = select_best_path (graphe, ensemble_chemin, ensemble_taille, ensemble_poids, delete_entry_node = False, delete_sink_node = False)
    return(graphe)

def solve_out_tips (graphe, sortie):
    ensemble_chemin = []
    ensemble_poids = []
    ensemble_taille = []
    starting_nodes = get_starting_nodes (graphe)
    for start in starting_nodes:
        for node in sortie:
            chemin = list(nx.all_simple_paths(graphe, start, node))
            for path in chemin:
                ensemble_chemin.append(path)
    for chemins in ensemble_chemin:
        ensemble_poids.append(path_weight(graphe, chemins))
        ensemble_taille.append(len(chemins))
    graphe = select_best_path (graphe, ensemble_chemin, ensemble_taille, ensemble_poids, delete_entry_node = False, delete_sink_node = False)
    return(graphe)


#==============================================================
# Main program
#==============================================================

def main():
    """
    Main program function
    """


    # Get arguments and using it

    args = get_arguments()
    dico = build_kmer_dict(args.fastq_file, args.kmer_size)

    graphe = build_graph(dico)


    graphe = simplify_bubbles (graphe)
    
    node_entree = []

    for node_A in graphe.nodes:
        if len(list(graphe.predecessors(node_A))) >= 2:
            graphe.predecessors(node)
            node_entree.append(list(node_A))
    if len(node_entree) != 0:
        graphe = solve_entry_tips(graphe, node_entree)

    node_sortie = []
    for node_B in graphe.nodes:
        if len(list(graphe.predecessors(node_B))) >= 2:
            node_sortie.append(list(node_B))

    if len(node_sortie) != 0:
        graphe = solve_entry_tips(graphe, node_sortie)

    entree = get_starting_nodes(graphe)
    sortie = get_sink_nodes(graphe)
    list_contig = get_contigs(graphe, entree, sortie)
    save_contigs (list_contig, args.output_file)

if __name__ == '__main__':
    main()


#https://docs.google.com/presentation/d/1tfCPs2C9FG4kNrerTdnsLnpK677Tq8f21ipq2jDRbMY/edit#slide=id.g44f5588588_0_92
#https://docs.google.com/document/d/1P4v3bHbSurD7RXA-ldVwmtNKGvWsnBae51RMGye_KLs/edit
#https://python.sdv.univ-paris-diderot.fr/cours-python.pdf
#https://networkx.github.io/documentation/stable/tutorial.html#accessing-edges-and-neighbors
#https://networkx.github.io/documentation/stable/reference/functions.html
#https://iut-info.univ-reims.fr/users/blanchard/ISN20181218/utilisation-de-la-bibliotheque-networkx.html
#https://networkx.github.io/documentation/networkx-1.10/reference/generated/networkx.classes.function.nodes.html#networkx.classes.function.nodes
#https://courspython.com/tuple.html
#https://fr.wikipedia.org/wiki/FASTA_(format_de_fichier)
#https://openclassrooms.com/forum/sujet/python-interpreter-n-retour-chariot-en-parametres-78721


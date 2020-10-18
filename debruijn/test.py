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
import matplotlib
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics
import matplotlib.pyplot as plt
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
    #print(dico)
    return dico


def build_graph(dico):
    G = nx.DiGraph()
    for key in dico:
        prefixe = str(key)[:-1]
        suffixe = str(key)[1:]
        G.add_edge(prefixe, suffixe, weight=int(dico[key]))
    for n, nbrs in G.adj.items():
        for nbr, eattr in nbrs.items():
            wt = eattr['weight']
            if wt > 0.5: pass#print(f"({n}, {nbr}, {wt})")
    
    return(G)

def get_starting_nodes (G):
    entree = []    
    for node in G.nodes:
        if list(G.predecessors(node)) == []:
            entree.append(node)
    return entree

def get_sink_nodes (G):
    sortie = []    
    for node in G.nodes:
        if list(G.successors(node)) == []:
            sortie.append(node)
    return sortie

def get_contigs (G,entree,sortie):
    list_contig = []
    for i in entree:
        for j in sortie:
            for path in nx.all_simple_paths(G, i, j):
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

def path_average_weight (G, chemin):
    poids = []
    for i in range(0, len(chemin)-1):
        poids.append(G[chemin[i]][chemin[i+1]]["weight"])
    return statistics.mean(poids)



def remove_paths (G, chemin, delete_entry_node, delete_sink_node):
    for i in chemin:
        for j in range(1, len(i)-1):
            G.remove_node(i[j])
        if delete_entry_node == True:
            G.remove_node(i[0])
        if delete_sink_node == True:
            G.remove_node(i[len(i)-1])
    return G

def select_best_path (G, chemin, long_chemin, poids_chemin, delete_entry_node = False, delete_sink_node = False):
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
        G = remove_paths (G, chemin, delete_entry_node, delete_sink_node)
        return G
    if len(index_taille) == 1:
        bon_chemin = index_taille[0]
        del chemin[index_taille[0]]
        chemin = conservation_nodes(bon_chemin, chemin)
        G = remove_paths (G, chemin, delete_entry_node, delete_sink_node)
        return G
    else : 
        hasard = random.randint(0,len(index_taille)-1)
        bon_chemin = chemin[hasard]
        del chemin[hasard]
        chemin = conservation_nodes(bon_chemin, chemin)
        G = remove_paths (G, chemin, delete_entry_node, delete_sink_node)
        return G

def conservation_nodes (conserve, suppr):
    definitive_list = []
    a_suppr = []
    for path in suppr:
        print(conserve)
        print(path)
        a_suppr.append(list(set(path) - set(conserve)))
    nouveau_suppr = []
    nouveau_suppr.append(conserve[0])
    for path_1 in a_suppr:
        for path_2 in path_1:
            nouveau_suppr.append(path_2)
    nouveau_suppr.append(conserve[-1])
    definitive_list.append(nouveau_suppr)
    print(definitive_list)
    return definitive_list



def solve_bubble (G, ancestor, descendant):
    pathways = nx.all_simple_paths(G, ancestor, descendant, cutoff=None)
    poids_chemin = []
    long_chemin = []
    chemin = []
    for path in pathways:
        chemin.append(path)
        poids_chemin.append(path_average_weight(G, path))
        long_chemin.append(len(path))
    select_best_path (G, chemin, long_chemin, poids_chemin, delete_entry_node = False, delete_sink_node = False)
    return G


def simplify_bubbles (G):
    ensemble_chemin = []
    ensemble_ancestre = []
    ensemble_descendant = []
    for node_A in G.nodes:
        for node_B in G.nodes:
            if node_A != node_B:
                chemin = nx.all_simple_paths(G, node_A, node_B)
                if len(list(chemin)) > 1:
                    chemin = nx.all_simple_paths(G, node_A, node_B)
                    ensemble_ancestre.append(node_A)
                    ensemble_descendant.append(node_B)
                    for path in chemin:
                        ensemble_chemin.append(list(chemin))
    print(ensemble_chemin)
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
    G = solve_bubble(G, ancestre, descendant)
    return G


def solve_entry_tips ():

#def solve_out_tips ():


#==============================================================
# Main program
#==============================================================

def main():
    """
    Main program function
    """

    # Get arguments and using it

    

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


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
            yield next(filin)
            next(filin)
            next(filin)


def cut_kmer(seq, k):
    for i in range(0, len(seq), k):
        yield seq[i:i+k]


def build_kmer_dict(fichier, k):
    dico = {}         
    for seq in read_fastq(fichier):
        seq = str(seq)[:-1]
        for kmer in cut_kmer(seq, k):
            if kmer in dico:
                dico[kmer] += 1
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
            for path in nx.all_simple_paths(G, i, j, cutoff=None):
                #print(path)
                #print("\n")
                contig = path[0]
                #print(len(path))
                for k in range(len(path)):
                    contig = contig + path[k][:1]
            #print(contig)
                a = (str(contig), len(contig))
                list_contig.append(a)
    print(list_contig)

    return list_contig

def save_contigs (list_contig, text):
    with open(text, "w") as filout:
        i = 1
        for contig, longueur in list_contig:
            filout.write("contig_{} len={}\n".format(i,longueur))
            i = i + 1
            ligne = fill(contig)
            filout.write("{}\n".format(ligne))

def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

#def std (liste):
    #return statistics.stdev(liste)
#def path_average_weight ():
#def remove_paths ():
#def select_best_path ():
#def solve_bubble ():
#def simplify_bubbles ():
#def solve_entry_tips ():
#def solve_out_tips ():


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
    G = build_graph(dico)
    entree = get_starting_nodes(G)
    sortie = get_sink_nodes(G)
    list_contig = get_contigs(G, entree, sortie)
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


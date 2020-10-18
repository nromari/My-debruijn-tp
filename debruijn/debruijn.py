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
import matplotlib
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics
import networkx as nx

__author__ = "Noura Ahmar-Erras"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Noura Ahmar-Erras"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Noura Ahmar-Erras"
__email__ = "n_a_e@hotmail.fr"
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

def read_fastq(file):

    """take a fastq file as argument
    and return a sequence generator
    """

    with open(file, "r") as fasta_file:
        for line in fasta_file:
            yield next(fasta_file).strip()
            next(fasta_file)
            next(fasta_file)

def cut_kmer(seq, k):

    """take a sequence and a k-mer lenght as argument
	and return a k-mer gennerator
	"""

    for i in range(len(seq)-k+1):
        yield seq[i:i+k]

def build_kmer_dict(file, k):

    """take a fasta file and a k-mer lenght as argument
    and return a k-mer count dictionary
    """
    dico_kmer = {}
    seq = read_fastq(file)
    kmer_liste = cut_kmer(seq, k)
    for kmer in kmer_liste:
        if kmer in dico_kmer:
            dico_kmer[kmer] += 1
        else:
            dico_kmer[kmer] = 1
    return dico_kmer

def build_graph(dico_kmer):

    """take a K-mer dictionary as argument
    and retrun a digraph
    """

    g = nx.DiGraph()
    for kmer in dico_kmer:
        suffix = kmer[1:]
        g.add_node(suffix)
        prefix = kmer[:-1]
        g.add_node(prefix)
        g.add_edge(prefix, suffix, weight=dico_kmer[kmer])
    return g

    """visualisation de g :
    	nx.draw(g, with_labels=True)
    	plt.show()
    """

def get_starting_nodes(g):

    """take a de Bruijn graph as argument
    and return starting nodes
    """

    entry_nodes = []
    for node in g.nodes():
        if len(list(g.predecessors(node))) == 0:
            entry_nodes.append(node)
    return entry_nodes

def get_sink_nodes(g):

    """take a de Bruijn graph as argument
    and return exit nodes
    """
    exit_nodes = []
    for node in g.nodes():
        if len(list(g.successors(node))) == 0:
            exit_nodes.append(node)
    return exit_nodes

def get_contigs(g, entry_nodes, exit_nodes):

    """take a de Bruijn graph, a starting nodes list and a exit nodes list
    as argument and return a tuple of contig
    """
 
    contigs_tuple = []
    for node_start in entry_nodes:
        for node_end in exit_nodes:
            for path in nx.all_simple_paths(g, node_start, node_end):
                contig = ""
                for i in enumerate(path):
                    if not contig:
                        contig += path[i]
                    else:
                        contig += path[i][-1]
                contigs_tuple.append((contig, len(contig)))
    return contigs_tuple

def fill(text, width=80):

    """Split text with a line return to respect fasta format
    """

    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def save_contigs(contigs_tuple, file_name):

    """Save contigs in a fasta format file
    """

    with open(file_name, "w") as filout:
        for index, contig in enumerate(contigs_tuple):
            filout.write("> contig_" + str(index + 1) + " len = " + str(contig[1]) + "\n")
            filout.write(fill(contig[0]))
            filout.write("\n")

def std(values):

    """standard deviation
	"""

    return statistics.stdev(values)

def path_average_weight(g, path):

    """take a garph and a path as argument
    return the average weight of the path
    """

    weight = 0
    for i in range(len(path)-1):
        weight += g.edges[path[i], path[i+1]]["weight"]    
    av_weight = weight/(len(path)-1)
    return av_weight

def remove_paths(g, path_liste, delete_entry_node, delete_sink_node):

	"""take a graphe and a pathway liste as argument
	return a cleaned graph
	"""

	for i in range(0, len(path_liste)):
	    if delete_entry_node == True:
	        g.remove_node(path_liste[i][0])
	    if delete_sink_node == True:
	        g.remove_node(path_liste[i][-1])
        g.remove_nodes_from(path_liste[i][1:-1])
	return g

def select_best_path(g, paths, len_paths, weight_paths,
                     delete_entry_node=False,
                     delete_sink_node=False):

    """arguments : graph
                   a path list
                   a length paths list
                   an average length list for each path
                   the boleans deletion of entry and sink nodes
        Return a cleaned graph.
    """

    cleaned_g = g

    wm_paths = []
    wm_weight = []
    wm_length = []
    lm_paths = []
    lm_weight = []
    lm_length = []
    undes_paths = []
    for i in range(0, len(paths)):
        if weight_paths[i] == max(weight_paths):
            wm_paths.append(paths[i])
            wm_weight.append(weight_paths[i])
            wm_length.append(len_paths[i])

    for i in range(0, len(wm_paths)):
        if wm_length[i] == max(wm_length):
            lm_length.append(wm_length[i])
            lm_paths.append(wm_paths[i])
            lm_weight.append(wm_weight[i])

    for path in paths:
        if path not in lm_paths:
            undes_paths.append(path)
    cleaned_g = remove_paths(cleaned_g, undes_paths, delete_entry_node, delete_sink_node)

    return cleaned_g
 
def solve_bubble(g, ancestor, descendant):

    """
    Take an ancestor node, a descendant node an a node as argument
    return: a bubble cleaned graph
    """

    cleaned_g = g

    good_paths_list = []
    len_paths_list = []
    weight_paths_list = []

    for path in nx.all_simple_paths(g, ancestor, descendant):
        good_paths_list.append(path)
        len_paths_list.append(len(path))
        weight_paths_list.append(path_average_weight(cleaned_g, path))

    bubble_cleaned_g = select_best_path(g, good_paths_list, len_paths_list, weight_paths_list,
	                                    delete_entry_node=False, delete_sink_node=False)

    return bubble_cleaned_g

def simplify_bubbles():
    """non réalisée
    """
    pass

def solve_entry_tips():
    """non réalisée
    """
    pass

def solve_out_tips():
    """non réalisée
    """
    pass
#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    dico_kmer = build_kmer_dict(args.fastq_file, args.kmer_size)
    g = build_graph(dico_kmer)
    entry_nodes = get_starting_nodes(g)
    exit_nodes = get_sink_nodes(g)
    contigs_tuple = get_contigs(g, entry_nodes, exit_nodes)
    save_contigs(contigs_tuple, args.output_file)
    g = remove_paths(g, path_liste, delete_entry_node, delete_sink_node)
    cleaned_g = select_best_path(g, paths, len_paths,weight_paths, delete_entry_node=False,
                                 delete_sink_node=False)
    bubble_cleaned_g = solve_bubble(g, ancestor, descendant)


if __name__ == '__main__':
    main()

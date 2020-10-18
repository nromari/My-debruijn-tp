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
import matplotlib
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics

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
            yield next(fasta_file)[:-1]
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
    
    contigs = []
    for node_start in entry_nodes:
        for node_end in exit_nodes:
            for path in nx.all_simple_paths(g, node_start, node_end):
                contig = ""
                for i in enumerate(path) :
                    if not contig:
                        contig += path[i]
                    else:
                        contig += path[i][-1]
                contigs.append( (contig, len(contig)) )
    return contigs

#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

if __name__ == '__main__':
    main()

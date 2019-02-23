#!/usr/bin/env python
#Algo PS #3
#Stephanie Wankowicz

#Processing Files

#read in sequences
import os
import pandas as pd
import numpy as np

def read_seq(file_name):
    '''
    INPUT:
    OUTPUT:
    Given a .fa file, this function will read in the sequence & put it into a list
    '''
    seq1=""
    for line in open(file_name):
        if not line.startswith('>'):
            seq1=seq1+line.rstrip()
    #seq1 = seq1
    return seq1

def read_scoring_mat(file_name):
  '''
  INPUT:
  OUTPUT:
  '''
  scoring_mat = pd.read_fwf(file_name,header = 0, comment='#')
  scoring_mat.index = scoring_mat.columns.values.tolist()
  return scoring_mat


def process_seq_file_list(file_name):
    """
    Given a list of sequence pair filepaths, extract sequences
    and return list
    set up like: 'sequences/prot-0004.fa sequences/prot-0008.fa'
    """
    lines = open(file_name).read().splitlines()
    seqs = []
    for line in open(file_name):
        seq1_name, seq2_name = line.split(" ")[0].rstrip(), line.rsplit(" ")[1].rstrip()
        # read and extract sequences
        seq1 = read_seq(seq1_name)
        seq2 = read_seq(seq2_name)
        seq1, seq2 = seq1.replace("x","*"), seq2.replace("x","*")
        seqs.append((seq1, seq2))
    return seqs

def randomize_matrix(base_mat,i,j):
    '''
    This funciton is going to grab a random value between -10 and 10 and put it into the matrix. The SW algo will then run
    INPUT: Input matrix, i&j are the spots on the matrix you want to change
    OUTPUT: New matrix to test SW algorithm on.
    '''
    new_score = np.random.choice(range(-10,10))
    base_mat.iloc[i,j]=new_score
    base_mat.iloc[j,i]=new_score
    return(base_mat)

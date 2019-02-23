#!/usr/bin/env python
#Algo PS #3
#Stephanie Wankowicz


import os
import pandas as pd
import pickle
import math
import numpy as np
import itertools


def find_match_score(i, j, seq1, seq2, scoring_matrix):
    '''
    This script is to find the value in the scoring matrix of interest.
    '''
    value1 = seq1[i]
    value2 = seq2[j]
    return float(scoring_matrix[value1][value2])

def build_matrix(seq1, seq2, gap_cost, gap_ext, scoring_matrix):
    """
    This function is going to create three matrices, one with match scores, one with gaps in seq1, one with gaps in seq2
    INPUTS:
    seq1: sequence 1 to be aligned
    seq2: sequence 2 to be aligned
    gap_cost: cost for opening a gap
    gap_ext: cost for extending a gap
    scoring_matrix: scoring matrix to generate alignment #this is already read in
    OUTPUT: names of aligned sequences and alignment score
    """
    #create empty matrix that is one larger than seq 1&2 (due to stop score)
    lower_matrix = np.full((len(seq1) + 1, len(seq2) + 1), float("-inf")) #gaps for first string
    middle_matrix = np.full((len(seq1) + 1, len(seq2) + 1), float("-inf")) #match/mismatch score
    upper_matrix = np.full((len(seq1) + 1, len(seq2) + 1), float("-inf")) #gaps for second string
    trace_matrix = np.zeros((len(seq1) + 1, len(seq2) + 1), np.int) #this is to keep track of where you are going.

    max_score = 0
    max_pos = None

    for i, j in itertools.product(range(1, len(seq1)+1), range(1, len(seq2)+1)): #iterate over each AA in each sequence
        match_score = find_match_score(i-1, j-1, seq1, seq2, scoring_matrix) #find the match and/or mismatch score for the item diagonally below the cell we are trying to fill in now
        #for each matrix, decide which value is the highest
        lower_matrix[i][j] = max(lower_matrix[i-1][j] - gap_ext, middle_matrix[i-1][j] - gap_cost - gap_ext) #gaps in second seq
        upper_matrix[i][j] = max(upper_matrix[i][j-1] - gap_ext, middle_matrix[i][j-1] - gap_cost - gap_ext) #gaps in first seq
        middle_matrix[i][j] = max(0, lower_matrix[i][j], upper_matrix[i][j], middle_matrix[i-1][j-1] + match_score)
        if middle_matrix[i][j] > max_score: #we are keeping track of the maximum score (used for all the questions). 
            max_score = middle_matrix[i][j] # filling in the new maximum score
            max_pos = i, j #this is going to tell us where we want to start tracing back
        if middle_matrix[i][j] == middle_matrix[i-1][j-1] + match_score or middle_matrix[i][j] == 0: #if we moved diagonal or are starting back at a new substring, indicate staying on middle
            trace_matrix[i][j]= 1
        elif middle_matrix[i][j] == upper_matrix[i][j]: #else if we are using a gap in the first sequence
            trace_matrix[i][j] = 2
        else: #else if we are using a gap in the second sequence
            trace_matrix[i][j] = 3
    return trace_matrix, max_pos, max_score, middle_matrix

#testing
#trace,start_pos,max_score,main_bitch=build_matrix(seq1, seq2, 2, 1, scoring_mat)


#figure out which
def traceback(trace_matrix, start_pos, seq1, seq2, main_matrix):
    '''
    This script will allow us to move along our trace matrix and output the sequences alignment for each sequence.
    '''
    align1 = [] #start off with blank lists
    align2 = []
    x, y = start_pos #start position 
    while main_matrix[x][y] != 0: #we are starting from the highest position (start position), and moving down to when the match/mismatch is 0
        if trace_matrix[x][y] == 1: #meaning we are staying in the main matrix (ie there was a match/mismatch)
            x -= 1 #move diagoinally down
            y -= 1 #move diagoinally down
            align1 += seq1[x] #add value from both sequences to output
            align2 += seq2[y] #add value from both sequences to output

        if trace_matrix[x][y] == 2: #meaning we had a gap in seq1
            y -= 1 #only move seq2 down
            align1 += '-' #indicate gap in first sequence
            align2 += seq2[y]#add value from second sequences to output

        if trace_matrix[x][y] == 3: #meaning we had a gap in seq2
            x -= 1 #only move seq2 down
            align1 += seq1[x] #add value from first sequences to output
            align2 += '-'#indicate gap in second sequence
    return ''.join(align1[::-1]), ''.join(align2[::-1]) #return each both reverse sequence

#print(traceback(trace,start_pos,seq1,seq2,main_bitch))

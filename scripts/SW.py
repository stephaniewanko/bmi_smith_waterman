#!/usr/bin/env python
#Algo PS #3
#Stephanie Wankowicz


import os
import pandas as pd
import pickle
import math
import numpy as np
import itertools

'''
Note: testing function- does the same sequence get optimal scoring function

'''

#os.chdir('/Users/stephaniewankowicz/Dropbox/BMI_203/HW3_due_02_23/')

#scoring_mat=read_scoring_mat('PAM100')

#os.chdir('/Users/stephaniewankowicz/Dropbox/BMI_203/HW3_due_02_23/sequences/')
#seq1=read_seq('prot-0069.fa')
#seq2=read_seq('prot-0008.fa')


def find_match_score(i, j, seq1, seq2, scoring_matrix):
    '''
    This function is looking up the match/mismatch score for a pair of AA as indicated by the i/j character of seq1/se2
    INPUT:
    OUTPUT:
    '''
    value1 = seq1[i]
    value2 = seq2[j]
    return float(scoring_matrix[value1][value2])

def build_matrix(seq1, seq2, gap_cost, gap_ext, scoring_matrix):
    """
    This function is going to create three matrices
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

    for i, j in itertools.product(range(1, len(seq1)+1), range(1, len(seq2)+1)): #iterate over each sequence
        match_score = find_match_score(i-1, j-1, seq1, seq2, scoring_matrix) #find the match and/or mismatch score for the item diagonally below the cell we are trying to fill in now
        #for each matrix, decide which value is the lowest
        lower_matrix[i][j] = max(lower_matrix[i-1][j] - gap_ext, middle_matrix[i-1][j] - gap_cost - gap_ext)
        upper_matrix[i][j] = max(upper_matrix[i][j-1] - gap_ext, middle_matrix[i][j-1] - gap_cost - gap_ext)
        middle_matrix[i][j] = max(0, lower_matrix[i][j], upper_matrix[i][j], middle_matrix[i-1][j-1] + match_score)
        if middle_matrix[i][j] > max_score:
            max_score = middle_matrix[i][j] # filling in the new maximum score
            max_pos = i, j #this is going to tell us where we want to start tracing back
        if middle_matrix[i][j] == middle_matrix[i-1][j-1] + match_score or middle_matrix[i][j] == 0: #if we moved diagonal or are starting back at a new substring, indicate staying on middle
            trace_matrix[i][j]= 1
        elif middle_matrix[i][j] == upper_matrix[i][j]: #else if we are
            trace_matrix[i][j] = 2
        else:
            trace_matrix[i][j] = 3
    return trace_matrix, max_pos, max_score, middle_matrix


#trace,start_pos,max_score,main_bitch=build_matrix(seq1, seq2, 2, 1, scoring_mat)


#figure out which
def traceback(trace_matrix, start_pos, seq1, seq2, main_matrix):
    '''
    INPUT: Position= max position, where we want to star the alignment
    '''
    align1 = []
    align2 = []
    x, y = start_pos #start position
    while main_matrix[x][y] != 0: #we are starting from the highest position, and moving down to when the match/mismatch is 0
        if trace_matrix[x][y] == 1: #meaning we are staying in the main matrix (ie there was a match/mismatch)
            x -= 1
            y -= 1
            align1 += seq1[x]
            align2 += seq2[y]

        if trace_matrix[x][y] == 2: #meaning we had a gap in seq1
            y -= 1
            align1 += '-'
            align2 += seq2[y]

        if trace_matrix[x][y] == 3: #meaning we had a gap in seq2
            x -= 1
            align1 += seq1[x]
            align2 += '-'
    return ''.join(align1[::-1]), ''.join(align2[::-1]) #return each both reverse sequence

#print(traceback(trace,start_pos,seq1,seq2,main_bitch))

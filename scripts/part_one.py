#!/usr/bin/env python
#Algo PS #3
#Stephanie Wankowicz

##ANSWERING PART ONE

from SW import *
from file_processing import *
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#load in scoring matrix # BLOM50
scoring_matrix = read_scoring_mat('../BLOSUM50')

#load in sequences
positive = process_seq_file_list('/Users/stephaniewankowicz/Dropbox/BMI_203/HW3_due_02_23_GIT/Pospairs.txt')
negative = process_seq_file_list('/Users/stephaniewankowicz/Dropbox/BMI_203/HW3_due_02_23_GIT/Negpairs.txt')


	# Create alignments for each positive and negative pair of alignments
#	run_alignment(pos_seqs, score_mats.get(best_mat), out_align=True, outfile="alignments_pos.txt")
#	run_alignment(neg_seqs, score_mats.get(best_mat), out_align=True, outfile="alignments_neg.txt")



'''
Question 1a:
Consider the false positive rate (proportion of negative pairs with scores that exceed a scorethreshold)
when the true positive rate (proportion of positive pairs with scores above thethreshold) is 0.7.
What's the best false positive rate that you can achieve with varying both gap opening (from 1 to 20)
and extension penalties (from 1 to 5) with the BLOSUM50matrix?
'''

'''
Quesiton 1b:
What is the best gap penalty combination gap opening (from 1 to 20) and extension penalties (from 1 to 5) with the BLOSUM50matrix?
'''

def vary_gap_penalties(pos_fa, neg_fa,max_gap, max_ext,score_mat):
    '''
    Function to determine the best gap and extension penalities by finding the
    smallest fpr with tpr of 0.7
    INPUT: positive and negative sequences, max gap and ext, dict of scores
    OUTPUT: file with optimal gap and extension
    '''
    fpr=100
    best_gap=0
    best_ext=0
    positive_score_list = np.zeros(len(pos_fa))
    negative_score_list = np.zeros(len(neg_fa))
    for i,j in itertools.product(range(1, max_gap+1), range(1, max_ext+1)): #i==max_gap  #j==max_ext
        print(i,j)
        positive_score_list = np.zeros(len(pos_fa))
        negative_score_list = np.zeros(len(neg_fa))
        #for each pair of gap and ext values, iterate over each pair of negative and positive sequences. Return the alignment score
        count=0
        for seq in pos_fa: #list of 2 sequences that go together
            trace,start_pos,max_score,main_bitch=build_matrix(seq[0], seq[1], i,j,score_mat)
            positive_score_list[count] = max_score
            count+=1
            #np.append(positive_score_list,max_score)
        count=0
        for seq in neg_fa:
            trace,start_pos,max_score,main_bitch=build_matrix(seq[0], seq[1], i,j,score_mat)
            negative_score_list[count] = max_score
            count+=1
            #negative_score_list=np.append(negative_score_list,max_score)
        #calculate TPR + FPR
        #sort positive_score_list
        sorted_pos_scores=np.sort(positive_score_list)[::-1]
        tpr_cutoff=0.7*len(sorted_pos_scores) #we want to keep the TPR at 0.7. Determine score cutoff at that level.
        score_threshold=sorted_pos_scores[int(tpr_cutoff)] #determine the score threshold
        tmp_fpr=len(negative_score_list[negative_score_list>=score_threshold])/len(negative_score_list) #determine the numebr of neg sequences above the score threshold
        if tmp_fpr<fpr: #keep track of the best fpr, gap, ext costs to output
            print('switch!')
            fpr=tmp_fpr
            best_gap=i
            best_ext=j
            print(fpr)
    print(fpr)
    print('Best Gap:')
    print(best_gap)
    print('Best Extension:')
    print(best_ext)
    return fpr, best_gap, best_ext
#vary_gap_penalties(positive, negative,20, 5,scoring_matrix)


'''
#Question 2a: Using the gap penalties you determined from question 1, which of the provided scoring matrices performs the best,
#in terms of false positive rate (at a true positive rate of 0.7)?
'''
def pick_best_matrix():
    matrix_list=['../BLOSUM50', '../BLOSUM62', '../MATIO', '../PAM100', '../PAM250']
    pos_max_scores=[]
    neg_max_scores=[]
    positive_scores={}
    negative_scores={}
    for i in matrix_list: #iterate over every matrix
        print(i)
        positive_score_list = np.zeros(len(positive))
        negative_score_list = np.zeros(len(negative))
        matrix=read_scoring_mat(i)
        count=0
        for j in positive:
            trace,start_pos,max_score,main_bitch=build_matrix(j[0], j[1], 4,3,matrix) #we found out the best gap=5, ext=3
            #traceback if you want to return the alignments
            #=traceback(trace_matrix, start_pos, seq1, seq2, main_matrix)
            positive_score_list[count] = max_score
            positive_scores[i]=pos_max_scores
            count+=1
        count=0
        for j in negative:
            trace,start_pos,max_score,main_bitch=build_matrix(j[0], j[1], 4,3,matrix) #we found out the best gap=5, ext=3
            negative_score_list[count] = max_score
            count+=1
        sorted_pos_scores=np.sort(positive_score_list)[::-1]
        tpr_cutoff=0.7*len(sorted_pos_scores) #we want to keep the TPR at 0.7. Determine score cutoff at that level.
        score_threshold=sorted_pos_scores[int(tpr_cutoff)] #determine the score threshold
        fpr=len(negative_score_list[negative_score_list>=score_threshold])/len(negative_score_list) #determine the numebr of neg sequences above the score threshold
        print(fpr)

#pick_best_matrix()
'''
#Question 2b: What are the performance rates of each of the matrices? (create ROC curve)
'''

def performance_matrix():
    matrix_list=['BLOSUM50' , 'BLOSUM62', 'MATIO', 'PAM100', 'PAM250']
    for i in matrix_list: #iterate over every matrix
        print(i)
        #positive_score_list = np.zeros(len(pos_fa))
        #negative_score_list = np.zeros(len(neg_fa))
        negative_score_df=pd.DataFrame()
        positive_score_df=pd.DataFrame()
        mat_name=('../'+i)
        print(mat_name)
        matrix=read_scoring_mat(mat_name)
        count=0
        for j in positive:
            trace,start_pos,max_score,main_bitch = build_matrix(j[0], j[1], 4,3,matrix) #we found out the best gap=5, ext=3
            positive_score_df.loc[count, 'Max_Score'] = max_score
            positive_score_df.loc[count, 'Pos_Neg'] = 'Positive'
            count+=1
        count=0
        for j in negative:
            trace,start_pos,max_score,main_bitch = build_matrix(j[0], j[1], 4,3,matrix) #we found out the best gap=5, ext=3
            negative_score_df.loc[count, 'Max_Score'] = max_score
            negative_score_df.loc[count, 'Pos_Neg'] = 'Negative'
            count+=1
        #combine the df
        combined = positive_score_df.append(pd.DataFrame(data = negative_score_df), ignore_index=True)
        combined_sorted=combined.sort_values(['Max_Score'], ascending=False, )
        combined_sorted=combined_sorted.reset_index()
        for index, row in combined_sorted.iterrows():
            subset=combined_sorted[0:(index+1)]
            #len(combined_sorted.loc[combined_sorted['Pos_Neg'] == 'Positive'].index)
            #print(combined_sorted.loc[,['Pos_Neg'] == 'Positive'])
            combined_sorted.loc[index,'TP']=len(subset.loc[subset['Pos_Neg'] == 'Positive'].index)/(len(combined_sorted.index))
            combined_sorted.loc[index,'FP']=len(subset.loc[subset['Pos_Neg'] == 'Negative'].index)/(len(combined_sorted.index))
            combined_sorted.loc[index,'FN']=(len(positive)-len(subset.loc[subset['Pos_Neg'] == 'Positive'].index))/(len(combined_sorted.index))
            combined_sorted.loc[index,'TN']=(len(negative)-len(subset.loc[subset['Pos_Neg'] == 'Negative'].index))/(len(combined_sorted.index))
        combined_sorted['FPR']=combined_sorted['FP']/(combined_sorted['FP']+combined_sorted['TN'])
        combined_sorted['TPR']=combined_sorted['TP']/(combined_sorted['TP']+combined_sorted['FN'])
        combined_sorted.to_csv('combined_sorted.csv', sep=',')
        #print(type(combined_sorted['FP'].values))
        output_name='combined_score_'+i
        output_name=combined_sorted
    plt.plot((combined_sorted_BLOSUM50['FPR'].values),combined_sorted_BLOSUM50['TPR'].values)
    plt.plot((combined_sorted_BLOSUM62['FPR'].values),combined_sorted_BLOSUM62['TPR'].values)
    plt.plot((combined_sorted_MATIO['FPR'].values),combined_sorted_MATIO['TPR'].values)
    plt.plot((combined_sorted_PAM100['FPR'].values),combined_sorted_PAM100['TPR'].values)
    plt.plot((combined_sorted_PAM250['FPR'].values),combined_sorted_PAM250['TPR'].values)
    plt.plot([0,1],[0,1], 'k--') #Identity Line
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC Curves')
    plt.xlim(0,1)
    plt.ylim(0,1)
    plt.axes().set_aspect('equal')
    plt.show()
        #combined_sort['TP']=/combined_sort.index
        #combined_sorted.to_csv(i + "_ROCinput.txt")
    return combined_sorted
performance_matrix()

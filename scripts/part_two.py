#!/usr/bin/env python
#Algo PS #3
#Stephanie Wankowicz

##ANSWERING PART TWO

from SW import *
from file_processing import *
import numpy as np
import math
import matplotlib.pyplot as plt

scoring_matrix = read_scoring_mat('../BLOSUM62') #best matrix

#load in sequences
positive = process_seq_file_list('/Users/stephaniewankowicz/Dropbox/BMI_203/HW3_due_02_23_GIT/Pospairs.txt')
negative = process_seq_file_list('/Users/stephaniewankowicz/Dropbox/BMI_203/HW3_due_02_23_GIT/Negpairs.txt')

'''
Pre-Questions: Create an alignment for each positive pair of sequences and each negative pair of sequences
'''

def return_alignment():
    with open("Part2_BLOSUM62_seqoutput.txt", "w") as f:
        for seq in positive: #list of 2 sequences that go together
            trace,start_pos,max_score,main_bitch=build_matrix(seq[0], seq[1], 5,3,scoring_matrix)
            seq1,seq2=traceback(trace, start_pos, seq[0], seq[1], main_bitch)
            f.write("{} {} {} Positive\n".format(seq1,seq2,max_score))
        for seq in negative: #list of 2 sequences that go together
            trace,start_pos,max_score,main_bitch=build_matrix(seq[0], seq[1], 5,3,scoring_matrix)
            seq1,seq2=traceback(trace, start_pos, seq[0], seq[1], main_bitch)
            f.write("{} {} {} Negative\n".format(seq1,seq2,max_score))
return_alignment()

'''
Question 1a: Devise an optimization algorithm to modify the values in a starting score matrix
such as to maximize the following objective function: sum of TP rates for FP rates of 0.0, 0.1, 0.2, and 0.3.
'''

def calc_tpr(fpr, positive_score_list, negative_score_list):
    """
    This is going to calculate the TPR of the positive list, given a FPR
    INPUT: FPR cutoff you want to use, List of scores from positive controls, list of scores from negative controls
    OUTPUT: TPR
    """
    sorted_neg_scores=np.sort(negative_score_list)[::-1]
    fpr_cutoff=fpr*len(sorted_neg_scores) #we want to keep the TPR at 0.7. Determine score cutoff at that level.
    score_threshold=sorted_neg_scores[int(fpr_cutoff)] #determine the score threshold
    tpr=len(positive_score_list[positive_score_list>=score_threshold])/len(positive_score_list) #determine the numebr of neg sequences above the score threshold
    return tpr

def objective_function(positive_score_list, negative_score_list):
    '''
    This is going to calculate the objective function, as defined in the assignment.
    INPUT: List of scores from positive controls, list of scores from negative controls
    OUTPUT: the optomization function
    '''
    opt_fun=[]
    for i in [0,0.1,0.2,0.3]:
        opt_fun.append(calc_tpr(i,positive_score_list, negative_score_list))
    return sum(opt_fun)


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

#for static alignment
def calc_total_score(seq, scoring_mat, gap, ext):
    '''
    INPUT: Pair of sequences, scoring matrix, gap + ext get_values
    OUTPUT: Score to be put in optimization function
    '''
    seq1 = seq[0]
    seq2 = seq[1]
    score = 0
    for n in range(len(seq1)):
        if seq1[n] != '-' and seq2[n] != '-':
            score += float(scoring_mat[seq1[n]][seq2[n]])
        else:
            if seq1[n] == '-':
                if seq1[n-1] == '-':
                    score -= ext
                else:
                    score -= gap
            if seq2[n] == '-':
                if seq2[n-1] == '-':
                    score -= ext
                else:
                    score -= gap
    return float(score)

def plot_ROC(positive_df, negative_df, label):
        combined = positive_df.append(pd.DataFrame(data = negative_df), ignore_index=True)
        combined_sorted=combined.sort_values(['Max_Score'], ascending=False)
        combined_sorted=combined_sorted.reset_index()
        print(combined_sorted)
        for index, row in combined_sorted.iterrows():
            subset=combined_sorted[0:(index+1)]
            combined_sorted.loc[index,'TP']=len(subset.loc[subset['Pos_Neg'] == 'Positive'].index)/(len(combined_sorted.index))
            combined_sorted.loc[index,'FP']=len(subset.loc[subset['Pos_Neg'] == 'Negative'].index)/(len(combined_sorted.index))
            combined_sorted.loc[index,'FN']=(len(positive)-len(subset.loc[subset['Pos_Neg'] == 'Positive'].index))/(len(combined_sorted.index))
            combined_sorted.loc[index,'TN']=(len(negative)-len(subset.loc[subset['Pos_Neg'] == 'Negative'].index))/(len(combined_sorted.index))
        combined_sorted['FPR']=combined_sorted['FP']/(combined_sorted['FP']+combined_sorted['TN'])
        combined_sorted['TPR']=combined_sorted['TP']/(combined_sorted['TP']+combined_sorted['FN'])
        print(combined_sorted.head())
        #combined_sorted.to_csv('combined_sorted.csv', sep=',')
        plt.plot(combined_sorted['FPR'].values,combined_sorted['TPR'].values, label=label)



def optimize_static(scoring_matrix):
    #get initial scoring value
    starting_mat=scoring_matrix
    count=0
    positive_score_list = np.zeros(len(positive))
    negative_score_list = np.zeros(len(negative))
    positive_seq={}
    negative_seq={}
    negative_score_df_original=pd.DataFrame()
    positive_score_df_original=pd.DataFrame()
    for seq in positive: #list of 2 sequences that go together
        trace,start_pos,max_score,main_bitch=build_matrix(seq[0], seq[1],4,3,scoring_matrix)
        positive_score_list[count] = max_score
        seq1,seq2=traceback(trace, start_pos, seq[0], seq[1], main_bitch)
        positive_seq[count] = np.array([seq1,seq2])
        positive_score_df_original.loc[count, 'Max_Score'] = max_score
        positive_score_df_original.loc[count, 'Pos_Neg'] = 'Positive'
        count+=1
    count=0
    for seq in negative:
        trace,start_pos,max_score,main_bitch=build_matrix(seq[0], seq[1],4,3,scoring_matrix)
        negative_score_list[count] = max_score
        seq1,seq2=traceback(trace, start_pos, seq[0], seq[1], main_bitch)
        negative_seq[count] = np.array([seq1,seq2])
        negative_score_df_original.loc[count, 'Max_Score'] = max_score
        negative_score_df_original.loc[count, 'Pos_Neg'] = 'Negative'
        count+=1
    obj_fun=objective_function(positive_score_list, negative_score_list)
    print('original object function')
    print(obj_fun)
    iterations=0
    while(iterations<500): #iterate 10 times
        #first, mutate the matrix
        i=np.random.choice(range(len(starting_mat.columns)))
        j=np.random.choice(range(len(starting_mat.columns)))
        iterations+=1
        testing_matrix=randomize_matrix(starting_mat,i,j)
        testing_positive_list_score=np.zeros(len(positive))
        testing_negative_list_score=np.zeros(len(negative))
        for key, value in positive_seq.items():
            testing_positive_list_score[key]=calc_total_score(value, testing_matrix,4,3)
        for key, value in negative_seq.items():
            testing_negative_list_score[key]=calc_total_score(value, testing_matrix,4,3)
        tmp_obj_fun=objective_function(testing_positive_list_score, testing_negative_list_score)
        if tmp_obj_fun>obj_fun:#if the objective function is better than before, keep the change, else return to previous.
            print('switch!')
            print(tmp_obj_fun)
            obj_fun=tmp_obj_fun
            starting_mat=testing_matrix
    #get values for ROC
    negative_score_df=pd.DataFrame()
    positive_score_df=pd.DataFrame()
    for key, value in positive_seq.items():
        positive_score_df.loc[key, 'Max_Score'] = calc_total_score(value, starting_mat,4,3)
        positive_score_df.loc[key, 'Pos_Neg'] = 'Positive'
    for key, value in negative_seq.items():
        negative_score_df.loc[key, 'Max_Score'] = calc_total_score(value, starting_mat,4,3)
        negative_score_df.loc[key, 'Pos_Neg'] = 'Negative'
    combined = positive_score_df.append(pd.DataFrame(data = negative_score_df), ignore_index=True)
    combined_sorted=combined.sort_values(['Max_Score'], ascending=False)
    combined_sorted=combined_sorted.reset_index()
    for index, row in combined_sorted.iterrows():
        subset=combined_sorted[0:(index+1)]
        combined_sorted.loc[index,'TP']=len(subset.loc[subset['Pos_Neg'] == 'Positive'].index)/(len(combined_sorted.index))
        combined_sorted.loc[index,'FP']=len(subset.loc[subset['Pos_Neg'] == 'Negative'].index)/(len(combined_sorted.index))
        combined_sorted.loc[index,'FN']=(len(positive)-len(subset.loc[subset['Pos_Neg'] == 'Positive'].index))/(len(combined_sorted.index))
        combined_sorted.loc[index,'TN']=(len(negative)-len(subset.loc[subset['Pos_Neg'] == 'Negative'].index))/(len(combined_sorted.index))
    combined_sorted['FPR']=combined_sorted['FP']/(combined_sorted['FP']+combined_sorted['TN'])
    combined_sorted['TPR']=combined_sorted['TP']/(combined_sorted['TP']+combined_sorted['FN'])
    combined_org = positive_score_df_original.append(pd.DataFrame(data = negative_score_df_original), ignore_index=True)
    combined_sorted_org=combined_org.sort_values(['Max_Score'], ascending=False)
    combined_sorted_org=combined_sorted_org.reset_index()
    for index, row in combined_sorted_org.iterrows():
        subset=combined_sorted_org[0:(index+1)]
        combined_sorted_org.loc[index,'TP']=len(subset.loc[subset['Pos_Neg'] == 'Positive'].index)/(len(combined_sorted_org.index))
        combined_sorted_org.loc[index,'FP']=len(subset.loc[subset['Pos_Neg'] == 'Negative'].index)/(len(combined_sorted_org.index))
        combined_sorted_org.loc[index,'FN']=(len(positive)-len(subset.loc[subset['Pos_Neg'] == 'Positive'].index))/(len(combined_sorted_org.index))
        combined_sorted_org.loc[index,'TN']=(len(negative)-len(subset.loc[subset['Pos_Neg'] == 'Negative'].index))/(len(combined_sorted_org.index))
    combined_sorted_org['FPR']=combined_sorted_org['FP']/(combined_sorted_org['FP']+combined_sorted_org['TN'])
    combined_sorted_org['TPR']=combined_sorted_org['TP']/(combined_sorted_org['TP']+combined_sorted_org['FN'])
    return testing_matrix, obj_fun


#new_matrix, obj_fun=optimize_static(scoring_matrix)


def optomize_realign(pos,neg, original_matrix, new_matrix):
    count=0
    negative_score_df_original=pd.DataFrame()
    positive_score_df_original=pd.DataFrame()
    for seq in positive: #list of 2 sequences that go together
        trace,start_pos,max_score,main_bitch=build_matrix(seq[0], seq[1],4,3,original_matrix)
        positive_score_df_original.loc[count, 'Max_Score'] = max_score
        positive_score_df_original.loc[count, 'Pos_Neg'] = 'Positive'
        count+=1
    count=0
    for seq in negative:
        trace,start_pos,max_score,main_bitch=build_matrix(seq[0], seq[1],4,3,original_matrix)
        negative_score_df_original.loc[count, 'Max_Score'] = max_score
        negative_score_df_original.loc[count, 'Pos_Neg'] = 'Negative'
        count+=1
    ### using new matrix
    negative_score_df=pd.DataFrame()
    positive_score_df=pd.DataFrame()
    for seq in positive: #list of 2 sequences that go together
        trace,start_pos,max_score,main_bitch=build_matrix(seq[0], seq[1],4,3,new_matrix)
        positive_score_df.loc[count, 'Max_Score'] = max_score
        positive_score_df.loc[count, 'Pos_Neg'] = 'Positive'
        count+=1
    count=0
    for seq in negative:
        trace,start_pos,max_score,main_bitch=build_matrix(seq[0], seq[1],4,3,new_matrix)
        negative_score_df.loc[count, 'Max_Score'] = max_score
        negative_score_df.loc[count, 'Pos_Neg'] = 'Negative'
        count+=1
    combined = positive_score_df.append(pd.DataFrame(data = negative_score_df), ignore_index=True)
    combined_sorted=combined.sort_values(['Max_Score'], ascending=False)
    combined_sorted=combined_sorted.reset_index()
    for index, row in combined_sorted.iterrows():
        subset=combined_sorted[0:(index+1)]
        combined_sorted.loc[index,'TP']=len(subset.loc[subset['Pos_Neg'] == 'Positive'].index)/(len(combined_sorted.index))
        combined_sorted.loc[index,'FP']=len(subset.loc[subset['Pos_Neg'] == 'Negative'].index)/(len(combined_sorted.index))
        combined_sorted.loc[index,'FN']=(len(positive)-len(subset.loc[subset['Pos_Neg'] == 'Positive'].index))/(len(combined_sorted.index))
        combined_sorted.loc[index,'TN']=(len(negative)-len(subset.loc[subset['Pos_Neg'] == 'Negative'].index))/(len(combined_sorted.index))
    combined_sorted['FPR']=combined_sorted['FP']/(combined_sorted['FP']+combined_sorted['TN'])
    combined_sorted['TPR']=combined_sorted['TP']/(combined_sorted['TP']+combined_sorted['FN'])
    plt.plot(combined_sorted['FPR'].values,combined_sorted['TPR'].values, label='Optomized Realign')
    combined_org = positive_score_df_original.append(pd.DataFrame(data = negative_score_df_original), ignore_index=True)
    combined_sorted_org=combined_org.sort_values(['Max_Score'], ascending=False)
    combined_sorted_org=combined_sorted_org.reset_index()
    for index, row in combined_sorted_org.iterrows():
        subset=combined_sorted_org[0:(index+1)]
        combined_sorted_org.loc[index,'TP']=len(subset.loc[subset['Pos_Neg'] == 'Positive'].index)/(len(combined_sorted_org.index))
        combined_sorted_org.loc[index,'FP']=len(subset.loc[subset['Pos_Neg'] == 'Negative'].index)/(len(combined_sorted_org.index))
        combined_sorted_org.loc[index,'FN']=(len(positive)-len(subset.loc[subset['Pos_Neg'] == 'Positive'].index))/(len(combined_sorted_org.index))
        combined_sorted_org.loc[index,'TN']=(len(negative)-len(subset.loc[subset['Pos_Neg'] == 'Negative'].index))/(len(combined_sorted_org.index))
    combined_sorted_org['FPR']=combined_sorted_org['FP']/(combined_sorted_org['FP']+combined_sorted_org['TN'])
    combined_sorted_org['TPR']=combined_sorted_org['TP']/(combined_sorted_org['TP']+combined_sorted_org['FN'])
    print(combined_sorted_org)
    plt.plot(combined_sorted_org['FPR'].values,combined_sorted_org['TPR'].values, label='Original Realign')
    plt.plot([0,1],[0,1], 'k--') #Identity Line
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.xlim(0,1)
    plt.legend()
    plt.ylim(0,1)
    plt.axes().set_aspect('equal')
    plt.show()
#optomize_realign(positive,negative, scoring_matrix, new_matrix)

#do it again with MATIO
scoring_matrix_MATIO = read_scoring_mat('../MATIO') #best matrix
new_matrix_MATIO, obj_fun=optimize_static(scoring_matrix_MATIO)
optomize_realign(positive,negative, scoring_matrix, new_matrix_MATIO)

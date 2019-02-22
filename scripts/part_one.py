from SW import *
from file_processing import *
import numpy as np

#load in scoring matrix # BLOM50
scoring_matrix = read_scoring_mat('../BLOSUM50')

#load in sequences
positive = process_seq_file_list('/Users/stephaniewankowicz/Dropbox/BMI_203/HW3_due_02_23_GIT/Pospairs.txt')
negative = process_seq_file_list('/Users/stephaniewankowicz/Dropbox/BMI_203/HW3_due_02_23_GIT/Negpairs.txt')


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
    print('Best Gap:')
    print(best_gap)
    print('Best Extension:')
    print(best_ext)
    return fpr, best_gap, best_ext
#vary_gap_penalties(positive, negative,20, 5,scoring_matrix)


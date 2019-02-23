
from SW import *
from file_processing import *

#test alignment
def test_align_score():
    score_mat = read_("score_matrices/BLOSUM50")
    s1 = "ABCDEFG"
    s2 = "CDEFG"
    trace_matrix, max_pos, max_score, middle_matrix=build_matrix(seq1, seq2, gap_cost, gap_ext, scoring_matrix)
    assert max_score==5

 def test_import_function():
    



import sys
sys.path.append('../scripts/')
from scripts import SW
from scripts import file_processing

#test alignment
def test_align_score():
    score_mat = file_processing.read_scoring_mat("../BLOSUM50")
    s1 = "ABCDEFG"
    s2 = "CDEFG"
    trace_matrix, max_pos, max_score, middle_matrix=SW.build_matrix(seq1, seq2, gap_cost, gap_ext, scoring_matrix)
    assert max_score==5

def test_optimization():
    base_mat=file_processing.read_scoring_mat("../BLOSUM62")
    new_matrix = file_processing.randomize_matrix(base_mat,i,j)
    assert new_matrix != base_mat
    
 
assert 1==1

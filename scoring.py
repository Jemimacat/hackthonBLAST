
# Nucleotide scoring function

def score_nt_seq(seq1,seq2):
    score = 0
    for i in range(len(seq1)):
        base1 = seq1[i]
        base2 = seq2[i]
        tmp = score_nt(base1,base2)
        score += tmp
    return score

def score_nt(base1,base2): # scoring function
    score = 0
    pairs = ['AC','GT','CA','TG']
    if base1 == base2:
        score += 2
    elif base1 + base2 in pairs:
        score += -5
    else:
        score += -7
    return score
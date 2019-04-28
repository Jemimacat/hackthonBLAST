
# Nucleotide scoring function

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
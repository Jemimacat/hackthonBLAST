from multiprocessing import Pool, Manager
from scoring import score_nt_seq

def extending_words_positions(word,query_seq,scores,one_query_dict,db_dict,max_gap=5,step=1):
    db_positions = {}
    positions = {}

    len_of_word = len(word)
    num_of_words = len(one_query_dict)
    tag_q = 0

    for q_pos in one_query_dict[word]:
        pos = q_pos
        for (db_record,score) in sorted(scores[word].items(), key=lambda d:d[1], reverse=True):
            if not tag_q in positions.keys():
                positions[tag_q] = {}
            positions[tag_q][pos] = 1
            if not tag_q in db_positions.keys():
                db_positions[tag_q] = {}

            for gene in db_dict[db_record].keys():

                if not gene in db_positions[tag_q].keys():
                    db_positions[tag_q][gene] = {}

                for db_pos in db_dict[db_record][gene]:
                    this_db_pos = db_pos
                    tag_db = 0

                    if not tag_db in db_positions[tag_q][gene].keys(): 
                        db_positions[tag_q][gene][tag_db] = {}
                    db_positions[tag_q][gene][tag_db][this_db_pos] = q_pos

                    while pos < num_of_words -1:
                        next_pos = pos + 1
                        next_word = query_seq[next_pos:next_pos+len_of_word]
                        if not next_word in scores.keys():
                            break   

                        for (next_db_record,next_score) in sorted(scores[next_word].items(), key=lambda d:d[1], reverse=True):
                            if gene in db_dict[next_db_record].keys():
                                db_pos_diff = 999999
                                next_db_pos = -1
                                for new_db_pos in db_dict[next_db_record][gene]:
                                    if abs(new_db_pos-this_db_pos) < db_pos_diff:
                                        db_pos_diff = abs(new_db_pos-this_db_pos)
                                        next_db_pos = new_db_pos
                                if db_pos_diff <= max_gap and next_pos > pos:
                                    positions[tag_q][next_pos] = 1
                                    db_positions[tag_q][gene][tag_db][next_db_pos] = next_pos
                                    break
                                if next_db_pos != -1:
                                    this_db_pos = new_db_pos
                        pos = next_pos
                    
                    pos = q_pos
                    this_db_pos = db_pos

                    while pos > 0:
                        prior_pos = pos - 1
                        prior_word = query_seq[prior_pos:prior_pos+len_of_word]
                        if not prior_word in scores.keys():
                            break  
                        for (prior_db_record,prior_score) in sorted(scores[prior_word].items(),key=lambda d:d[1], reverse=True):
                            if gene in db_dict[prior_db_record].keys():
                                db_pos_diff = 999999
                                prior_db_pos = - 1
                                for new_db_pos in db_dict[prior_db_record][gene]:
                                    if abs(new_db_pos-this_db_pos) < db_pos_diff:
                                        db_pos_diff = abs(new_db_pos-this_db_pos)
                                        prior_db_pos = new_db_pos
                                if db_pos_diff <= max_gap and prior_pos < pos:
                                    positions[tag_q][prior_pos] = 1
                                    db_positions[tag_q][gene][tag_db][prior_db_pos] = prior_pos
                                    break
                                if prior_db_pos != -1:
                                    this_db_pos = new_db_pos
                        pos = prior_pos

                    tag_db += 1
        
        tag_q += 1

    return positions, db_positions


def extending_words(word,query_seq,positions,db_positions,genes):
    segment_hits = {}

    len_of_word = len(word)
    for i in positions.keys():       
        this_db_positions = db_positions[i]
        for gene in this_db_positions.keys():           
            ref_seq = genes[gene]
            for j in this_db_positions[gene].keys():
                db_pos_list = sorted(this_db_positions[gene][j].keys())
                q_pos_list = []
                q_segment = ''
                db_segment = ''
                for k in range(len(db_pos_list)):
                    p = db_pos_list[k]
                    q = this_db_positions[gene][j][p]
                    if k > 0:
                        if p - db_pos_list[k-1] == 1: 
                            if q - q_pos_list[-1] == 1:  ## continue segments
                                q_segment += query_seq[q]
                                db_segment += ref_seq[p]
                            else:  ## insertion
                                for x in range(q_pos_list[-1]+1,q+1):
                                   q_segment += query_seq[x]
                                   db_segment += '-'
                        else:  ## deletion
                            for x in range(db_pos_list[k-1],p+1):
                                q_segment += '-'
                                db_segment += ref_seq[x]
                    q_pos_list.append(q)
                q_segment += query_seq[q_pos_list[-1]+1:q_pos_list[-1]+1+len_of_word]
                db_segment += ref_seq[db_pos_list[-1]+1:db_pos_list[-1]+1+len_of_word]
                layer = ''
                score = score_nt_seq(q_segment,db_segment)
                consensus = 0
                q_length = 0
                for c in range(len(q_segment)):
                    if q_segment[c] != '-':
                        q_length += 1
                    if q_segment[c] == db_segment[c]:
                        consensus += 1
                        layer += ' '
                    else:
                        layer += '+'
                if not q_length in segment_hits.keys():
                    segment_hits[q_length] = {}
                if not consensus in segment_hits[q_length].keys():
                    segment_hits[q_length][consensus] = {}
                if not score in segment_hits[q_length][consensus].keys():
                    segment_hits[q_length][consensus][score] = {}
                if not gene in segment_hits[q_length][consensus][score].keys():
                    segment_hits[q_length][consensus][score][gene] = {}
                if not q_segment in segment_hits[q_length][consensus][score][gene].keys():
                    segment_hits[q_length][consensus][score][gene][q_segment] = {}
                if not db_segment in segment_hits[q_length][consensus][score][gene][q_segment].keys():
                    segment_hits[q_length][consensus][score][gene][q_segment][db_segment] = {}
               
                segment_hits[q_length][consensus][score][gene][q_segment][db_segment][(layer,q_pos_list[0],q_pos_list[-1]+len_of_word-1,db_pos_list[0],db_pos_list[-1]+len_of_word)] = 1

    return segment_hits
                        
                            
                        


                                











# def one_query_seeds_concat(word,query_seq,scores,one_query_dict,db_dict):

#     (db_record,score) = sorted(scores[word].items(), key=lambda d:d[1], reverse=True)[0]
#     db_pos = db_dict[db_record]
#     for q_pos in one_query_dict[word]:

    



    



    


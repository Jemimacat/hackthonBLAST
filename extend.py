from multiprocessing import Pool, Manager

def build_hit_matrix():
    



def xxxxxx(word,query_seq,scores,one_query_dict,db_dict,max_gap=5,step=1):
    concat_seed = {}
    positions = {}

    len_of_word = len(word)
    num_of_words = len(one_query_dict)
    for (db_record,score) in sorted(scores[word].items(), key=lambda d:d[1], reverse=True):
        for gene in db_dict[db_record].keys():
            for db_pos in db_dict[db_record][gene]:
                for q_pos in one_query_dict[word]:
                    pos = q_pos
                    while pos < num_of_words -1:
                        position = find_neiborghood(pos,word,query_seq,gene,db_pos,db_dict,dire=1,max_gap=max_gap,step=step)
                        positions.update(position)
                        pos += 1
                    last = pos

                    pos = q_pos
                    while pos > 0:
                        find_neiborghood(pos,word,query_seq,gene,db_pos,db_dict,dire=0,max_gap=max_gap,step=step)
                        pos -= 1
                    first = pos
                    if not query_seq[first,last+len_of_word] in concat_seed:
                        concat_seed[query_seq[first,last+len_of_word]] = {}
                    if not gene in concat_seed[query_seq[first,last+len_of_word]]:
                        concat_seed[query_seq[first,last+len_of_word]][gene] = []
                    concat_seed[query_seq[first,last+len_of_word]][gene].append((first,last+len_of_word))
    return concat_seed

def find_neiborghood(pos,word,query_seq,gene,db_pos,db_dict,dire=1,max_gap=5,step=1):
    ## find neiborghood hits. dire = 1 for foreward, 0 for backward
    positions = {}

    len_of_word = len(word)
    if dire == 1:
        new_word = query_seq[pos+step,pos+len_of_word+step]
    elif dire == 0:
        new_word = query_seq[pos-step,pos+len_of_word-step]

    for (db_record,score) in sorted(scores[new_word].items(), key=lambda d:d[1], reverse=True):
        if gene in db_dict[db_record]:
            for pos in db_dict[db_record][gene]:
                abs_diff = abs(db_pos - pos)
                if abs_diff <= max_gap:
                    positions[pos] = db_record

    return positions,new_word,
            









def one_query_seeds_concat(word,query_seq,scores,one_query_dict,db_dict):

    (db_record,score) = sorted(scores[word].items(), key=lambda d:d[1], reverse=True)[0]
    db_pos = db_dict[db_record]
    for q_pos in one_query_dict[word]:

    



    



    


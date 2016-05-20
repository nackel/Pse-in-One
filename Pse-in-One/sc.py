# -*- coding: utf-8 -*-
"""
Created on Sat May 13 15:35:42 2016
@version:0.1.2
@author: Nackel
"""

import re
import time
from itertools import combinations_with_replacement, permutations, product
import numpy as np
from data import index_list
from util_sc import get_rnasc_data

    
def get_kmer_lst(letter, k):
    """Generate a list of all possible k-mer pattern.
    
    :param letter: a list that contains all the possible letters in the sequence
    :param k:the length of k-mer
    :return: a kmer list
    """
    kmerlst = []
    letter_set = set(letter)
    letter = [''.join(i) for i in letter_set]
    
    partkmers = list(combinations_with_replacement(letter, k))

    for element in partkmers:
        elelst = set(permutations(element, k))
        strlst = [''.join(ele) for ele in elelst]
        kmerlst += strlst
    kmerlst = np.sort(kmerlst)
    return list(kmerlst)


def delete_free_base(sequence, sstructure):
    """Delete free base based on secondary structure to produce a new sequence and secondary structure. New sequence and secondary structure is a substring of the original sequence and secondary structure.
    :param sequence: an RNA sequence.
    :param sstructure: its corresponding secondary structure.
    :return: a new sequence and sstructure,string.
    """
    left_pos = sstructure.index('(')
    right_pos = sstructure.rindex(')')
    return sequence[left_pos:right_pos+1], sstructure[left_pos:right_pos+1]
    
    
def delete_loop(sequence, sstructure):
    """Delete loop(hairpin) based on secondary structure to produce a new sequence and secondary structure. New sequence and secondary structure is a substring of the original sequence and secondary structure.
   
    :param sequence: an RNA sequence.
    :param sstructure: its corresponding secondary structure.
    :return: a new sequence and sstructure,string.
    """
    loop_re = r'(\(\.+\))'
    loop_list = re.findall(loop_re, sstructure)
    for loop in loop_list:
        pos = sstructure.index(loop)
        length = len(loop)
        sstructure_dl = sstructure[:pos+1] + sstructure[pos+length-1:]
        sequence_dl = sequence[:pos+1] + sequence[pos+length-1:]
    return sequence_dl, sstructure_dl
    
    
#======================Complete process in Triplet=============================
def get_triplet_matrix(filename):
    '''This is a complete process in triplet,aim to gernerate feature vectors.
     
       The FASTA format of the input file is as follows:    
       >sequence name
       An RNA sequence should be consist of AGCU
       Second structure
 
    :param input_file_name: f: HANDLE to input. open(<file>)
    :return: Feature matrix through Triplet
    '''
    letter = ["(","."]
    alphabet = 'AGCU'     #Don't change the alphabetical, or the order of features will change.
    with open(filename) as f:
        seqsslst= get_rnasc_data(f)
    tripletdict = get_triplet_dict(letter, 3, alphabet)
    features = []
    for seqss in seqsslst:
        vector = get_triplet_vector(seqss.sequence,seqss.sstruc,tripletdict)
        features.append(vector)
    return features
 
     
def get_triplet_vector(sequence, sstructure,patterndict):    
    '''This is a process in triplet,aim to gernerate feature vector.
     
    :param sequence: an RNA sequence,string.
    :param sstructure: The corresponding secondary structure, string.
    :param patterndict: All the features, dictionary.
    :return: Feature vector through Triplet.
    '''
    elelen = len(patterndict)
    vector=np.zeros((1,elelen))
    sequence, sstructure = delete_free_base(sequence, sstructure)
    sequence, sstructure = delete_loop(sequence, sstructure)
    
    for i in range(len(sequence)):
        letter =sequence[i]
        middle = sstructure[i]
        if i == 0:
            near_left = "."
            near_right = sstructure[i+1]
        elif i == len(sequence)-1:
            near_left = sstructure[i-1]
            near_right = "."
        else:
            near_left = sstructure[i-1]
            near_right = sstructure[i+1]
        #rectify the empty loop structure
        if middle == '(' and near_right == ')':
            near_right = '.'
        if middle == ')' and near_left == '(':
            near_left = '.'
             
        letter_sstruc_comb = letter+near_left+middle+near_right
        letter_sstruc_comb_r = letter_sstruc_comb.replace(')', '(')
        position=patterndict.get(letter_sstruc_comb_r)
        vector[0, position] += 1
        #print letter_sstruc_comb ,position
    #return list (vector[0])
    return list(vector[0]/sum(vector[0]))    
     
     
def get_triplet_dict(letter, k, alphabet=index_list.RNA):
    """Generate a dictionary of all possible triplet pattern.
    :param letter: a list that contains all the possible characters in the secondary structure. eg:['.','(']
    :param k: the length of k-mer
    :param alphabet: a string that contains all the possible characters in the sequence.
    :return: a triplet dictionary
    """
    kmerlst = get_kmer_lst(letter ,k)
    kmerlst.reverse()
    tripletlst = [''.join(ele) for ele in product(list(alphabet), kmerlst)]
    #tripletlst = np.sort(tripletlst)
    tripletdict = {tripletlst[i]: i for i in range(len(tripletlst))}
    return tripletdict
#==============================================================================

 
def main(args):
    if args.method == "triplet":
        res = get_triplet_matrix(args.inputfile)
    else:
        print("Method error!")
        
     # Write correspond res file.
    if args.f == 'tab':
        from util import write_tab

        write_tab(res, args.outputfile)
    elif args.f == 'svm':
        from util import write_libsvm
        write_libsvm(res, [args.l] * len(res), args.outputfile)
    elif args.f == 'csv':
        from util import write_csv
        write_csv(res, args.outputfile)   
    

if __name__ == '__main__':
#==============================================================================
#     import argparse
#     from argparse import RawTextHelpFormatter
# 
#     parse = argparse.ArgumentParser(description="This is a kmer module for generate kmer vector.",
#                                     formatter_class=RawTextHelpFormatter)
#     parse.add_argument('inputfile',
#                        help="The input file, in valid FASTA format.")
#     parse.add_argument('outputfile',
#                        help="The outputfile stored results.")
#     parse.add_argument('alphabet', choices=['DNA', 'RNA', 'Protein'],
#                        help="The alphabet of sequences.")
#     parse.add_argument('method', type=str,
#                        help="The method name of structure composition.")
#     #parse.add_argument('k', type=int, choices=range(1, 7),
#      #                  help="The k value of kmer.")
#    # parse.add_argument('alphabet', choices=['DNA', 'RNA', 'PROTEIN'],
#     #                   help="The alphabet of sequences.")
#     #parse.add_argument('-r', default=0, type=int, choices=[1, 0],
#      #                  help="Whether need to reverse complement.\n"
#       #                      "1 means True, 0 means False. (default = 0)")
#     parse.add_argument('-f', default='tab', choices=['tab', 'svm', 'csv'],
#                        help="The output format (default = tab).\n"
#                             "tab -- Simple format, delimited by TAB.\n"
#                             "svm -- The libSVM training data format.\n"
#                             "csv -- The format that can be loaded into a spreadsheet program.")
#     parse.add_argument('-l', default='+1', choices=['+1', '-1'],
#                        help="The libSVM output file label.")
# 
#     args = parse.parse_args()
# 
#     print(args)
#     print("Calculating...")
#     start_time = time.time()
#     main(args)
#     
#     print("Used time: %ss" % (time.time() - start_time))
#     print("Done.")
#==============================================================================

#==========================Triplet test========================================
#    letter = ['(', '.']
#    alphabet ="AGCU"
#    sequence = 'CUUUCUACACAGGUUGGGAUCGGUUGCAAUGCUGUGUUUCUGUAUGGUAUUGCACUUGUCCCGGCCUGUUGAGUUUGG'
#    sstructure="..(((...((((((((((((.(((.(((((((((((......)))))))))))))).)))))))))))).)))....."
#    patterndic= get_triplet_dict(letter, 3, alphabet)
#    vector =get_triplet_vector(sequence, sstructure, patterndic)
#    lst=[">hsa-let-7c MI0000064", 'CUUUCUACACAGGUUGGGAUCGGUUGCAAUGCUGUGUUUCUGUAUGGUAUUGCACUUGUCCCGGCCUGUUGAGUUUGG', '..(((...((((((((((((.(((.(((((((((((......)))))))))))))).)))))))))))).))).....']
#    is_rnasc_list(lst)
#==============================================================================
    list_pattern = ['C', 'G', 'U', 'A-U', 'U-A', 'G-C', 'C-G', 'G-U', 'U-G']
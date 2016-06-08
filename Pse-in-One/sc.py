# -*- coding: utf-8 -*-
"""
Created on Sat May 13 15:35:42 2016
@version:0.2.1
@author: Nackel
"""
import sys
import re
import time
from itertools import combinations_with_replacement, permutations, product
import numpy as np
from data import index_list
from util_sc import is_rnasc_list
from util_sc import get_rnasc_data
from util_sc import get_corresp_sequence
import os
    
def get_kmer_lst(letterlst, k):
    """Generate a list of all possible k-mer pattern.
    
    :param letter: a list that contains all the possible letters in the sequence
    :param k:the length of k-mer
    :return: a kmer list
    """
    kmerlst = []
    letter_set = set(letterlst)
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
    :param sequence: An RNA sequence.
    :param sstructure: Its corresponding secondary structure.
    :return: A new sequence and sstructure,string.
    """
    left_pos = sstructure.index('(')
    right_pos = sstructure.rindex(')')
    return sequence[left_pos:right_pos+1], sstructure[left_pos:right_pos+1]
    
    
def delete_loop(sequence, sstructure):
    """Delete loop(hairpin) based on secondary structure to produce a new sequence and secondary structure. New sequence and secondary structure is a substring of the original sequence and secondary structure.
   
    :param sequence: An RNA sequence.
    :param sstructure: Its corresponding secondary structure.
    :return: A new sequence and sstructure,string.
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
 
    :param filename: Name of inputfile.
    :return: Feature matrix through Triplet.
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
    vector=np.zeros((1,len(patterndict)))
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
        position = patterndict.get(letter_sstruc_comb_r)
        vector[0, position] += 1
        #print letter_sstruc_comb ,position
    #return list (vector[0])
    return [round(f,3) for f in list(vector[0]/sum(vector[0]))]
     
     
def get_triplet_dict(letter, k, alphabet=index_list.RNA):
    """Generate a dictionary of all possible triplet pattern.
    :param letter: A list that contains all the possible characters in the secondary structure. eg:['.','(']
    :param k: The length of k-mer
    :param alphabet: A string that contains all the possible characters in the sequence.
    :return: A triplet dictionary.
    """
    kmerlst = get_kmer_lst(letter ,k)
    kmerlst.reverse()
    tripletlst = [''.join(ele) for ele in product(list(alphabet), kmerlst)]
    #tripletlst = np.sort(tripletlst)
    tripletdict = {tripletlst[i]: i for i in range(len(tripletlst))}
    return tripletdict
    
    
#=========================PseKNC===============================================
def get_pseknc_matrix(filename,k):
    '''This is a complete process in PseKNC,aim to gernerate feature matrix.
     
       The FASTA format of the input file is as follows:    
       >sequence name
       An RNA sequence should be consist of AGCU
       Second structure
 
    :param filename: Name of input file.
    :return: Feature matrix through PseKNC.
    '''
    
    alphabet = 'ACGU'
    letter = list(alphabet)
    with open(filename) as f:
        seqsslst= get_rnasc_data(f)
    psekncdict = get_pseknc_dict(letter, k)
    features = []
    for seqss in seqsslst:
        vector = get_pseknc_vector(seqss.sequence,seqss.sstruc,psekncdict, k)
        features.append(vector)
    return features


def get_pseknc_dict(letter_list, k):
    """Generate a dictionary of all possible PseKNC pattern.
    :param letter: A list that contains all the possible characters in an RNA sequence. eg:['A','C','G','U']
    :param k: The length of K-tuple nucleotide composition.
    :return: A PseKNC pattern dictionary.
    """
    pseknclst = []
    part_psessc = list(combinations_with_replacement(letter_list, k))
    for element in part_psessc:
        elelst=set(permutations(element, k))
        pseknclst +=elelst
    pseknclst.sort()
    psekncdict = {pseknclst[i]:i for i in range(len(pseknclst))}
    return psekncdict


def get_pseknc_vector(sequence, sstructure, k, letter_list = ['A', 'C', 'G', 'U', 'A-U', 'U-A', 'G-C', 'C-G', 'G-U', 'U-G']):
    '''This is a process in PseKNC, aim to gernerate feature vector.
     
    :param sequence: An RNA sequence,string.
    :param sstructure: The corresponding secondary structure, string.
    :param psekncdict: All the features, dictionary.
    :param k: The length of K-tuple nucleotide composition.
    :param letter_list: default  ['A', 'C', 'G', 'U']
    :return: Feature vector through PseKNC.
    '''
    psekncdict = get_pseknc_dict(letter_list, k)
    vector=np.zeros((1, len(psekncdict)))
    correspseq = get_corresp_sequence(sequence,sstructure)
    pattern = zip(list(sequence),list(correspseq))
    
    for i in xrange(len(pattern)-k+1):
        stem=[]
        for x,y in pattern[i:i+k]:
            if x == '.' or y == '.':
                if x == '.':
                    stem.append(y)
                else:
                    stem.append(x)
            else:
                stem.append(x + '-' + y)
        stem_tuple= tuple(stem)
        position = psekncdict.get(stem_tuple)
        vector[0, position] += 1
        #print stem_tuple,position
    #return vector[0]
    return list(vector[0]/sum(vector[0]))
    #return [round(f,4) for f in list(vector[0]/sum(vector[0]))]
    
#=========================PseSSC===============================================    
def get_psessc_matrix(filename, k, r, w, pattern_list = ['A', 'C', 'G', 'U', 'A-U', 'U-A', 'G-C', 'C-G', 'G-U', 'U-G']):
    '''This is a complete process in PseSSC, aim to gernerate feature matrix.
     
       The FASTA format of the input file is as follows:    
       >Sequence name
       An RNA sequence should be consist of AGCU
       Second structure
 
    :param filename: Name of input file.
    :param k: The number of k adjacent structure statuses.
    :param r: The highest counted rank (or tier) of the structural correlation along a RNA chain.
    :param w: The wight of theta, from 0.1 to 1.
    :param pattern_list: Structure statuses, default:['A', 'C', 'G', 'U', 'A-U', 'U-A', 'G-C', 'C-G', 'G-U', 'U-G'].
    :return: Feature matrix through PseSSC.
    '''
    with open(filename) as f:
        seqsslst= get_rnasc_data(f)
    features = []
    for seqss in seqsslst:
        vector = get_psessc_vector(seqss.sequence,seqss.sstruc, k, r, w, pattern_list)
        features.append(vector)
    return features
    
    
def get_psessc_vector(sequence, sstructure, k, r, w, pattern_list = ['A', 'C', 'G', 'U', 'A-U', 'U-A', 'G-C', 'C-G', 'G-U', 'U-G']):
    '''This is a complete process in PseSSC, aim to gernerate feature vector.
    :param sequence: an RNA sequence,string.
    :param sstructure: The corresponding secondary structure, string.
    :param k: The number of k adjacent structure statuses.
    :param r: The highest counted rank (or tier) of the structural correlation along a RNA chain.
    :param w: The wight of theta, from 0.1 to 1.
    :param pattern_list: Structure statuses, default:['A', 'C', 'G', 'U', 'A-U', 'U-A', 'G-C', 'C-G', 'G-U', 'U-G'].
    :return: Feature vector through PseSSC.'''
    #psekncdict = get_pseknc_dict(pattern_list, k)
    psekncvec = get_pseknc_vector(sequence,sstructure, k, pattern_list)
    
    psesscvec_tmp = np.array(psekncvec)
    for i in range(1, r+1):
        psesscvec_tmp = np.hstack((psesscvec_tmp, w * calculate_theta(sequence, sstructure, i)))
    psesscvec = psesscvec_tmp / sum(psesscvec_tmp)
    
    #return psesscvec
    return [round(f,4) for f in psesscvec]
       
def calculate_theta(sequence, sstructure, j):
    '''calculate theta
    :param sequence: an RNA sequence,string.
    :param sstructure: The corresponding secondary structure, string.
    :param j: the counted rank (or tier) of the structural correlation along a RNA chain.
    :return: theta.
    '''
    if j >= len(sequence):
        error_info = 'r should be less than the length of the sequence.'
        sys.stderr.write(error_info)
    else:
        correspseq = get_corresp_sequence(sequence,sstructure)
        pattern = zip(list(sequence),list(correspseq))
        stem=[]
        for x,y in pattern:
            if x == '.' or y == '.':
                if x == '.':
                    stem.append(y)
                else:
                    stem.append(x)
            else:
                stem.append(x + '-' + y)
        freevalue_vector = []
        for i in stem:
            if i == 'A-U' or i == 'U-A':
               freevalue_vector.append(-2)
            elif i == 'C-G' or i =='G-C':
                freevalue_vector.append(-3)
            elif i == 'U-G' or i =='G-U':
                freevalue_vector.append(-1)
            else:
                freevalue_vector.append(0)
                
        s=0.0
        for i in range(len(freevalue_vector)-j):
            s += (freevalue_vector[i] - freevalue_vector[i+j]) ** 2
            #print i,i+j
            #print (freevalue_vector[i] - freevalue_vector[i+r])
        #print s,len(freevalue_vector)-r
        return s / (len(freevalue_vector)-j)
  
#=========================PseDPC===============================================
    
def get_psedpc_vector(sequence, sstructure, d, r, w, pattern_list = ['A', 'C', 'G', 'U', 'A-U', 'U-A', 'G-C', 'C-G', 'G-U', 'U-G']):
    '''This is a complete process in PseSSC, aim to gernerate feature vector.
    :param sequence: An RNA sequence,string.
    :param sstructure: The corresponding secondary structure, string.
    :param d: Tthe distance between structure statuses.
    :param r: The highest counted rank (or tier) of the structural correlation along a RNA chain.
    :param w: The wight of theta, from 0.1 to 1.
    :param pattern_list: Structure statuses, default:['A', 'C', 'G', 'U', 'A-U', 'U-A', 'G-C', 'C-G', 'G-U', 'U-G'].'''
    correspseq = get_corresp_sequence(sequence, sstructure)
    pattern = zip(list(sequence), list(correspseq))
    vector = []
    vector = np.array(vector, ndmin = 2)
    for i in range(d+1):
        if i != 0:
            k = 2
        else:
            k = 1
        psedpcdict = get_pseknc_dict(pattern_list, k)
        vec_tmp = np.zeros((1, len(psedpcdict)))
        
        
        for j in xrange(len(pattern)-i):
            stem=[]
            if i != 0:
                for x, y in pattern[j], pattern[j+i]:
                    
                    if x == '.' or y == '.':
                        if x == '.':
                            stem.append(y)
                        else:
                            stem.append(x)
                    else:
                        stem.append(x + '-' + y)
                stem_tuple = tuple(stem)
                #print stem_tuple
            else:
                for x,y in [pattern[j]]:
                    if x == '.' or y == '.':
                        if x == '.':
                            stem.append(y)
                        else:
                            stem.append(x)
                    else:
                        stem.append(x + '-' + y)
                stem_tuple = tuple(stem)
                #print stem_tuple
            position = psedpcdict.get(stem_tuple)
            vec_tmp[0, position] += 1                
        vector = np.hstack((vector, vec_tmp))
        
    psedpcvec_tmp = vector[0]
    for i in range(1, r+1):
        psedpcvec_tmp = np.hstack((psedpcvec_tmp,(w * calculate_theta(sequence, sstructure, i))))
    psedpcvec = psedpcvec_tmp / (sum(psedpcvec_tmp[10+100*d:])+1)
   # psedpcvec = psesscvec_tmp / sum(psesscvec_tmp)
    #return psedpcvec
    return [round(f,4) for f in psedpcvec]
    
    
def get_psedpc_matrix(filename, d, r, w, pattern_list = ['A', 'C', 'G', 'U', 'A-U', 'U-A', 'G-C', 'C-G', 'G-U', 'U-G']):
    '''This is a complete process in PseSSC, aim to gernerate feature vector.
     
       The FASTA format of the input file is as follows:    
       >sequence name
       An RNA sequence should be consist of AGCU
       Second structure
    :param filename: Name of input file.
    :param d: The distance between structure statuses.
    :param r: The highest counted rank (or tier) of the structural correlation along a RNA chain.
    :param w: The wight of theta, from 0.1 to 1.
    :param pattern_list: Structure statuses, default:['A', 'C', 'G', 'U', 'A-U', 'U-A', 'G-C', 'C-G', 'G-U', 'U-G'].'''
    with open(filename) as f:
        seqsslst= get_rnasc_data(f)
    features = []
    for seqss in seqsslst:
        vector = get_psedpc_vector(seqss.sequence, seqss.sstruc, d, r, w, pattern_list)
        features.append(vector)
    return features
    
    
def main(args):
    #TODO:args.method will be finished
    #TODO:args.inputfile, name 
    
    
    
    
    if args.method == 'Triplet':
        res = get_triplet_matrix(args.inputfile)
    elif args.method == 'PseSSC':
        if args.k is None:
            print "parameters k is required. The default value of k is 2."
            args.k = 2
        if args.r is None:
            print "parameters r is required. The default value of r is 2."
            args.r = 2
        if args.w is None:
            print "parameters w is required. The default value of w is 0.1."
            args.w = 0.1
        res = get_psessc_matrix(args.inputfile, args.k, args.r, args.w)
    elif args.method == 'PseDPC':
        if args.d is None:
            print "parameters d is required. The default value of d is 0."
            args.d = 0
        if args.r is None:
            print "parameters r is required. The default value of r is 2."
            args.r = 2
        if args.w is None:
            print "parameters w is required. The default value of w is 0.1."
            args.w = 0.1
        res = get_psedpc_matrix(args.inputfile, args.d, args.r, args.w)
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
     import argparse
     from argparse import RawTextHelpFormatter
 
     parse = argparse.ArgumentParser(description="This is a structure composition module for generate feature vectors.",
                                     formatter_class=RawTextHelpFormatter)
     parse.add_argument('inputfile',
                        help="The input file, in valid FASTA format.")
     parse.add_argument('outputfile',
                        help="The outputfile stored results.")
     parse.add_argument('alphabet', choices = ['DNA', 'RNA', 'Protein'],
                        help="The alphabet of sequences.")
     parse.add_argument('method', type=str,
                        help="The method name of structure composition.")
     parse.add_argument('-k', type=int, 
                        help="The number of k adjacent structure statuses. default=2")
     parse.add_argument('-d', type=int,
                        help="The distance between structure statuses. default=0")
     parse.add_argument('-r', type=int,
                        help="The value of lambda, represents the highest counted rank (or tier) of the structural correlation along a RNA chain. default=2")
     parse.add_argument('-w', type=float,
                        help="weight factor used to adjust the effect of the correlation factors. default=0.1")
     
     #parse.add_argument('k', type=int, choices=range(1, 7),
      #                  help="The k value of kmer.")
    # parse.add_argument('alphabet', choices=['DNA', 'RNA', 'PROTEIN'],
     #                   help="The alphabet of sequences.")
     #parse.add_argument('-r', default=0, type=int, choices=[1, 0],
      #                  help="Whether need to reverse complement.\n"
       #                      "1 means True, 0 means False. (default = 0)")
     parse.add_argument('-f', default='tab', choices=['tab', 'svm', 'csv'],
                        help="The output format (default = tab).\n"
                             "tab -- Simple format, delimited by TAB.\n"
                             "svm -- The libSVM training data format.\n"
                             "csv -- The format that can be loaded into a spreadsheet program.")
     parse.add_argument('-l', default='+1', choices=['+1', '-1'],
                        help="The libSVM output file label.")
 
     args = parse.parse_args()
 
     print(args)
     print("Calculating...")
     start_time = time.time()
     main(args)
     
     print("Used time: %ss" % (time.time() - start_time))
     print("Done.")
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
    
#==========================PseSSC test==========================================
#     pattern_list = ['A', 'C', 'G', 'U', 'A-U', 'U-A', 'G-C', 'C-G', 'G-U', 'U-G']
#     letter_list= ['A', 'C', 'G', 'U']
#     sequence = "GCAUCCGGGUUGAGGUAGUAGGUUGUAUGGUUUAGAGUUACACCCUGGGAGUUAACUGUACAACCUUCUAGCUUUCCUUGGAGC"
#     sstructure = '((.((((((..(((.(((.(((((((((((((..((.(..((...))..).))))))))))))))).))).)))..))))))))'
# 
#     k = 1
#     r = 1
#     w = 1
#     psesscvec = get_psessc_vector(sequence, sstructure, k, r, w)
#     psesscvec_f = get_psessc_matrix("test.txt", k, r, w, pattern_list)
#     
#     psekncdict = get_pseknc_dict(letter_list,3)
#     psekncvec = get_pseknc_vector(sequence, sstructure,1)
#==============================================================================
#    pattern_list = ['A', 'C', 'G', 'U', 'A-U', 'U-A', 'G-C', 'C-G', 'G-U', 'U-G']
#    sequence = "GCAUCCGGGUUGAGGUAGUAGGUUGUAUGGUUUAGAGUUACACCCUGGGAGUUAACUGUACAACCUUCUAGCUUUCCUUGGAGC"
#    sstructure = '((.((((((..(((.(((.(((((((((((((..((.(..((...))..).))))))))))))))).))).)))..))))))))'
#    d=1
#    r=1
#    w=1
#    psedpc_vec=get_psedpc_vector(sequence,sstructure,d,r,w,pattern_list)
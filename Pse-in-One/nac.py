"""
Created on Sat May 13 15:35:42 2016
@version:0.2.1./pyc
@author: Fule Liu, Nackel, luo
"""

import sys
import re
import time
import multiprocessing

from numpy import array
from itertools import combinations, combinations_with_replacement, permutations
import numpy as np

from util import frequency
from util import get_data
from data import index_list


#===========================Kmer===================================================
def make_kmer_list(k, alphabet):
    if k < 0:
        print("Error, k must be an inter and larger than 0.")
 
    kmers = []
    for i in range(1, k + 1):
        if len(kmers) == 0:
            kmers = list(alphabet)
        else:
            new_kmers = []
            for kmer in kmers:
                for c in alphabet:
                    new_kmers.append(kmer + c)
            kmers = new_kmers
 
    return kmers
 
 
def find_revcomp(sequence, revcomp_dictionary):
    # Save time by storing reverse complements in a hash.
    if sequence in revcomp_dictionary:
        return revcomp_dictionary[sequence]
 
    # Make a reversed version of the string.
    rev_sequence = list(sequence)
    rev_sequence.reverse()
    rev_sequence = ''.join(rev_sequence)

    return_value = ""
    for letter in rev_sequence:
        if letter == "A":
            return_value += "T"
        elif letter == "C":
            return_value += "G"
        elif letter == "G":
            return_value += "C"
        elif letter == "T":
            return_value += "A"
        elif letter == "N":
            return_value += "N"
        else:
            error_info = ("Unknown DNA character (%s)\n" % letter)
            sys.exit(error_info)
 
    # Store this value for future use.
    revcomp_dictionary[sequence] = return_value
 
    return return_value
 
 
def _cmp(a, b):
    return (a > b) - (a < b)
 
 
def make_revcomp_kmer_list(kmer_list):
    revcomp_dictionary = {}
    new_kmer_list = [kmer for kmer in kmer_list if _cmp(kmer, find_revcomp(kmer, revcomp_dictionary)) <= 0]
    return new_kmer_list
 
 
def make_kmer_vector(k, alphabet, filename, revcomp=False):
    """Generate kmer vector."""
    with open(filename) as f:
        seq_list = get_data(f, alphabet=alphabet)

        if revcomp and re.search(r'[^acgtACGT]', ''.join(alphabet)) is not None:
            sys.exit("Error, Only DNA sequence can be reverse compliment.")
 
        vector = []
        kmer_list = make_kmer_list(k, alphabet)
        for seq in seq_list:
            count_sum = 0
 
            # Generate the kmer frequency dict.
            kmer_count = {}
            for kmer in kmer_list:
                temp_count = frequency(seq, kmer)
                if not revcomp:
                    if kmer not in kmer_count:
                        kmer_count[kmer] = 0
                    kmer_count[kmer] += temp_count
                else:
                    rev_kmer = find_revcomp(kmer, {})
                    if kmer <= rev_kmer:
                        if kmer not in kmer_count:
                            kmer_count[kmer] = 0
                        kmer_count[kmer] += temp_count
                    else:
                        if rev_kmer not in kmer_count:
                            kmer_count[rev_kmer] = 0
                        kmer_count[rev_kmer] += temp_count
 
                count_sum += temp_count
 
            # Normalize.
            if not revcomp:
                count_vec = [kmer_count[kmer] for kmer in kmer_list]
            else:
                revc_kmer_list = make_revcomp_kmer_list(kmer_list)
                count_vec = [kmer_count[kmer] for kmer in revc_kmer_list]
            count_vec = [round(float(e)/count_sum, 8) for e in count_vec]

            vector.append(count_vec)

    return vector
#==============================================================================
    
    
#==================getting (k,m)-mismatch profile==============================
def getMismatchProfileMatrix(filename, alphabet, k, m):
    alphabet = list(alphabet)
    p = len(alphabet)
    with open(filename) as f:
        seq_list = get_data(f, alphabet=alphabet)
        kmerdict = getKmerDict(alphabet, k)
        features = []
        if m==0 and m < k:
            for sequence in seq_list:
                vector=getSpectrumProfileVector(sequence, kmerdict, p, k)
                features.append(vector)
        elif m > 0 and m < k:
            for sequence in seq_list:
                vector=getMismatchProfileVector(sequence, alphabet, kmerdict, p, k)
                features.append(vector)   
    return array(features)
        
def getKmerDict(alphabet, k):
    kmerlst = []
    partkmers = list(combinations_with_replacement(alphabet, k))
    for element in partkmers:
        elelst = set(permutations(element, k))
        strlst = [''.join(ele) for ele in elelst]
        kmerlst += strlst
    kmerlst = np.sort(kmerlst)
    kmerdict = {kmerlst[i]:i for i in range(len(kmerlst))}
    return kmerdict


def getSpectrumProfileVector(sequence, kmerdict, p, k):    
    vector = np.zeros((1, p**k))
    n = len(sequence)
    for i in range(n-k+1):
        subsequence=sequence[i:i+k]
        position=kmerdict.get(subsequence)
        vector[0,position] += 1
    return list(vector[0])


def getMismatchProfileVector(sequence, alphabet, kmerdict, p, k): 
    n = len(sequence)
    vector = np.zeros((1, p**k))
    for i in range(n-k+1):
        subsequence = sequence[i:i+k]
        position = kmerdict.get(subsequence)
        vector[0, position]+=1
        for j in range(k):
            substitution = subsequence
            for letter in list(set(alphabet)^set(subsequence[j])):
                substitution = list(substitution)
                substitution[j] = letter
                substitution = ''.join(substitution)
                position = kmerdict.get(substitution)
                vector[0,position] += 1    
    return list(vector[0])

#==============================================================================

#=================getting (k, delta)-subsequence profile=======================
def getSubsequenceProfileByParallel(filename, alphabet, k, delta):
    alphabet = list(alphabet)
    with open(filename) as f:
        seq_list = get_data(f, alphabet=alphabet)
        cpu_num = multiprocessing.cpu_count()   
        batches = constructPartitions(seq_list, cpu_num)
        pool = multiprocessing.Pool(cpu_num)
        results = []
        
        for batch in batches:
            temp=pool.apply_async(getSubsequenceProfile, (batch, alphabet, k, delta))
            results.append(temp)
        pool.close()
        pool.join()
        i = 1
        for temp in results:
            temp_X = temp.get()
            if len(temp_X) != 0:
                if i == 1:
                    X = temp_X
                else:
                    X = np.vstack((X,temp_X))
                i += 1
        return X
    
def constructPartitions(seq_list, cpu_num):
    seqs_num = len(seq_list)
    batch_num = seqs_num//cpu_num
    batches = []
    for i in range(cpu_num-1):
        batch = seq_list[i*batch_num:(i+1)*batch_num]
        batches.append(batch)
    batch=seq_list[(cpu_num-1)*batch_num:]
    batches.append(batch)
    return batches
    
def getSubsequenceProfile(seq_list, alphabet, k, delta):
    kmerdict = getKmerDict(alphabet, k)
    X=[]
    for sequence in seq_list:
        vector = getSubsequenceProfileVector(sequence, kmerdict, k, delta)
        X.append(vector)
    X=array(X)    
    return X

def getSubsequenceProfileVector(sequence, kmerdict, k, delta):      
    vector = np.zeros((1,len(kmerdict)))
    sequence = array(list(sequence))
    n = len(sequence)
    index_lst = list(combinations(range(n), k))
    for subseq_index in index_lst:
        subseq_index = list(subseq_index)
        subsequence = sequence[subseq_index]
        position = kmerdict.get(''.join(subsequence))     
        subseq_length = subseq_index[-1] - subseq_index[0] + 1
        subseq_score = 1 if subseq_length == k else delta**subseq_length    
        vector[0,position] += subseq_score
    #return list(vector[0])
    return [round(f, 4) for f in list(vector[0])]
#==============================================================================


def main(args):
    # Set revcomp parameter.
    if args.r != 1:
        args.r = False
    elif args.r == 1 and args.alphabet != 'DNA':
        print("Error, the -r parameter can only be used in DNA.")
    elif args.r == 1 and args.alphabet == 'DNA':
        args.r = True

    # Set alphabet parameter.
    if args.alphabet == 'DNA':
        args.alphabet = index_list.DNA
    elif args.alphabet == 'RNA':
        args.alphabet = index_list.RNA
    elif args.alphabet == 'Protein':
        args.alphabet = index_list.PROTEIN
        
    if args.method.upper() == 'KMER': 
        if args.k is None:
            print "parameters k is required. The default value of k is 2."
            args.k = 2
        if args.r is None:
            print "parameters r is required. The default value of r is 0."
            args.r = 0
        res = make_kmer_vector(k=args.k, alphabet=args.alphabet, filename=args.inputfile, revcomp=args.r)
    elif args.method.upper() == "MISMATCH":
        if args.k is None:
            print "parameters k is required. The default value of k is 3."
            args.k = 3
        if args.m is None:
            print "parameters m is required. The default value of m is 1."
            args.m = 1
        if args.m >= args.k:
            print "parameters m should be less than parameter k."
        else:
            res = getMismatchProfileMatrix(args.inputfile, args.alphabet, args.k, args.m)
    elif args.method.upper() == "SUBSEQUENCE":
        if args.delta is None:
            print "parameters delta is required. The default value of delta is 1."
            args.delta = 1
        elif args.delta > 1 or args.delta < 0:
            print "delta should be greater than or equal to 0 and less than or equal to 1."
        
        if args.k is None:
            print "parameters k is required. The default value of k is 3."
            args.k = 3 
        res = getSubsequenceProfileByParallel(filename=args.inputfile, alphabet=args.alphabet, k=args.k, delta=args.delta)
            
    else:
            print("Method error!")

    # Write correspond res file.
    if args.f == 'svm':
        from util import write_libsvm
        write_libsvm(res, [args.l] * len(res), args.outputfile)
    elif args.f == 'tab':
        from util import write_tab
        write_tab(res, args.outputfile)
    elif args.f == 'csv':
        from util import write_csv
        write_csv(res, args.outputfile)


if __name__ == '__main__':
    import argparse
    from argparse import RawTextHelpFormatter

    parse = argparse.ArgumentParser(description="This is kmer module for generate nucleic acid compositio vector.",
                                    formatter_class=RawTextHelpFormatter)
    parse.add_argument('inputfile',
                       help="The input file in FASTA format.")
    parse.add_argument('outputfile',
                       help="The output file stored results.")
    parse.add_argument('alphabet', choices=['DNA', 'RNA', 'Protein'],
                       help="The sequence type.")
    parse.add_argument('method', type=str,
                        help="The method name of nucleic acid composition. {Kmer,mismatch,subsequence}")
    parse.add_argument('-k', type=int, 
                       help="For Kmer, mismatch, subsequence methods. The k value of kmer.")
    parse.add_argument('-m', type=int, default=1,
                       help="For mismatch. The max value inexact matching. (m<k)")
    parse.add_argument('-delta', type=float, default=1,
                       help="For subsequence method. The value of penalized factor. (0<=delta<=1)")
    
    parse.add_argument('-r', type=int, choices=[1, 0],
                       help="Whether consider the reverse complement or not.\n"
                            "1 means True, 0 means False. (default = 0)")
    parse.add_argument('-f', default='tab', choices=['tab', 'svm', 'csv'],
                       help="The output format (default = tab).\n"
                            "tab -- Simple format, delimited by TAB.\n"
                            "svm -- The libSVM training data format.\n"
                            "csv -- The format that can be loaded into a spreadsheet program.")
    parse.add_argument('-l', default='+1', choices=['+1', '-1'],
                       help="The libSVM output file label.")

    args = parse.parse_args()

    print("Calculating...")
    start_time = time.time()
    main(args)
    print("Done.")
    print("Used time: %ss" % (time.time() - start_time))
    
    
    
    
#===========================test misnatch==========================================
    #matrix = getMismatchProfileMatrix("test_s.txt", index_list.RNA, 3,0)
#    f=open("test_s.txt")
#    alphabet = index_list.DNA
#    k = 3
#    delta = 1
#    matrix_subseq=getSubsequenceProfileByParallel("test_s.txt", alphabet, k, delta)
#==============================================================================
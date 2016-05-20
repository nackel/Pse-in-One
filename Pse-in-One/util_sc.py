# -*- coding: utf-8 -*-
"""
Created on Fri May 20 10:05:20 2016
@version:0.0.1
@author: Nackel
"""

import sys
import re
from data import index_list
from util import is_under_alphabet

class SeqSS:
    """An object that contains RNA sequence and its secondary structure"""
    def __init__(self, name, sequence, sstruc, no, mfe):
        self.name = name
        self.sequence = sequence.upper()
        self.sstruc = sstruc            #its secondary structure
        self.no = no
        self.length = len(sequence) if len(sequence)==len(sstruc) else -1
        self.MFE = mfe          #float

    def __str__(self):
        """Output seqSS when 'print' method is called."""
        return "%s\tNo:%s\tlength:%s\tMFE:%.2f\n%s\n%s" % (self.name, str(self.no), str(self.length), self.MFE, self.sequence, self.sstruc)
        
        
def get_rnasc_data(input_data):
    """Get rnasc data from file or list with check.

    :param input_data: type file or list
    :param desc: with this option, the return value will be a SeqSS object list(it only works in file object).
    :return: a SeqSS object or shutdown.
    """
    if hasattr(input_data, 'read'):
        return read_rnasc_fasta(input_data)
    elif isinstance(input_data, list):
        input_data = is_rnasc_list(input_data)
        if input_data is not False:
            return input_data
        else:
            sys.exit(0)
    else:
        error_info = 'Sorry, the parameter in get_rnasc_data method must be list or file type.'
        sys.exit(error_info)
        
        
def is_rnasc_list(rnasc_list): 
    """Judge the rnasc_list is within the scope of alphabet and change the lowercase to capital.
    The format of the input list is as follows:
        [">sequence name","RNA sequence should be consist of AGCU","second structure"]"""
    seqss_list = []
    count = 0    # record line numbers
    no = 0       # record sequence numbers
    alphabet = index_list.RNA
    name,sequence,sstruc = '','',''
    seqss_name = r'>(.+)'
    seqss_sequence = r'([ACUGacug]+)'
    seqss_sstructure = r'([\(\)\.]+)'
    seqss_MFE = r'-\d+\.\d+'
    mfe = ""
    for line in rnasc_list:
            
        if not line:
            break
            
        if count %3 == 0:
            name = re.findall(seqss_name, line)[0]
        elif count%3 == 1:
            res = is_under_alphabet(line, alphabet)
            if res is not True:
                error_info = 'Error, sequence ' + str(no) \
                             + ' has character ' + str(res) + '.(The character must be ' + alphabet + ').'
                sys.stderr.write(error_info)
                return False
            else:
                sequence = re.findall(seqss_sequence, line)[0]
                
        elif count%3 == 2:
            sstruc = re.findall(seqss_sstructure, line)[0]
            if re.findall(seqss_MFE, line):
                mfe = float(re.findall(seqss_MFE, line)[0])
            else:
                mfe = 0                        #MFE is no match.
            if sequence != "" and sstruc != "":
                if is_rnasc_fasta(SeqSS(name, sequence, sstruc, no, mfe)):
                    seqss_list.append(SeqSS(name,sequence,sstruc,no, mfe))
                    no += 1
                else:
                    sys.exit()    
        count += 1
    return seqss_list
    
    
def read_rnasc_fasta(f):
    """Read a fasta file.
    
    :param f: HANDLE to input. e.g. sys.stdin, or open(<file>)
    
    Return SeqSS obj list.
    """
    seqss_list = []
    count = 0    # record line numbers
    no = 0       # record sequence numbers
    name,sequence,sstruc = '','',''
    seqss_name = r'>(.+)'
    seqss_sequence = r'([ACUG]+)'
    seqss_sstructure = r'([\(\)\.]+)'
    seqss_MFE = r'-\d+\.\d+'
    mfe = ""
    for line in f:
            
        if not line:
            break
            
        if count %3 == 0:
            name = re.findall(seqss_name, line)[0]
        elif count%3 == 1:
            sequence = re.findall(seqss_sequence, line)[0]
        elif count%3 == 2:
            sstruc = re.findall(seqss_sstructure, line)[0]
            if re.findall(seqss_MFE, line):
                mfe = float(re.findall(seqss_MFE, line)[0])
            else:
                mfe = 0  #MFE is no match.
            if is_rnasc_fasta(SeqSS(name, sequence, sstruc, no, mfe)):
                seqss_list.append(SeqSS(name,sequence,sstruc,no, mfe))
                no += 1
            else:
                sys.exit()
        count += 1
    return seqss_list
    
   
def is_rnasc_fasta(seqss):
    """Judge the SeqSS object is in FASTA format.
    Five situations will:
    1. No seqss name.
    2. SeqSS name is illegal.
    3. No sequence.
    4. The length of sequence is not equal to the length of its secondary structure
    5. The RNA sequence was not composed by AUCG
    :param seqss: SeqSS object.
    """
    alphabet ="AUCG"
    if not seqss.name:
        error_info = 'Error, sequence ' + str(seqss.no) + ' has no sequence name.'
        sys.stderr.write(error_info)
        return False
    if 0 == seqss.length:
        error_info = 'Error, sequence ' + str(seqss.no) + ' is null.'
        sys.stderr.write(error_info)
        return False
    if len(seqss.sequence) != len(seqss.sstruc):
        error_info = 'Error, the length of sequence ' + str(seqss.no) + ' is not equal to the length of its secondary structure.'
        sys.stderr.write(error_info)
        return False
        
    res = is_under_alphabet(seqss.sequence, alphabet)
    if res is not True:
        error_info = 'Error, sequence ' + str(seqss.no) \
                         + ' has character ' + str(res) + '.(The character must be ' + alphabet + ').'
        sys.stderr.write(error_info)
        return False
                         
    return True


def get_rnasc_sequences(f):
    """Read the fasta file.

    Input: f: HANDLE to input. e.g. sys.stdin, or open(<file>)

    Return a list of sequences.
    """
    seqss_sequences_list = []
    seqss_list = read_rnasc_fasta(f)
    for seqss in seqss_list:
        seqss_sequences_list.append(seqss.sequence)
    return seqss_sequences_list
    
    
def get_rnasc_sstructures(f):
    """Read the fasta file.

    Input: f: HANDLE to input. e.g. sys.stdin, or open(<file>)

    Return a secondary structure list.
    """
    seqss_sstrucs_list = []
    seqss_list = read_rnasc_fasta(f)
    for seqss in seqss_list:
        seqss_sstrucs_list.append(seqss.sstruc)
    return seqss_sstrucs_list
    


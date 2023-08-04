#!/usr/bin/env python
# _*_coding:utf-8_*_

import re, os, sys

def read_nucleotide_sequences(file):
    if os.path.exists(file) == False:
        print('Error: file %s does not exist.' % file)
        sys.exit(1)
    with open(file) as f:
        records = f.read()
    if re.search('>', records) == None:
        print('Error: the input file %s seems not in FASTA format!' % file)
        sys.exit(1)
    records = records.split('>')[1:]
    fasta_sequences = []
    for fasta in records:
        array = fasta.split('\n')
        header, sequence = array[0].split()[0], re.sub('[^ACGTU-]', '-', ''.join(array[1:]).upper())
        header_array = header.split('|')
        name = header_array[0]
        label = header_array[1] if len(header_array) >= 2 else '0'
        label_train = header_array[2] if len(header_array) >= 3 else 'training'
        sequence = re.sub('U', 'T', sequence)
        fasta_sequences.append([name, sequence, label, label_train])
    return fasta_sequences

def read_protein_sequences(file):
    if os.path.exists(file) == False:
        print('Error: file %s does not exist.' % file)
        sys.exit(1)
    with open(file) as f:
        records = f.read()
    if re.search('>', records) == None:
        print('Error: the input file %s seems not in FASTA format!' % file)
        sys.exit(1)
    records = records.split('>')[1:]
    fasta_sequences = []
    for fasta in records:
        array = fasta.split('\n')
        header, sequence = array[0].split()[0], re.sub('[^ACDEFGHIKLMNPQRSTVWY-]', '-', ''.join(array[1:]).upper())
        header_array = header.split('|')
        name = header_array[0]
        label = header_array[1] if len(header_array) >= 1 else '0'
        label_train = header_array[2] if len(header_array) >= 2 else 'training'
        fasta_sequences.append([name, sequence, label, label_train])
    return fasta_sequences

def readFasta(file):
    if os.path.exists(file) == False:
        print('Error: "' + file + '" does not exist.')
        sys.exit(1)

    with open(file) as f:
        records = f.read()

    if re.search('>', records) == None:
        print('The input file seems not in fasta format.')
        sys.exit(1)

    records = records.split('>')[1:]
    myFasta = []
    for fasta in records:
        array = fasta.split('\n')
        name, sequence = array[0].split()[0], re.sub('[^ARNDCQEGHILKMFPSTWYV-]', '-', ''.join(array[1:]).upper())
        myFasta.append([name, sequence])
    return myFasta
#!/usr/bin/env python
#_*_coding:utf-8_*_

#!/usr/bin/env python
#_*_coding:utf-8_*_

import re

def check_fasta_with_equal_length(fastas):
    status = True
    lenList = set()
    for i in fastas:
        lenList.add(len(i[1]))
    if len(lenList) == 1:
        return True
    else:
        return False

def get_min_sequence_length(fastas):
    minLen = 10000
    for i in fastas:
        if minLen > len(i[1]):
            minLen = len(i[1])
    return minLen

def get_min_sequence_length_1(fastas):
    minLen = 10000
    for i in fastas:
        if minLen > len(re.sub('-', '', i[1])):
            minLen = len(re.sub('-', '', i[1]))
    return minLen

import re

def check_fasta_with_equal_length(fastas):
    status = True
    lenList = set()
    for i in fastas:
        lenList.add(len(i[1]))
    if len(lenList) == 1:
        return True
    else:
        return False

def get_min_sequence_length(fastas):
    minLen = 10000
    for i in fastas:
        if minLen > len(i[1]):
            minLen = len(i[1])
    return minLen

def get_min_sequence_length_1(fastas):
    minLen = 10000
    for i in fastas:
        if minLen > len(re.sub('-', '', i[1])):
            minLen = len(re.sub('-', '', i[1]))
    return minLen
#!/usr/bin/env python
#_*_coding:utf-8_*_

import re

def checkFasta(fastas):
	status = True
	lenList = set()
	for i in fastas:
		lenList.add(len(i[1]))
	if len(lenList) == 1:
		return True
	else:
		return False

def minSequenceLength(fastas):
	minLen = 10000
	for i in fastas:
		if minLen > len(i[1]):
			minLen = len(i[1])
	return minLen

def minSequenceLengthWithNormalAA(fastas):
	minLen = 10000
	for i in fastas:
		if minLen > len(re.sub('-', '', i[1])):
			minLen = len(re.sub('-', '', i[1]))
	return minLen

import pandas as pd
import numpy as np
import tensorflow as tf
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import KFold
import math
#import read_fasta_sequences
#import check_sequences
import itertools
from collections import Counter
from keras.layers import Input, Dense, Conv1D, Flatten, AveragePooling1D,MaxPooling1D,  Dropout, \
    Reshape, normalization,LeakyReLU
from tensorflow.keras.layers import BatchNormalization
#from keras.layers.embeddings import Embedding
from sklearn.preprocessing import scale
from keras.layers import GRU,Bidirectional,LSTM
from keras.models import Model
import keras.backend as K
import h5py
from sklearn import metrics
import matplotlib.pyplot as plt
import random
from keras.regularizers import l1, l2
from tensorflow.python.keras.optimizers import TFOptimizer
#from DNARNAFeatureExtraction import*

def ANF(fastas, **kw):
    if check_sequences.check_fasta_with_equal_length == False:
        print('Error: for "ANF" encoding, the input fasta sequences should be with equal length. \n\n')
        return 0

    AA = 'ACGT'
    encodings = []
    #header = ['#','label']
    #for i in range(1, len(fastas[0][1]) + 1):
    #    header.append('ANF.' + str(i))
    #encodings.append(header)

    for i in fastas:
        name, sequence  = i[0], i[1]
        code = []
        for j in range(len(sequence)):
            code.append(sequence[0: j + 1].count(sequence[j]) / (j + 1))
        encodings.append(code)
    return np.array(encodings,dtype=float)
def binary(fastas, **kw):
    if check_sequences.check_fasta_with_equal_length == False:
        print('Error: for "BINARY" encoding, the input fasta sequences should be with equal length. \n\n')
        return 0

    AA = 'ACGT'
    encodings = []
    #header = ['#']
    #for i in range(1, len(fastas[0][1]) * 4 + 1):
    #    header.append('BINARY.F'+str(i))
    #encodings.append(header)

    for i in fastas:
        name, sequence = i[0], i[1]
        code = []
        for aa in sequence:
            if aa == '-':
                code = code + [0, 0, 0, 0]
                continue
            for aa1 in AA:
                tag = 1 if aa == aa1 else 0
                code.append(tag)
        encodings.append(code)
    return np.array(encodings,dtype=float)

def CKSNAP(fastas, gap, **kw):
    if gap < 0:
        print('Error: the gap should be equal or greater than zero' + '\n\n')
        return 0

    if check_sequences.get_min_sequence_length(fastas) < gap + 2:
        print('Error: all the sequence length should be larger than the (gap value) + 2 = ' + str(gap + 2) + '\n\n')
        return 0

    AA = kw['order'] if kw['order'] != None else 'ACGT'
    encodings = []
    aaPairs = []
    for aa1 in AA:
        for aa2 in AA:
            aaPairs.append(aa1 + aa2)

    #header = ['#']
    #for g in range(gap + 1):
    #    for aa in aaPairs:
     #       header.append(aa + '.gap' + str(g))
    #encodings.append(header)

    for i in fastas:
        name, sequence = i[0], i[1]
        code = []
        for g in range(gap + 1):
            myDict = {}
            for pair in aaPairs:
                myDict[pair] = 0
            sum = 0
            for index1 in range(len(sequence)):
                index2 = index1 + g + 1
                if index1 < len(sequence) and index2 < len(sequence) and sequence[index1] in AA and sequence[
                    index2] in AA:
                    myDict[sequence[index1] + sequence[index2]] = myDict[sequence[index1] + sequence[index2]] + 1
                    sum = sum + 1
            for pair in aaPairs:
                code.append(myDict[pair] / sum)
        encodings.append(code)
    return np.array(encodings,dtype=float)
def DNC(fastas, **kw):
    base = 'ACGT'

    encodings = []
    dinucleotides = [n1 + n2 for n1 in base for n2 in base]
    #header = ['#', 'label'] + dinucleotides
    #encodings.append(header)

    AADict = {}
    for i in range(len(base)):
        AADict[base[i]] = i

    for i in fastas:
        name, sequence = i[0], re.sub('-', '', i[1])
        code = []
        tmpCode = [0] * 16
        for j in range(len(sequence) - 2 + 1):
            tmpCode[AADict[sequence[j]] * 4 + AADict[sequence[j+1]]] = tmpCode[AADict[sequence[j]] * 4 + AADict[sequence[j+1]]] +1
        if sum(tmpCode) != 0:
            tmpCode = [i/sum(tmpCode) for i in tmpCode]
        code = code + tmpCode
        encodings.append(code)
    return np.array(encodings,dtype=float)
def EIIP(fastas, **kw):
    if check_sequences.check_fasta_with_equal_length == False:
        print('Error: for "EIIP" encoding, the input fasta sequences should be with equal length. \n\n')
        return 0

    AA = 'ACGT'

    EIIP_dict ={
        'A': 0.1260,
        'C': 0.1340,
        'G': 0.0806,
        'T': 0.1335,
        '-': 0,
    }

    encodings = []
    #header = ['#', 'label']
    #for i in range(1, len(fastas[0][1]) + 1):
    #    header.append('F'+str(i))
    #encodings.append(header)

    for i in fastas:
        name, sequence = i[0], i[1]
        code = []
        for aa in sequence:
            code.append(EIIP_dict.get(aa, 0))
        encodings.append(code)
    return np.array(encodings,dtype=float)
def ENAC(fastas, window=5, **kw):
    if check_sequences.check_fasta_with_equal_length == False:
        print('Error: for "ENAC" encoding, the input fasta sequences should be with equal length. \n\n')
        return 0

    if window < 1:
        print('Error: the sliding window should be greater than zero' + '\n\n')
        return 0

    # if check_sequences.get_min_sequence_length(fastas) < window:
        # print('Error: all the sequence length should be larger than the sliding window :' + str(window) + '\n\n')
        # return 0

    AA = kw['order'] if kw['order'] != None else 'ACGT'
    encodings = []
    #header = ['#']
    #for w in range(1, len(fastas[0][1]) - window + 2):
    #    for aa in AA:
     #       header.append('SW.' + str(w) + '.' + aa)
    #encodings.append(header)

    for i in fastas:
        name, sequence = i[0], i[1]
        code = []
        for j in range(len(sequence)):
            if j < len(sequence) and j + window <= len(sequence):
                count = Counter(sequence[j:j + window])
                for key in count:
                    count[key] = count[key] / len(sequence[j:j + window])
                for aa in AA:
                    code.append(count[aa])
        encodings.append(code)
    return np.array(encodings,dtype=float)
def kmerArray(sequence, k):
    kmer = []
    for i in range(len(sequence) - k + 1):
        kmer.append(sequence[i:i + k])
    return kmer


def Kmer(fastas, k=2, type="RNA", upto=False, normalize=True, **kw):
    encoding = []
    header = []
    NA = 'ACGT'
    if type in ("DNA", 'RNA'):
        NA = 'ACGT'
    else:
        NA = 'ACDEFGHIKLMNPQRSTVWY'

    if k < 1:
        print('Error: the k-mer value should larger than 0.')
        return 0

    if upto == True:
        for tmpK in range(1, k + 1):
            for kmer in itertools.product(NA, repeat=tmpK):
                header.append(''.join(kmer))
        encoding.append(header)
        for i in fastas:
            name, sequence= i[0], re.sub('-', '', i[1])
            count = Counter()
            for tmpK in range(1, k + 1):
                kmers = kmerArray(sequence, tmpK)
                count.update(kmers)
                if normalize == True:
                    for key in count:
                        if len(key) == tmpK:
                            count[key] = count[key] / len(kmers)
            code = []
            for j in range(2, len(header)):
                if header[j] in count:
                    code.append(count[header[j]])
                else:
                    code.append(0)
            encoding.append(code)
    else:
        for kmer in itertools.product(NA, repeat=k):
            header.append(''.join(kmer))
        #encoding.append(header)
        for i in fastas:
            name, sequence = i[0], re.sub('-', '', i[1])
            kmers = kmerArray(sequence, k)
            count = Counter()
            count.update(kmers)
            if normalize == True:
                for key in count:
                    count[key] = count[key] / len(kmers)
            code = []
            for j in range(2, len(header)):
                if header[j] in count:
                    code.append(count[header[j]])
                else:
                    code.append(0)
            encoding.append(code)
    return np.array(encoding,dtype=float)
chemical_property = {
    'A': [1, 1, 1],
    'C': [0, 1, 0],
    'G': [1, 0, 0],
    'T': [0, 0, 1],
    'U': [0, 0, 1],
    '-': [0, 0, 0],
}

def NCP(fastas, **kw):
    if check_sequences.check_fasta_with_equal_length == False:
        print('Error: for "NCP" encoding, the input fasta sequences should be with equal length. \n\n')
        return 0

    AA = 'ACGT'
    encodings = []
    #header = ['#', 'label']
    #for i in range(1, len(fastas[0][1]) * 3 + 1):
    #    header.append('NCP.F'+str(i))
    #encodings.append(header)

    for i in fastas:
        name, sequence = i[0], i[1]
        code = []
        for aa in sequence:
            code = code + chemical_property.get(aa, [0, 0, 0])
        encodings.append(code)
    return np.array(encodings,dtype=float)
def TriNcleotideComposition(sequence, base):
    trincleotides = [nn1 + nn2 + nn3 for nn1 in base for nn2 in base for nn3 in base]
    tnc_dict = {}
    for triN in trincleotides:
        tnc_dict[triN] = 0
    for i in range(len(sequence) - 2):
        tnc_dict[sequence[i:i + 3]] += 1
    for key in tnc_dict:
       tnc_dict[key] /= (len(sequence) - 2)
    return tnc_dict

def PseEIIP(fastas, **kw):
    for i in fastas:
        if re.search('[^ACGT-]', i[1]):
            print('Error: illegal character included in the fasta sequences, only the "ACGT-" are allowed by this PseEIIP scheme.')
            return 0

    base = 'ACGT'

    EIIP_dict = {
        'A': 0.1260,
        'C': 0.1340,
        'G': 0.0806,
        'T': 0.1335,
    }

    trincleotides = [nn1 + nn2 + nn3 for nn1 in base for nn2 in base for nn3 in base]
    EIIPxyz = {}
    for triN in trincleotides:
        EIIPxyz[triN] = EIIP_dict[triN[0]] + EIIP_dict[triN[1]] + EIIP_dict[triN[2]]

    encodings = []
    #header = ['#', 'label'] + trincleotides
    #encodings.append(header)

    for i in fastas:
        name, sequence = i[0], re.sub('-', '', i[1])
        code = []
        trincleotide_frequency = TriNcleotideComposition(sequence, base)
        code = code + [EIIPxyz[triN] * trincleotide_frequency[triN] for triN in trincleotides]
        encodings.append(code)
    return np.array(encodings,dtype=float)
def CalculateMatrix(data, order):
    matrix = np.zeros((len(data[0]) - 2, 24))
    for i in range(len(data[0]) - 2): # position
        for j in range(len(data)):
            if re.search('-', data[j][i:i+3]):
                pass
            else:
                matrix[i][order[data[j][i:i+3]]] += 1
    return matrix


def PSTNPds(fastas, **kw):
    if check_fasta_with_equal_length == False:
        print('Error: for "PSTNP" encoding, the input fasta sequences should be with equal length. \n\n')
        return 0

    for i in fastas:
        if re.search('[^ACGT-]', i[1]):
            print('Error: illegal character included in the fasta sequences, only the "ACGT[U]" are allowed by PSTNPds encoding scheme.')
            return 0

    for i in fastas:
        i[1] = re.sub('T', 'A', i[1])
        i[1] = re.sub('G', 'C', i[1])

    encodings = []
    #header = ['#','label']
    #for pos in range(len(fastas[0][1])-2):
    #    header.append('Pos.%d' %(pos+1))
    #encodings.append(header)

    positive = []
    negative = []
    positive_key = []
    negative_key = []
    for i in fastas:
        if i[3] == 'training':
            if i[2] == '1':
                positive.append(i[1])
                positive_key.append(i[0])
            else:
                negative.append(i[1])
                negative_key.append(i[0])

    nucleotides = ['A', 'C','T','G']
    trinucleotides = [n1 + n2 + n3 for n1 in nucleotides for n2 in nucleotides for n3 in nucleotides]
    order = {}
    for i in range(len(trinucleotides)):
        order[trinucleotides[i]] = i

    matrix_po = CalculateMatrix(positive, order)
    matrix_ne = CalculateMatrix(negative, order)

    positive_number = len(positive)
    negative_number = len(negative)

    for i in fastas:
        if i[3] == 'training':
            name, sequence = i[0], i[1]
            code = []
            for j in range(len(sequence) - 2):
                if re.search('-', sequence[j: j + 3]):
                    code.append(0)
                else:
                    p_num, n_num = positive_number, negative_number
                    po_number = matrix_po[j][order[sequence[j: j+3]]]
                    if i[0] in positive_key and po_number > 0:
                        po_number -= 1
                        p_num -= 1
                    ne_number = matrix_ne[j][order[sequence[j: j+3]]]
                    if i[0] in negative_key and ne_number > 0:
                        ne_number -= 1
                        n_num -= 1
                    code.append(po_number/p_num - ne_number/n_num)
                    # print(sequence[j: j+3], order[sequence[j: j+3]], po_number, p_num, ne_number, n_num)
            encodings.append(code)
    return np.array(encodings,dtype=float)
def CalculateMatrix(data, order):
    matrix = np.zeros((len(data[0]) - 2, 64))
    for i in range(len(data[0]) - 2): # position
        for j in range(len(data)):
            if re.search('-', data[j][i:i+3]):
                pass
            else:
                matrix[i][order[data[j][i:i+3]]] += 1
    return matrix


def PSTNPss(fastas, **kw):
    if check_fasta_with_equal_length == False:
        print('Error: for "PSTNP" encoding, the input fasta sequences should be with equal length. \n\n')
        return 0

    for i in fastas:
        if re.search('[^ACGU-]', i[1]):
            print('Error: illegal character included in the fasta sequences, only the "ACGT[U]" are allowed by this encoding scheme.')
            return 0

    encodings = []
    #header = ['#', 'label']
    #for pos in range(len(fastas[0][1])-2):
    #    header.append('Pos.%d' %(pos+1))
    #encodings.append(header)

    # print(fastas[0])

    positive = []
    negative = []
    positive_key = []
    negative_key = []
    for i in fastas:
        if i[3] == 'training':
            if i[2] == '1':
                positive.append(i[1])
                positive_key.append(i[0])
            else:
                negative.append(i[1])
                negative_key.append(i[0])

    nucleotides = ['A', 'C', 'G', 'T']
    trinucleotides = [n1 + n2 + n3 for n1 in nucleotides for n2 in nucleotides for n3 in nucleotides]
    order = {}
    for i in range(len(trinucleotides)):
        order[trinucleotides[i]] = i

    matrix_po = CalculateMatrix(positive, order)
    matrix_ne = CalculateMatrix(negative, order)

    positive_number = len(positive)
    negative_number = len(negative)
    for i in fastas:
        if i[3] == 'training':
            name, sequence = i[0], i[1]
            code = []
            for j in range(len(sequence) - 2):
                if re.search('-', sequence[j: j+3]):
                    code.append(0)
                else:
                    p_num, n_num = positive_number, negative_number
                    po_number = matrix_po[j][order[sequence[j: j+3]]]
                    if i[0] in positive_key and po_number > 0:
                        po_number -= 1
                        p_num -= 1
                    ne_number = matrix_ne[j][order[sequence[j: j+3]]]
                    if i[0] in negative_key and ne_number > 0:
                        ne_number -= 1
                        n_num -= 1
                    code.append(po_number/p_num - ne_number/n_num)
                    # print(sequence[j: j+3], order[sequence[j: j+3]], po_number, p_num, ne_number, n_num)
            encodings.append(code)

    return np.array(encodings,dtype=float)
def kmerArray(sequence, k):
    kmer = []
    for i in range(len(sequence) - k + 1):
        kmer.append(sequence[i:i + k])
    return kmer


def RC(kmer):
    myDict = {
        'A': 'T',
        'C': 'G',
        'G': 'C',
        'T': 'A'
    }
    return ''.join([myDict[nc] for nc in kmer[::-1]])


def generateRCKmer(kmerList):
    rckmerList = set()
    myDict = {
        'A': 'T',
        'C': 'G',
        'G': 'C',
        'T': 'A'
    }
    for kmer in kmerList:
        rckmerList.add(sorted([kmer, ''.join([myDict[nc] for nc in kmer[::-1]])])[0])
    return sorted(rckmerList)


def RCKmer(fastas, k, upto=False, normalize=True, **kw):
    encoding = []
    header = []
    NA = 'ACGT'

    if k < 1:
        print('Error: the k-mer value should larger than 0.')
        return 0

    if upto == True:
        for tmpK in range(1, k + 1):
            tmpHeader = []
            for kmer in itertools.product(NA, repeat=tmpK):
                tmpHeader.append(''.join(kmer))
            header = header + generateRCKmer(tmpHeader)
        myDict = {}
        for kmer in header[2:]:
            rckmer = RC(kmer)
            if kmer != rckmer:
                myDict[rckmer] = kmer
        #encoding.append(header)
        for i in fastas:
            name, sequence = i[0], re.sub('-', '', i[1])
            count = Counter()
            for tmpK in range(1, k + 1):
                kmers = kmerArray(sequence, tmpK)
                for j in range(len(kmers)):
                    if kmers[j] in myDict:
                        kmers[j] = myDict[kmers[j]]
                count.update(kmers)
                if normalize == True:
                    for key in count:
                        if len(key) == tmpK:
                            count[key] = count[key] / len(kmers)
            code = []
            for j in range(2, len(header)):
                if header[j] in count:
                    code.append(count[header[j]])
                else:
                    code.append(0)
            encoding.append(code)
    else:
        tmpHeader = []
        for kmer in itertools.product(NA, repeat=k):
            tmpHeader.append(''.join(kmer))
        header = header + generateRCKmer(tmpHeader)
        myDict = {}
        for kmer in header[2:]:
            rckmer = RC(kmer)
            if kmer != rckmer:
                myDict[rckmer] = kmer

        #encoding.append(header)
        for i in fastas:
            name, sequence = i[0], re.sub('-', '', i[1])
            kmers = kmerArray(sequence, k)
            for j in range(len(kmers)):
                if kmers[j] in myDict:
                    kmers[j] = myDict[kmers[j]]
            count = Counter()
            count.update(kmers)
            if normalize == True:
                for key in count:
                    count[key] = count[key] / len(kmers)
            code = []
            for j in range(2, len(header)):
                if header[j] in count:
                    code.append(count[header[j]])
                else:
                    code.append(0)
            encoding.append(code)
    return np.array(encoding,dtype=float)
def TNC(fastas, **kw):
    AA = 'ACGT'
    encodings = []
    triPeptides = [aa1 + aa2 + aa3 for aa1 in AA for aa2 in AA for aa3 in AA]
    #header = ['#', 'label'] + triPeptides
    #encodings.append(header)

    AADict = {}
    for i in range(len(AA)):
        AADict[AA[i]] = i

    for i in fastas:
        name, sequence = i[0], re.sub('-', '', i[1])
        code = []
        tmpCode = [0] * 64
        for j in range(len(sequence) - 3 + 1):
            tmpCode[AADict[sequence[j]] * 16 + AADict[sequence[j+1]]*4 + AADict[sequence[j+2]]] = tmpCode[AADict[sequence[j]] * 16 + AADict[sequence[j+1]]*4 + AADict[sequence[j+2]]] +1
        if sum(tmpCode) != 0:
            tmpCode = [i/sum(tmpCode) for i in tmpCode]
        code = code + tmpCode
        encodings.append(code)
    return np.array(encodings,dtype=float)
from sklearn.metrics import roc_auc_score

ALPHABET='ACGT'

def readDNAFasta(file):
	with open(file) as f:
		records = f.read()
	if re.search('>', records) == None:
		print('Error,the input DNA sequence must be fasta format.')
		sys.exit(1)
	records = records.split('>')[1:]
	myFasta = []
	for fasta in records:
		array = fasta.split('\n')
		name, sequence = array[0].split()[0], re.sub('[^ACGT-]', '-', ''.join(array[1:]).upper())
		myFasta.append([name, sequence])
	return myFasta

def MonoDiKGap_vector(input_data,g):   
    fastas=readDNAFasta(input_data)
    vector=[] 
    header=['#']
    for f in range((g)*32):
        header.append('MonoDi.'+str(f))
    vector.append(header)
    sample=[]
    for i in fastas:
        name, sequence = i[0], re.sub('-', '', i[1])
        sample = [name]
        each_vec=MonoDiKGap(sequence,g)
        sample=sample+each_vec
        vector.append(sample)
    return np.array(vector,dtype=float)
     
########################################################################

kw = {'order': 'ACGT'}
fastas = read_nucleotide_sequences('/kaggle/input/promotersdeep/P-Non.txt')
#feat0 = ANF(fastas, **kw)
feat0 = RCKmer(fastas, k=3,  **kw)
#feat1 = binary(fastas, **kw)
#feat2 = CKSNAP(fastas, 4, **kw)
#feat3 = DNC(fastas, **kw)
#feat4 = EIIP(fastas, **kw)
#feat3 = ENAC(fastas, window=5, **kw)
#feat4 = Kmer(fastas, k=2,  **kw)
#feat4 = NCP(fastas, **kw)
#feat8 = PseEIIP(fastas, **kw)
feat5 = PSTNPds(fastas, **kw)
feat6 = PSTNPss(fastas, **kw)
#feat1 = RCKmer(fastas, k=3,  **kw)
#feat5 = MonoDiKGap_vector(fastas,1)
#feat6 = MonoDiKGap_vector(fastas,2)
#feat7 = MonoDiKGap_vector(fastas,3)
#feat8 = MonoDiKGap_vector(fastas,4)
#feat7=TNC(fastas, **kw)
feat12=TNC(fastas, **kw)
data=pd.read_csv('/kaggle/input/promotersdeep/Train_promoter_vecs.csv',header=None)
feat9=np.array(data)
data=pd.read_csv('/kaggle/input/promotersdeep/1_interval_promoter_train_PSTNP.csv',header=None)
feat10=np.array(data)
data=pd.read_csv('/kaggle/input/promotersdeep/2_interval_promoter_train_PSTNP.csv',header=None)
feat11=np.array(data)
all_feat=np.hstack((feat0 ,feat9,feat6,feat5,feat12,feat10,feat11))#,feat6,feat7))


def precision(y_true, y_pred):
    # Calculates the precision
    true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
    predicted_positives = K.sum(K.round(K.clip(y_pred, 0, 1)))
    precision = true_positives / (predicted_positives + K.epsilon())
    return precision

def recall(y_true, y_pred):
    # Calculates the recall
    true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
    possible_positives = K.sum(K.round(K.clip(y_true, 0, 1)))
    recall = true_positives / (possible_positives + K.epsilon())
    return recall

def f1(test_Y, pre_test_y):
    #calculate the F1-score
    Precision = precision(test_Y, pre_test_y)
    Recall = recall(test_Y, pre_test_y)
    f1 = 2 * ((Precision * Recall) / (Precision + Recall + K.epsilon()))
    return f1 

def TP(test_Y,pre_test_y):
    #calculate numbers of true positive samples
    TP = K.sum(K.round(K.clip(test_Y * pre_test_y, 0, 1)))#TP
    return TP

def FN(test_Y,pre_test_y):
     #calculate numbers of false negative samples
    TP = K.sum(K.round(K.clip(test_Y * pre_test_y, 0, 1)))#TP
    P=K.sum(K.round(K.clip(test_Y, 0, 1)))
    FN = P-TP #FN=P-TP
    return FN

def TN(test_Y,pre_test_y):
    #calculate numbers of True negative samples
    TN=K.sum(K.round(K.clip((test_Y-K.ones_like(test_Y))*(pre_test_y-K.ones_like(pre_test_y)), 0, 1)))#TN
    return TN

def FP(test_Y,pre_test_y):
    #calculate numbers of False positive samples
    N = (-1)*K.sum(K.round(K.clip(test_Y-K.ones_like(test_Y), -1, 0)))#N
    TN=K.sum(K.round(K.clip((test_Y-K.ones_like(test_Y))*(pre_test_y-K.ones_like(pre_test_y)), 0, 1)))#TN
    FP=N-TN
    return FP

probas_cnn=[]
tprs_cnn = []
sepscore_cnn = []
def calculate_performace(test_num, pred_y, labels):
    tp = 0
    fp = 0
    tn = 0
    fn = 0
    for index in range(test_num):
        if labels[index] == 1:
            if labels[index] == pred_y[index]:
                tp = tp + 1
            else:
                fn = fn + 1
        else:
            if labels[index] == pred_y[index]:
                tn = tn + 1
            else:
                fp = fp + 1
    acc = float(tp + tn) / test_num
    precision = float(tp) / (tp + fp)
    sensitivity = float(tp) / (tp + fn)
    specificity = float(tn) / (tn + fp)
    MCC = float(tp * tn - fp * fn) / (np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)))
    return acc, precision, sensitivity, specificity, MCC
def transfer_label_from_prob(proba):
    label = [1 if val >= 0.5 else 0 for val in proba]
    return label
vs = 50
max_length=50
sepscore_cnn = []
all_prob = {}
all_prob[0] = []
all_y_pred=[]
def dnn_model(train_X, train_Y, test_X, test_Y, lr, epoch, batch_size):
    train_X = np.expand_dims(train_X, 2)
    test_X = np.expand_dims(test_X, 2)
    inputs = Input(shape = (train_X.shape[1], train_X.shape[2]))
    x = Conv1D(32, kernel_size = 3, strides = 1, padding = 'same', activation = 'relu')(inputs)
    x = MaxPooling1D(pool_size = 2, strides = 2, padding = 'valid')(x)
    x = Conv1D(32, kernel_size = 3, strides = 1, padding = 'same', activation = 'relu')(x)
  #  x = MaxPooling1D(pool_size = 2, strides = 2, padding = 'same')(x)
  #  x = LSTM(units=128,   name='Lstm',  return_sequences=True)(x)
    x = Dropout(0.5)(x)
    x = Bidirectional(LSTM(units=128, name='Bilstm',return_sequences=True))(x)
    x = Flatten()(x)
  #  x = Dropout(0.5)(x)
   # x = Dense(32, activation = 'relu',kernel_regularizer = l2(1e-7))(x)
   # x = Dense(16, activation = 'relu',kernel_regularizer = l2(1e-7))(x)
    x = Dense(8, activation = 'relu',kernel_regularizer = l2(1e-7))(x)
    predictions = Dense(1, activation = 'sigmoid')(x) 
    model = Model(inputs = inputs, outputs = predictions)
    print("model")
    model.compile(optimizer = 'Adam',
                  loss = 'mean_squared_error',
                  metrics = ['acc'])#,precision,recall,f1,TP,FN,TN,FP]),RMSProp, SGD,Nadam,
    print("compile")
    result = model.fit(train_X, train_Y, epochs = 80, batch_size = 64, validation_data = (test_X, test_Y), shuffle = True)
    model.save('GP50_model') #save model
    pre_test_y = model.predict(test_X, batch_size = 64)
    pre_train_y = model.predict(train_X, batch_size = 64)
    test_auc = metrics.roc_auc_score(test_Y, pre_test_y)
    train_auc = metrics.roc_auc_score(train_Y, pre_train_y)
    all_prob[0] = all_prob[0] + [val for val in pre_test_y]
    y_pred_LSTM = transfer_label_from_prob(pre_test_y)
    np.savetxt('yscore.csv', np.asarray(all_prob[0]), delimiter=',', fmt='%f')
    np.savetxt('TestLabel.csv', np.asarray(test_Y), delimiter=',', fmt='%f')
    acc, precision, sensitivity, specificity, MCC = calculate_performace(len(test_Y), y_pred_LSTM, test_Y)
    print("train_auc: ", train_auc)
    print("test_auc: ", test_auc)
    print('DeepAntibiotics:acc=%f,precision=%f,sensitivity=%f,specificity=%f,MCC=%f,roc_auc=%f'
          % (acc, precision, sensitivity, specificity, MCC, test_auc))
    sepscore_cnn.append([acc, precision, sensitivity, specificity, MCC,test_auc])

    return test_auc,sepscore_cnn,result 

data=all_feat
X=scale(data)
[m1,n1]=np.shape(data) 
print(m1,n1)
label1=np.ones((int(m1/2),1))
label2=np.zeros((int(m1/2),1))
Y=np.append(label1,label2) 
#training_data=np.reshape(data,(9126,83,6))
#Y = Y.reshape((Y.shape[0], -1))
print (X)
print ("X.shape: ", X.shape)
print ("Y.shape: ", Y.shape)

lr = 0.00001 #learning rate
epoch = 80
batch_size = 64
kf = KFold(n_splits = 10, shuffle = True, random_state = 5)
#kf = KFold(n_splits = 10, shuffle = False)
kf = kf.split(X)

test_aucs = []
for i, (train_fold, validate_fold) in enumerate(kf):
    print("\n\ni: ", i)
    #test_auc = dnn_model(X[train_fold], Y[train_fold], X[validate_fold], Y[validate_fold], lr, epoch, batch_size)
    #test_auc = test_aucs.append(test_auc)
    test_auc,sepscore_cnn,history = dnn_model(X[train_fold], Y[train_fold], X[validate_fold], Y[validate_fold], lr, epoch, batch_size)
    test_auc = test_aucs.append(test_auc)
    #sepscore_cnn=sepscore_cnn.append((sepscore_cnn))
scores=np.array(sepscore_cnn)
result1=np.mean(scores,axis=0)
H1=result1.tolist()
sepscore_cnn.append(H1)
result=sepscore_cnn
data_csv = pd.DataFrame(data=result)
data_csv.to_csv('E 80--BS 64 L 0.00001 +TCN.csv')
plt.plot(history.history['acc'])
plt.plot(history.history['val_acc'])
plt.plot(history.history['loss'])
plt.plot(history.history['val_loss'])
plt.title('Model accuracy')
plt.ylabel('Accuracy_loss')
plt.xlabel('Epoch')
# plt.legend(['Train', 'Test'], loc='upper left')
# plt.savefig("LSTM_RESULTS/Train_Val_accuracy.png")
# plt.figure(2)
# plt.title('Model_loss')
# plt.ylabel('Loss')
# plt.xlabel('Epochs')
plt.legend(['Train', 'test'], loc='best')
plt.savefig("Train_Test_loss_accuracy.png")
plt.show()
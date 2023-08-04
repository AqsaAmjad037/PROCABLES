#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import sys, os, platform
import re
pPath = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(pPath)
father_path = os.path.abspath(
    os.path.dirname(pPath) + os.path.sep + ".") + r'\pubscripts' if platform.system() == 'Windows' else os.path.abspath(
    os.path.dirname(pPath) + os.path.sep + ".") + r'/pubscripts'
sys.path.append(father_path)
import check_sequences
import read_fasta_sequences
import save_file
import pandas as pd
import  numpy as np
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
    if check_sequences.check_fasta_with_equal_length == False:
        print('Error: for "PSTNP" encoding, the input fasta sequences should be with equal length. \n\n')
        return 0

    for i in fastas:
        if re.search('[^ACGT-]', i[1]):
            print('Error: illegal character included in the fasta sequences, only the "ACGT[U]" are allowed by this encoding scheme.')
            return 0

    encodings = []
    header = []
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

    return encodings
if __name__ == '__main__':
    '''
    parser = argparse.ArgumentParser(usage="it's usage tip.",
                                     description="Generating CKSAAP feature vector for nucleotide sequences")
    parser.add_argument("--file", required=True, help="input fasta file")
    parser.add_argument("--gap", type=int, default=5, help="the k-space value for CKSNAP descriptor")
    parser.add_argument("--format", choices=['csv', 'tsv', 'svm', 'weka'], default='svm', help="the output format")
    parser.add_argument("--out", help="the generated descriptor file")
    args = parser.parse_args()
    output = args.out if args.out != None else 'encoding.txt'
    '''
    kw = {'order': 'ACGT'}
    fastas = read_fasta_sequences.read_nucleotide_sequences('Enhancer_nonEnhancer_Layer2_traing.txt')
    encodings = PSTNPss(fastas, **kw)
    save_file.save_file(encodings, 'csv')
    data_ = pd.DataFrame(data=encodings)
    data_.to_csv('Enhancer_nonEnhancer_Layer2_traing_PSTNPss.csv')

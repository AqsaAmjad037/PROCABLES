#!/usr/bin/env python
# coding: utf-8


from gensim.models import Word2Vec
from gensim.models.word2vec import LineSentence
import pandas as pd
import numpy as np
import os
import sys
import math
import random
import warnings
from sklearn import preprocessing
import sklearn.preprocessing
from gensim import corpora, models, similarities

#build corpus
def DNA2Sentence(dna, K):

    sentence = ""
    length = len(dna)

    for i in range(length - K + 1):
        sentence += dna[i: i + K] + " "

    #delete extra space
    sentence = sentence[0 : len(sentence) - 1]
    return sentence

def Get_Unsupervised(fname,gname,kmer):
	f = open(fname,'r')
	g = open(gname,'w')
	k = kmer
	for i in f:
		if '>' not in i:
			i = i.strip('\n').upper()
			line = DNA2Sentence(i,k)
			g.write(line+'\n')
	f.close()
	g.close()

Get_Unsupervised('positive_train.fasta','Train_Pos2Unsuper',2)#postive samples contained in 'pos.fasta';'pos2Un' is outputFile of corpus;'2' is the size of Kmer
Get_Unsupervised('negative_train.fasta','Train_Neg2Unsuper',2)#negative samples contained in 'neg.fasta';'neg2Un' is outputFile of corpus;'2' is the size of Kmer

#combine two corpus and generate final corpus
with open('Train_2Unsuper','ab') as f:
        f.write(open('Train_Pos2Unsuper','rb').read())
        f.write(open('Train_Neg2Unsuper','rb').read())

#get model
def getWord_model(word,num_features,min_count):
	word_model = ""
	if not os.path.isfile("Train_model1"):
		sentence = LineSentence("Train_2Unsuper",max_sentence_length = 15000)
		print ("Start Training Word2Vec model...")
		# Set values for various parameters
		num_features = int(num_features)	  # Word vector dimensionality
		min_word_count = int(min_count)	  # Minimum word count
		num_workers = 20		 # Number of threads to run in parallel
		context = 20			# Context window size
		downsampling = 1e-3	 # Downsample setting for frequent words

		# Initialize and train the model
		print ("Training Word2Vec model...")
		word_model = Word2Vec(sentence, workers=num_workers, size=num_features, min_count=min_word_count, window=context, sample=downsampling, seed=1,iter = 50)
		word_model.init_sims(replace=False)
		word_model.save("Train_model1")
		#print word_model.most_similar("CATAGT")
	else:
		print ("Loading Word2Vec model...")
		word_model = Word2Vec.load("Train_model1")
		#word_model.init_sims(replace=True)
	return word_model

getWord_model(2,200,1)

def combine(Posfile, Negfile, Combfile):
    f1 = open(Posfile)
    f2 = open(Negfile)
    g = open(Combfile,'w')
    g.write('lable\tseq\n')
    for i in f1:
        if '>'not in i:
            g.write('1\t'+i)
    for i in f2:
        if '>'not in i:
            g.write('0\t'+i)
    f1.close()
    f2.close()
    g.close()

combine('positive_train.fasta','negative_train.fasta','Train_Comb.fasta')#'pos.fasta' contains positive samples with fasta format; 'neg.fasta' contains negative samples with fasta format; 'all.fasta' is a combination file of positive and negative samples.

#obtain feature file with .npy format
def getDNA_split(DNAdata,word):
	DNAlist1 = []
	#DNAlist2 = []
	counter = 0
	for DNA in DNAdata["seq"]:
		#if counter % 100 == 0:
			#print ("DNA %d of %d\r" % (counter, 2*len(DNAdata)))
			#sys.stdout.flush()

		DNA = str(DNA).upper()
		DNAlist1.append(DNA2Sentence(DNA,word).split(" "))#[['ACG', 'CGT', 'GTC'],['ACG', 'CGT', 'GTC'],['ACG', 'CGT', 'GTC']]

		counter += 1
	return DNAlist1

#	print()
#	return DNAlist1,DNAlist2

#def getAvgFeatureVecs(DNAdata1,DNAdata2,model,num_features):
def getAvgFeatureVecs(DNAdata1,model,num_features):
	counter = 0
	DNAFeatureVecs = np.zeros((len(DNAdata1),num_features), dtype="float32")
	for DNA in DNAdata1:
		if counter % 1000 == 0:
			print ("DNA %d of %d\r" % (counter, len(DNAdata1)))
			sys.stdout.flush()

		DNAFeatureVecs[counter][0:num_features] = np.mean(model[DNA],axis = 0)
		counter += 1
	print()
	counter = 0
	return DNAFeatureVecs

def DNA2Sentence(dna, K):
	sentence = ""
	length = len(dna)

	for i in range(length - K + 1):
		sentence += dna[i: i + K] + " "

	#delete extra space
	sentence = sentence[0 : len(sentence) - 1]
	return sentence

data = pd.read_csv('Train_Comb.fasta',sep = "\t",error_bad_lines=False)
#datawords1,datawords2 = getDNA_split(data,4)
datawords1 = getDNA_split(data,2)#size of kmer
#print(datawords1)

word_model = Word2Vec.load("Train_model1")#'model_2' is above generated in 'get model'
#dataDataVecs = getAvgFeatureVecs(datawords1,datawords2,word_model,100)
dataDataVecs = getAvgFeatureVecs(datawords1,word_model,200)
print (dataDataVecs.shape)
np.save("222_vecs.npy",dataDataVecs) #outputFile

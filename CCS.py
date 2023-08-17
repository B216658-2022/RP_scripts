#user/bin/python

import re
import math
import random
import shutil
import numpy as np
import numpy.matlib
import pandas as pd
import os, sys, subprocess
import matplotlib.pyplot as plt
from pandas import DataFrame,Series
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import silhouette_score, calinski_harabasz_score, davies_bouldin_score

def get_chrcla(filename):
	os.system(f'mkdir {filename}.each.fold')
	print(f'Processing {filename}, please wait')
	df = pd.read_csv(filename, sep='\t', header=None)
	a = df[0]
	b = a.unique()
	c = []
	for i in b:
		match = re.search(r'(.+?)_(.+?)', i)
		if match is not None:
			c.append(match.group(1))
		else:
			c.append(i)
	d = set(c)
	e = sorted(d, key=lambda x: (int(x[3:]) if x[3:].isdigit() else float('inf'), x))
	chr = choice(filename,e)
	f = []
	for j in chr:
		pattern = r'{}$|{}\D'.format(j, j)
		for k in b:
			match = re.search(pattern, k)
			if match is not None:
				f.append(k)
	for m in range(len(f)):
		row_tem = df[df[0] == f[m]]
		row_tem.to_csv(f'./{filename}.each.fold/narr_{f[m]}.txt', sep='\t', index=False, header=False)
	return chr


def choice(filename,chr_all):
	print(f'Classifying the {filename} file based on chromosome, it will be divided into the following categories: {chr_all}')	
	print(f'whether only keep autosome and sex chromosome.')
	choice = input("Please enter Y/N \n ")
	if choice == "Y":
		chr_sel = []
		for i in chr_all:
			match = re.search(r'chr(\d+|X|Y)',i)
			if match is not None:
				chr_sel.append(i)
		chr = chr_sel
	else :
		chr = chr_all
	print(f'Classifying the {filename} file based on chromosome, it will be divided into the following categories: {chr}')
	print(f'Classifying {filename}, please wait.')
	return chr


def data_matrix(file_name):
	matrix = []
	with open(file_name) as my_file:
		for eachline in my_file:
			fields = eachline.strip().split('\t')
			matrix.append(float(fields[6]))
	matrix = np.array(matrix)
	print(len(matrix))
	return matrix

def compare(X):
	clustering = AgglomerativeClustering(linkage='average').fit(X)
	labels = clustering.labels_
	label_uni = np.unique(labels)
	DB_HC = davies_bouldin_score(X, labels)
        return labels,DB_HC

def filt(data):
	random.shuffle(data)
	X = np.array(data).reshape(-1, 1)
	labels,DB = compare(X)
	s1 = pd.Series(data)
	s2 = pd.Series(labels)
	df = pd.DataFrame( { 'Data' : s1, 'Label' : s2, })
	# bulid datafram to story the data	
	label_uni = np.unique(labels)
	print()
	value_label={}
	# build a dic write down the labels and find out which label is represent the low peaks.
	for i in range(len(label_uni)):
		value = df[df['Label']==label_uni[i]]
		min_values = value.min().values[0]
		value_label[f'{label_uni[i]}'] = min_values
	# accodning to the labels to divide all the data in two groups: high_peak and unsure_peak.
	min_label = min(value_label, key=value_label.get)
	max_label = max(value_label, key=value_label.get)
	high_peak_df = df[df['Label']==int(max_label)]
	high_peak = np.array(high_peak_df['Data'])
	unsure_peak_df = df[df['Label']==int(min_label)]
	unsure_peak = np.array(unsure_peak_df['Data'])
	plt.figure(figsize=(60, 4))
	d = range(len(data)+1)
	plt.plot(d[1:],data)
	plt.ylabel('high')
	plt.xlabel('position')
	plt.show()

	plt.figure(figsize=(60, 4))
	plt.scatter(np.arange(len(X)), X, c=labels)
	plt.show()
	return high_peak_df,high_peak,unsure_peak_df,unsure_peak,DB


def circle(matrix,filename):
	matrix = np.sort(matrix)
	data = matrix[-30000:]
	merged_df = pd.DataFrame()
	high_peak_df,high_peak,unsure_peak_df,unsure_peak,DB = filt(data)
	print(DB)
	n = 0
	if DB > 0.50:
		merged_df = pd.concat([high_peak_df, unsure_peak_df], ignore_index=True)
	else:
		while DB <= 0.50:
			n += len(high_peak)
			a = -(30000+n)
			b = -n
			matrix = np.sort(matrix)
			c = matrix[a:b]
			d = set(unsure_peak)
			if len(unsure_peak) <= 2:
				merged_df = pd.concat([merged_df, high_peak_df, unsure_peak_df], ignore_index=True)
				return merged_df
			elif len(d) == 1 :
				merged_df = pd.concat([merged_df, high_peak_df, unsure_peak_df], ignore_index=True)
				return merged_df
			else :
				merged_df = pd.concat([merged_df, high_peak_df], ignore_index=True)
				high_peak_df,high_peak,unsure_peak_df,unsure_peak,DB = filt(c)
				print(DB)

	return merged_df


def highpeaks(chr,filename):
	#os.system(f'mkdir ./{filename}.fold') 
	os.system(f'mkdir ./{filename}.each.fold/peaks.file')
	for i in chr:
		print(f'Finding peaks for {i}.')
		matrix = data_matrix(f'./{filename}.each.fold/narr_{i}.txt')
		if len(matrix)<=2:
			row_tem = pd.read_csv(f'./{filename}.each.fold/narr_{i}.txt', sep='\t', header=None)
		else:
			merged_df = circle(matrix,f'./{filename}.each.fold/narr_{i}.txt')
			min_value = merged_df['Data'].min()
			df = pd.read_csv(f'./{filename}.each.fold/narr_{i}.txt', sep='\t', header=None)
			row_tem = df[df[6] >= min_value]
		row_tem.to_csv(f'./{filename}.each.fold/peaks.file/high_peaks_{i}.txt', sep='\t', index=False, header=False)
		row_tem.to_csv(f'./{filename}.each.fold/peaks.file/high_peaks_all.txt', sep='\t', mode='a', index=False, header=False)


### Fill in the name of the narrowPeak file that needs to be analyzed, and then run the script to complete it
                

chr = get_chrcla('file_name')
highpeaks(chr,'file_name') 

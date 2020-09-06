import scipy
import scipy.stats
import sys
import re
import os
import numpy as np
import math
from operator import itemgetter
from scipy.spatial.distance import squareform, pdist
from scipy.stats import percentileofscore
#import smfishHmrf.reader as reader
import pandas as pd
import logging

def get_distance_per_FD_2(mr_dissimilarity_FD, num_cell, clust, outcome=[1,2]):
	c1 = np.where(clust==1)[0]
	c2 = np.where(clust==2)[0]
	within_dist = mr_dissimilarity_FD[np.ix_(c1, c1)]
	across_dist = mr_dissimilarity_FD[np.ix_(c1, c2)]
	mm_vec = (np.sum(within_dist, axis=1) - within_dist.diagonal()) / float(within_dist.shape[0] - 1)
	mn_vec = np.mean(across_dist, axis=1)
	sil_vec = (mn_vec - mm_vec)/np.max(np.concatenate(([mn_vec], [mm_vec])), axis=0)
	avg_clust1_sil = np.mean(sil_vec)
	return avg_clust1_sil

#accepts a similarity matrix
def get_distance_per_FD_3(similarity_FD, num_cell, clust, outcome=[1,2]):
	c1 = np.where(clust==1)[0]
	within_dist = similarity_FD[np.ix_(c1, c1)]
	mm_vec = (np.sum(within_dist, axis=1) - within_dist.diagonal()) / float(within_dist.shape[0] - 1)
	sum_rbp = np.mean(mm_vec)
	return sum_rbp

def rank_transform_matrix(mat, rbp_p = 0.99, reverse=True, matrix_type="dissim", logger=logging):
	dim1 = mat.shape[0]
	dim2 = mat.shape[1]
	rank_forward = np.empty([dim1, dim2])

	#sys.stderr.write("Start ranking forward...\n")
	logger.info("Start ranking forward...")
	#sys.stderr.flush()

	for c1 in range(dim1):
		rd = scipy.stats.rankdata(mat[c1,:])
		if reverse==True:
			rd = dim2 - rd + 1
		rank_forward[c1, :] = rd
		if c1%1000==0:
			#sys.stderr.write("Done %d of %d\n" % (c1, dim1))
			logger.debug("Done %d of %d" % (c1, dim1))
			#sys.stderr.flush()

	#sys.stderr.write("Finished ranking forward...\n")
	#sys.stderr.flush()

	rank_backward = np.empty([dim1, dim2])
	#sys.stderr.write("Start ranking backward...\n")
	logger.info("Start ranking backward...")
	#sys.stderr.flush()

	for c1 in range(dim2):
		rd = scipy.stats.rankdata(mat[:,c1])
		if reverse==True:
			rd = dim1 - rd + 1
		rank_backward[:, c1] = rd
		if c1%1000==0:
			#sys.stderr.write("Done %d of %d\n" % (c1, dim2))
			logger.info("Done %d of %d" % (c1, dim2))
			#sys.stderr.flush()

	mutual_rank_rbp = np.empty([dim1, dim2])
	mutual_rank = np.empty([dim1, dim2])

	#sys.stderr.write("Calculate mutual rank...\n")
	logger.info("Calculate mutual rank...")
	#sys.stderr.flush()

	ma = np.sqrt(np.multiply(rank_forward, rank_backward))
	if matrix_type=="dissim":
		#sys.stderr.write("Calculate exponential transform...\n")
		logger.info("Calculate exponential transform...")
		#sys.stderr.flush()
		dissimilarity = np.subtract(1, np.power(rbp_p, np.subtract(ma, 1)))
		#sys.stderr.write("Finished exponential transform...\n")
		#sys.stderr.flush()
		#sys.stderr.write("Finished dissimilarity...\n")
		#sys.stderr.flush()
		return dissimilarity
	if matrix_type=="sim":
		#sys.stderr.write("Calculate exponential transform...\n")
		logger.info("Calculate exponential transform...")
		#sys.stderr.flush()
		similarity = np.power(rbp_p, np.subtract(ma, 1))
		#sys.stderr.write("Finished exponential transform...\n")
		#sys.stderr.flush()
		#sys.stderr.write("Finished similarity...\n")
		#sys.stderr.flush()
		return similarity	
	if verbose:
		sys.stderr.write("Error: matrix type is not recognized\n")
		sys.stderr.flush()
	return None

#matrix type is either dissim or sim
def random_pattern(matrix=None, matrix_type="dissim", num_cell=None, sizes=None, trials_per_gene=100, seed=-1, outdir="result.5000", run_gpd=False, logger=logging):
	ncell = num_cell
	ntrial = trials_per_gene
	if seed!=-1 and seed>=0:
		np.random.seed(seed)
	rand_sil = {}
	ind = 0
	repeats = {}
	while ind<len(sizes):
		aSize = sizes[ind]
		rand_clust = np.zeros((ncell), dtype="int32")
		rand_clust[0:aSize] = 1
		rand_clust[aSize:] = 2
		rand_sil.setdefault(aSize, [])
		#sys.stderr.write("Getting through size %d, %d of %d\n" % (aSize, ind, len(sizes)))
		logger.info("Getting through size %d, %d of %d" % (aSize, ind, len(sizes)))
		#sys.stderr.flush()
		repeats.setdefault(aSize, 0)
		repeats[aSize]+=1
		trials = []
		for i in range(trials_per_gene):
			a_trial = np.copy(rand_clust)
			np.random.shuffle(a_trial)
			trials.append(a_trial)
		if matrix_type=="dissim":
			for i in range(trials_per_gene):
				if i%100==0:
					#sys.stderr.write("%d\n" % i)
					logger.info("%d" % i)
					#sys.stderr.flush()
				r_avg_clust1_sil = get_distance_per_FD_2(matrix, ncell, trials[i], outcome=[1,2])
				rand_sil[aSize].append(r_avg_clust1_sil)
		elif matrix_type=="sim":
			for i in range(trials_per_gene):
				if i%100==0:
					#sys.stderr.write("%d\n" % i)
					logger.info("%d" % i)
					#sys.stderr.flush()
				r_avg_clust1_sil = get_distance_per_FD_3(matrix, ncell, trials[i], outcome=[1,2])
				rand_sil[aSize].append(r_avg_clust1_sil)
		fw = open("%s/%d" % (outdir, aSize), "w")
		for ix in rand_sil[aSize]:
			fw.write("%.10e\n" % ix)
		fw.close()
		re_run = False
		if run_gpd:
			os.system("cd '%s' && Rscript do_gpd.R %d" % (outdir, aSize))
			f = open("%s/par.%d" % (outdir, aSize))
			scale = float(f.readline().rstrip("\n").split("\t")[1])
			shape = float(f.readline().rstrip("\n").split("\t")[1])
			f.close()
			if shape>0:
				re_run = True
		if re_run:
			if repeats[aSize]>=10:
				#sys.stderr.write("Repeats reached, going to next one...\n")
				#sys.stderr.flush()
				logger.info("Repeats reached, going to next one...")
			else:
				rand_sil[aSize] = []
				ind -= 1
		ind+=1
	return rand_sil

def calc_silhouette_per_gene_approx(genes=None, expr=None, matrix=None, matrix_type="dissim", examine_top=0.1, seed=-1):
	if genes is None or expr is None or matrix is None:
		sys.stderr.write("Need genes, expr, similarity/dissimilarity matrix\n")
		sys.stderr.flush()
		return ;
	if seed!=-1 and seed>=0:
		np.random.seed(seed)
	sys.stderr.write("Started 2 " + "\n")
	sys.stderr.flush()
	sil = []
	ncell = expr.shape[1]
	ex = int((1.0-examine_top)*100.0)
	for ig,g in enumerate(genes):
		cutoff = np.percentile(expr[ig,:], ex)
		clust = np.zeros((ncell), dtype="int32")
		gt_eq = np.where(expr[ig,:]>=cutoff)[0]
		lt = np.where(expr[ig,:]<cutoff)[0]
		if gt_eq.shape[0]>int(ncell*examine_top):
			num_filter = gt_eq.shape[0] - int(ncell*examine_top)
			ss = np.random.choice(gt_eq, size=num_filter, replace=False)
			clust[gt_eq] = 1
			clust[lt] = 2
			clust[ss] = 2
		elif gt_eq.shape[0]<int(ncell*examine_top):
			num_filter = int(ncell*examine_top) - gt_eq.shape[0]
			ss = np.random.choice(lt, size=num_filter, replace=False)
			clust[gt_eq] = 1
			clust[lt] = 2
			clust[ss] = 1
		else:
			clust[gt_eq] = 1
			clust[lt] = 2
		if ig%100==0:
			sys.stderr.write("%s %d / %d\n" % (g, ig, len(genes)))
			sys.stderr.flush()
		if matrix_type=="dissim":
			avg_clust1_sil = get_distance_per_FD_2(matrix, ncell, clust, outcome=[1,2])
		elif matrix_type=="sim":
			avg_clust1_sil = get_distance_per_FD_3(matrix, ncell, clust, outcome=[1,2])
		sil.append((g, -1, avg_clust1_sil))
	res = []
	for ig,g in enumerate(genes):
		this_avg = sil[ig][1]
		this_sil = sil[ig][2]
		res.append((g, this_sil))
	#res.sort(lambda x,y:cmp(x[1], y[1]), reverse=True)
	res.sort(key=itemgetter(1), reverse=True)
	return res
	
def calc_silhouette_per_gene(genes=None, expr=None, matrix=None, matrix_type="dissim", examine_top=0.1, seed=-1, logger=logging):
	if genes is None or expr is None or matrix is None:
		sys.stderr.write("Need genes, expr, dissim/sim matrix\n")
		sys.stderr.flush()
		return ;
	if seed!=-1 and seed>=0:
		np.random.seed(seed)
	#sys.stderr.write("Started 2 " + "\n")
	logger.info("Started invoking silhouette function..")
	#sys.stderr.flush()
	sil = []
	ncell = expr.shape[1]
	ex = int((1.0-examine_top)*100.0)
	for ig,g in enumerate(genes):
		cutoff = np.percentile(expr[ig,:], ex)
		clust = np.zeros((ncell), dtype="int32")
		gt_eq = np.where(expr[ig,:]>=cutoff)[0]
		lt = np.where(expr[ig,:]<cutoff)[0]
		if cutoff==0:
			gt_eq = np.where(expr[ig,:]>cutoff)[0]
			lt = np.where(expr[ig,:]<=cutoff)[0]
		clust[gt_eq] = 1
		clust[lt] = 2
		if ig%100==0:
			#sys.stderr.write("%s %d / %d\n" % (g, ig, len(genes)))
			logger.info("%s %d / %d" % (g, ig, len(genes)))
			#sys.stderr.flush()
		if matrix_type=="dissim":
			avg_clust1_sil = get_distance_per_FD_2(matrix, ncell, clust, outcome=[1,2])
		elif matrix_type=="sim":
			avg_clust1_sil = get_distance_per_FD_3(matrix, ncell, clust, outcome=[1,2])
		sil.append((g, -1, avg_clust1_sil))
	res = []
	for ig,g in enumerate(genes):
		this_avg = sil[ig][1]
		this_sil = sil[ig][2]
		res.append((g, this_sil))
	res.sort(key=itemgetter(1), reverse=True)
	return res


def python_spatial_genes(spatial_locations, expression_matrix, metric = "euclidean", rbp_p = 0.95, examine_top = 0.3):
    
	Xcen =  spatial_locations
	mat = expression_matrix
	genes = []

	for g in range(mat.index.shape[0]):
		genes.append(str(mat.index[g]))
	expr = np.copy(mat.values)
    
	ncell = Xcen.shape[0] 
	sys.stderr.write("Calculate all pairwise Euclidean distance between cells using their physical coordinates\n")
	sys.stderr.flush()
	euc = squareform(pdist(Xcen, metric=metric))
	sys.stderr.write("Rank transform euclidean distance, and then apply exponential transform\n")
	sys.stderr.flush()
	dissim = rank_transform_matrix(euc, reverse=False, rbp_p=rbp_p)
	sys.stderr.write("Compute silhouette metric per gene\n")
	sys.stderr.flush()
	res = calc_silhouette_per_gene(genes=genes, expr=expr, dissim=dissim, examine_top=examine_top)
    
	return res
   
 

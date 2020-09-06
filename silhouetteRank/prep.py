import sys
import os
import re
import numpy as np
import subprocess
import math
import scipy
import silhouetteRank
import silhouetteRank.spatial_genes as spatial_genes
from shutil import copyfile
from operator import itemgetter
from scipy.spatial.distance import squareform, pdist
from scipy.stats import percentileofscore
from sklearn.metrics import roc_auc_score
import argparse
import logging

def read_matrix(f_expr="expression.txt", f_Xcen = "Xcen.good", logger=logging):
	#f_expr = "expression.txt"
	#f_Xcen = "Xcen.good"

	#sys.stderr.write("Reading gene expression...\n")
	logger.info("Reading gene expression...")
	#sys.stderr.flush()

	f = open(f_expr)
	h = f.readline().rstrip("\n").split("\t")[1:]
	num_cell = len(h)
	num_gene = 0
	for l in f:
		l = l.rstrip("\n")
		num_gene+=1
	f.close()

	expr = np.empty((num_gene, num_cell), dtype="float32")
	genes = []
	f = open(f_expr)
	f.readline()
	ig = 0
	for l in f:
		l = l.rstrip("\n")
		ll = l.split("\t")
		gene = ll[0]
		genes.append(gene)
		expr[ig,:] = [float(v) for v in ll[1:]]
		ig+=1
	f.close()

	#sys.stderr.write("Reading cell coordinates...\n")
	logger.info("Reading cell coordinates...")
	#sys.stderr.flush()
	f = open(f_Xcen)
	Xcen = np.empty((num_cell, 2), dtype="float32")
	f.readline()
	ic = 0
	for l in f:
		l = l.rstrip("\n")
		ll = l.split("\t")
		Xcen[ic,:] = [float(ll[1]), float(ll[2])]
		ic+=1
	f.close()
	return expr, genes, Xcen
 
def read_frequency(expr=None, genes=None, Xcen=None, frequency_file=None, read_from_file=True, outdir="", examine_top=0.05, num_query_sizes=10):
	num_cell = Xcen.shape[0]
	ncell = num_cell
	ex = int((1.0-examine_top)*100.0)
	pattern_size = {}
	for ig,g in enumerate(genes):
		cutoff = np.percentile(expr[ig,:], ex)
		clust = np.zeros((ncell), dtype="int32")
		gt_eq = np.where(expr[ig,:]>=cutoff)[0]
		lt = np.where(expr[ig,:]<cutoff)[0]
		t_size = gt_eq.shape[0]
		if cutoff==0:
			gt_eq = np.where(expr[ig,:]>cutoff)[0]
			lt = np.where(expr[ig,:]<=cutoff)[0]
			t_size = gt_eq.shape[0]
		clust[gt_eq] = 1
		clust[lt] = 2
		pattern_size.setdefault(t_size, 0)
		pattern_size[t_size]+=1

	freq = []
	if read_from_file==False:	
		fw = open("%s/gene.freq.good.txt" % outdir, "w")
		for s in pattern_size:
			fw.write("%d %d\n" % (pattern_size[s], s))
		fw.close()
		sizes = []
		for s in pattern_size:
			sizes.append(s)
			#frequency
			for i in range(pattern_size[s]):
				freq.append(s)
	else:
		f_frequency = frequency_file
		f = open(f_frequency)
		sizes = []
		for l in f:
			l = l.rstrip("\n")
			ll = l.split(" ")
			sizes.append(int(ll[1]))
			for i in range(int(ll[0])):
				freq.append(int(ll[1]))
		f.close()

	if len(sizes)<=num_query_sizes:
		#sys.stderr.write("Setting num_query_sizes to be %d instead of %d...\n" % (len(sizes), num_query_sizes))
		#sys.stderr.flush()
		num_query_sizes = len(sizes) 
	cmd = "cd '%s' && Rscript do_kmeans.R gene.freq.good.txt 1 %d 1000 km_centroid.txt km_label.txt" % (outdir, num_query_sizes)
	#print(cmd)
	os.system(cmd)

	f = open("%s/km_centroid.txt" % outdir)
	cc = []
	for l in f:
		l = l.rstrip("\n")
		ll = l.split(" ")
		cc.append(int(float(ll[1])))
	f.close()

	labels = []
	f = open("%s/km_label.txt" % outdir)
	for l in f:
		l = l.rstrip("\n")
		ll = l.split(" ")
		labels.append(int(ll[1]) - 1)
	f.close()

	by_cluster = {}
	for ind,c in enumerate(labels):
		by_cluster.setdefault(cc[c], [])
		by_cluster[cc[c]].append(freq[ind])

	if not os.path.isdir(outdir):
		os.mkdir(outdir)
	if not os.path.isdir("%s/good" % outdir):
		os.mkdir("%s/good" % outdir)
	
	for ind,v in enumerate(cc):
		fw = open("%s/good/%d" % (outdir, ind), "w")
		t_str = ";".join(["%d" % bc for bc in set(by_cluster[v])])
		fw.write("%d %d %s\n" % (1, v, t_str))
		fw.close()
	return list(by_cluster.keys())

def do_one(args):
	matrix_type = args.matrix_type
	rbp_ps = args.rbp_ps
	examine_tops = args.examine_tops


	result = subprocess.call("Rscript --version 2> /dev/null", shell=True)
	if result==127:
		sys.stderr.write("Rscript is not found")
		sys.stderr.flush()
		sys.exit(1)	

	if not os.path.isdir(args.output):
		os.mkdir(args.output)

	log_file = "%s/%s" % (args.output, args.log_file)
	logger = logging.getLogger("prep")
	logger.setLevel(logging.DEBUG)

	if not logger.hasHandlers():
		handler = logging.FileHandler(log_file)
		handler.setFormatter(logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s"))
		logger.addHandler(handler)
		if args.verbose:
			logger.addHandler(logging.StreamHandler())
	
	t_overwrite = args.overwrite_input_bin
	#if there is no binary input file, create one
	check_required = ["expr.npy", "Xcen.npy", "genes.npy"]
	for rbp_p in args.rbp_ps:
		check_required.append("t_matrix_%s_%.2f.npy" % (args.matrix_type, rbp_p))
	for cr in check_required:
		if not os.path.isfile("%s/%s" % (args.output, cr)):
			t_overwrite = True
			break

	expr, Xcen, genes, t_matrix = None, None, None, None
	if t_overwrite:
		expr, genes, Xcen = read_matrix(f_expr=args.expr, f_Xcen=args.centroid, logger=logger)
		#sys.stderr.write("Calculate all pairwise Euclidean distance between cells using their physical coordinates\n")
		logger.info("Calculate all pairwise Euclidean distance between cells using their physical coordinates")
		#sys.stderr.flush()
		euc = squareform(pdist(Xcen, metric="euclidean"))
		#sys.stderr.write("Rank transform euclidean distance, and then apply exponential transform\n")
		#sys.stderr.flush()
		for rbp_p in args.rbp_ps:
			#sys.stderr.write("For rbp_p %.2f: =====\n" % rbp_p)
			logger.info("For rbp_p %.2f:" % rbp_p)
			#sys.stderr.flush()
			t_matrix = spatial_genes.rank_transform_matrix(euc, reverse=False, rbp_p=rbp_p, matrix_type=matrix_type, logger=logger)
			np.save("%s/t_matrix_%s_%.2f.npy" % (args.output, args.matrix_type, rbp_p), t_matrix)
		np.save("%s/expr.npy" % args.output, expr)
		np.save("%s/Xcen.npy" % args.output, Xcen)
		np.save("%s/genes.npy" % args.output, genes)
	else:
		#sys.stderr.write("Using existing input binaries...\n")
		logger.info("Using existing input binaries...")
		#sys.stderr.flush()
		expr = np.load("%s/expr.npy" % args.output)
		Xcen = np.load("%s/Xcen.npy" % args.output)
		genes = np.load("%s/genes.npy" % args.output)
		#t_matrix = np.load("%s/t_matrix_%s_%.2f.npy" % (args.output, args.matrix_type, args.rbp_p))	
	
	ncell = Xcen.shape[0] 
	for rbp_p in args.rbp_ps:
		for examine_top in args.examine_tops:
			outdir = "%s/result_sim_5000_%.2f_%.3f" % (args.output, rbp_p, examine_top)
			if matrix_type=="dissim":
				outdir = "%s/result_5000_%.2f_%.3f" % (args.output, rbp_p, examine_top)
			if not os.path.isdir(outdir):
				os.mkdir(outdir)
			source_path = os.path.dirname(silhouetteRank.__file__)
			if not os.path.isfile("%s/do_kmeans.R" % outdir):
				copyfile("%s/do_kmeans.R" % source_path, "%s/do_kmeans.R" % outdir)
			if not os.path.isfile("%s/do_gpd.R" % outdir):
				copyfile("%s/do_gpd.R" % source_path, "%s/do_gpd.R" % outdir)
			if not os.path.isfile("%s/qval.R" % outdir):
				copyfile("%s/qval.R" % source_path, "%s/qval.R" % outdir)
			sizes = read_frequency(expr=expr, genes=genes, Xcen=Xcen, frequency_file=None, \
			read_from_file=False, outdir=outdir, examine_top=examine_top, num_query_sizes=args.query_sizes)
	
if __name__=="__main__":
	parser = argparse.ArgumentParser(description="evaluate.exact.2b.py: calculate silhouette score for randomly distributed spatial patterns", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("-x", "--file-expr", dest="expr", type=str, required=True, help="expression matrix. Will use input binary expr.npy (if exists) to speed up reading.")
	parser.add_argument("-c", "--file-centroid", dest="centroid", type=str, required=True, help="cell coordinate. Will use input binary Xcen.npy (if exists) to speed up reading.")
	parser.add_argument("-w", "--overwrite-input-binary", dest="overwrite_input_bin", action="store_true", help="overwrite input binaries")
	parser.add_argument("-r", "--rbp-ps", dest="rbp_ps", nargs="+", type=float, default=[0.95, 0.99], help="p parameter of RBP")
	parser.add_argument("-e", "--examine-tops", dest="examine_tops", nargs="+", type=float, default=[0.005, 0.010, 0.050, 0.100, 0.300], help="top proportion of cells per gene to be 1's (expressed)")
	parser.add_argument("-m", "--matrix-type", dest="matrix_type", type=str, choices=["sim", "dissim"], help="whether to calculate similarity matrix or dissimilarity matrix", default="dissim")
	parser.add_argument("-o", "--output-dir", dest="output", type=str, default=".", help="output directory")
	parser.add_argument("-l", "--log-filename", dest="log_file", type=str, default="master.prep.log", help="log file name (no path), will be generated in same directory as --output-dir")
	parser.add_argument("-v", "--verbose", dest="verbose", action="store_true", help="print verbose messages to console")
	parser.add_argument("-q", "--query-sizes", dest="query_sizes", type=int, default=10, help="query sizes (advanced user setting)")

	args = parser.parse_args()
	do_one(args)

import sys
import os
import re
import numpy as np
import math
import scipy
import silhouetteRank
import silhouetteRank.spatial_genes as spatial_genes
from operator import itemgetter
from scipy.spatial.distance import squareform, pdist
from scipy.stats import percentileofscore
from sklearn.metrics import roc_auc_score
import argparse
import logging

def read(f_expr="expression.txt", f_Xcen="Xcen.good", logger=logging):
	#f_expr = "expression.txt"
	#f_Xcen = "Xcen.good"
	#sys.stderr.write("Reading gene expression...\n")
	#sys.stderr.flush()
	logger.info("Reading gene expression...")
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
	#sys.stderr.flush()
	logger.info("Reading cell coordinates...")
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
	return expr, Xcen, genes

def do_one(args):
	matrix_type = args.matrix_type
	rbp_p = args.rbp_p

	if not os.path.isdir(args.output):
		os.mkdir(args.output)

	logdir="%s/logs" % args.output
	if not os.path.isdir(logdir):
		os.mkdir(logdir)
	log_file = "%s/real_%.2f_%.3f.out" % (logdir, args.rbp_p, args.examine_top)
	logger = logging.getLogger("real_%.2f_%.3f" % (args.rbp_p, args.examine_top))
	logger.setLevel(logging.DEBUG)

	if not logger.hasHandlers():
		handler = logging.FileHandler(log_file, "w")
		handler.setFormatter(logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s"))
		logger.addHandler(handler)
		if args.verbose:
			logger.addHandler(logging.StreamHandler())

	t_overwrite = args.overwrite_input_bin
	#if there is no binary input file, create one
	check_required = ["expr.npy", "Xcen.npy", "genes.npy", "t_matrix_%s_%.2f.npy" % (args.matrix_type, args.rbp_p)]
	for cr in check_required:
		if not os.path.isfile("%s/%s" % (args.output, cr)):
			t_overwrite = True
			break

	expr, Xcen, genes, t_matrix = None, None, None, None
	if t_overwrite:
		expr, Xcen, genes = read(f_expr=args.expr, f_Xcen=args.centroid, logger=logger)
		logger.info("Calculate all pairwise Euclidean distance between cells using their physical coordinates")
		#sys.stderr.flush()
		euc = squareform(pdist(Xcen, metric="euclidean"))
		logger.info("Rank transform euclidean distance, and then apply exponential transform")
		#sys.stderr.flush()
		t_matrix = spatial_genes.rank_transform_matrix(euc, reverse=False, rbp_p=rbp_p, matrix_type=matrix_type, logger=logger)
		np.save("%s/t_matrix_%s_%.2f.npy" % (args.output, args.matrix_type, args.rbp_p), t_matrix)
		np.save("%s/expr.npy" % args.output, expr)
		np.save("%s/Xcen.npy" % args.output, Xcen)
		np.save("%s/genes.npy" % args.output, genes)
	else:
		logger.info("Using existing input binaries...")
		#sys.stderr.flush()
		expr = np.load("%s/expr.npy" % args.output)
		Xcen = np.load("%s/Xcen.npy" % args.output)
		genes = np.load("%s/genes.npy" % args.output)
		t_matrix = np.load("%s/t_matrix_%s_%.2f.npy" % (args.output, args.matrix_type, args.rbp_p))	

	logger.info("Compute silhouette metric per gene")
	#sys.stderr.flush()
	examine_top = args.examine_top
	res = spatial_genes.calc_silhouette_per_gene(genes=genes, expr=expr, matrix=t_matrix, matrix_type=matrix_type, examine_top=examine_top, logger=logger)
	if matrix_type=="sim":
		f_name = "%s/silhouette.sim.exact.rbp.%.2f.top.%.3f.txt" % (args.output, rbp_p, examine_top)
	else:
		f_name = "%s/silhouette.exact.rbp.%.2f.top.%.3f.txt" % (args.output, rbp_p, examine_top)
		
	fw = open(f_name, "w")	
	for ind,v in enumerate(res):
		fw.write("%d\t%s\t%.10f\n" % (ind, v[0], v[1]))
	fw.close()

def main():
	parser = argparse.ArgumentParser(description="evaluate.2b.py: calculate silhouette score for spatial patterns", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("-x", "--file-expr", dest="expr", type=str, required=True, help="expression matrix. Will use input binary expr.npy (if exists) to speed up reading.")
	parser.add_argument("-c", "--file-centroid", dest="centroid", type=str, required=True, help="cell coordinate. Will use input binary Xcen.npy (if exists) to speed up reading.")
	parser.add_argument("-w", "--overwrite-input-binary", dest="overwrite_input_bin", action="store_true", help="overwrite input binaries")
	parser.add_argument("-e", "--examine-top", dest="examine_top", type=float, default=0.05, help="top proportion of cells per gene to be 1's (expressed)")

	parser.add_argument("-r", "--rbp-p", dest="rbp_p", type=float, default=0.95, help="p parameter of RBP")
	parser.add_argument("-m", "--matrix-type", dest="matrix_type", type=str, choices=["sim", "dissim"], help="whether to calculate similarity matrix or dissimilarity matrix", default="dissim")
	parser.add_argument("-o", "--output-dir", dest="output", type=str, default=".", help="output directory")
	parser.add_argument("-v", "--verbose", dest="verbose", action="store_true", help="print verbose messages to console")
	
	args = parser.parse_args()
	do_one(args)

if __name__=="__main__":
	main()

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
	logger = logging.getLogger("prep_fast")
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
		logger.info("Calculate all pairwise Euclidean distance between cells using their physical coordinates")
		euc = squareform(pdist(Xcen, metric="euclidean"))
		for rbp_p in args.rbp_ps:
			logger.info("For rbp_p %.2f:" % rbp_p)
			t_matrix = spatial_genes.rank_transform_matrix(euc, reverse=False, rbp_p=rbp_p, matrix_type=matrix_type, logger=logger)
			np.save("%s/t_matrix_%s_%.2f.npy" % (args.output, args.matrix_type, rbp_p), t_matrix)
		np.save("%s/expr.npy" % args.output, expr)
		np.save("%s/Xcen.npy" % args.output, Xcen)
		np.save("%s/genes.npy" % args.output, genes)
	else:
		logger.info("Using existing input binaries...")
		expr = np.load("%s/expr.npy" % args.output)
		Xcen = np.load("%s/Xcen.npy" % args.output)
		genes = np.load("%s/genes.npy" % args.output)
	
	ncell = Xcen.shape[0] 
	for rbp_p in args.rbp_ps:
		for examine_top in args.examine_tops:
			outdir = "%s/result_fast_sim_5000_%.2f_%.3f" % (args.output, rbp_p, examine_top)
			if matrix_type=="dissim":
				outdir = "%s/result_fast_5000_%.2f_%.3f" % (args.output, rbp_p, examine_top)
			if not os.path.isdir(outdir):
				os.mkdir(outdir)
			source_path = os.path.dirname(silhouetteRank.__file__)
			if not os.path.isfile("%s/do_kmeans.R" % outdir):
				copyfile("%s/do_kmeans.R" % source_path, "%s/do_kmeans.R" % outdir)
			if not os.path.isfile("%s/do_gpd.R" % outdir):
				copyfile("%s/do_gpd.R" % source_path, "%s/do_gpd.R" % outdir)
			if not os.path.isfile("%s/qval.R" % outdir):
				copyfile("%s/qval.R" % source_path, "%s/qval.R" % outdir)
	
if __name__=="__main__":
	parser = argparse.ArgumentParser(description="prep_fast.py", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("-x", "--file-expr", dest="expr", type=str, required=True, help="expression matrix. Will use input binary expr.npy (if exists) to speed up reading.")
	parser.add_argument("-c", "--file-centroid", dest="centroid", type=str, required=True, help="cell coordinate. Will use input binary Xcen.npy (if exists) to speed up reading.")
	parser.add_argument("-w", "--overwrite-input-binary", dest="overwrite_input_bin", action="store_true", help="overwrite input binaries")
	parser.add_argument("-r", "--rbp-ps", dest="rbp_ps", nargs="+", type=float, default=[0.95, 0.99], help="p parameter of RBP")
	parser.add_argument("-e", "--examine-tops", dest="examine_tops", nargs="+", type=float, default=[0.005, 0.010, 0.050, 0.100, 0.300], help="top proportion of cells per gene to be 1's (expressed)")
	parser.add_argument("-m", "--matrix-type", dest="matrix_type", type=str, choices=["sim", "dissim"], help="whether to calculate similarity matrix or dissimilarity matrix", default="dissim")
	parser.add_argument("-o", "--output-dir", dest="output", type=str, default=".", help="output directory")
	parser.add_argument("-l", "--log-filename", dest="log_file", type=str, default="master.prep.fast.log", help="log file name (no path), will be generated in same directory as --output-dir")
	parser.add_argument("-v", "--verbose", dest="verbose", action="store_true", help="print verbose messages to console")

	args = parser.parse_args()
	do_one(args)

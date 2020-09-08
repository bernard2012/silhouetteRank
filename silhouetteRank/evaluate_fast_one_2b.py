import sys
import os
import re
import numpy as np
import subprocess
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

def do_one(args):
	matrix_type = args.matrix_type
	rbp_p = args.rbp_p
	examine_top = args.examine_top
	result = subprocess.call("Rscript --version 2> /dev/null", shell=True)
	if result==127:
		sys.stderr.write("Rscript is not found")
		sys.stderr.flush()
		sys.exit(1)	

	if not os.path.isdir(args.output):
		os.mkdir(args.output)
	
	logdir="%s/logs.fast" % args.output
	if not os.path.isdir(logdir):
		os.mkdir(logdir)
	log_file = "%s/%.2f_%.3f.out" % (logdir, args.rbp_p, args.examine_top)
	logger = logging.getLogger("random_fast_%.2f_%.3f" % (args.rbp_p, args.examine_top))
	logger.setLevel(logging.DEBUG)

	if not logger.hasHandlers():
		handler = logging.FileHandler(log_file, "w")
		handler.setFormatter(logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s"))
		logger.addHandler(handler)
		if args.verbose:
			logger.addHandler(logging.StreamHandler())

	check_required = ["expr.npy", "Xcen.npy", "genes.npy", "t_matrix_%s_%.2f.npy" % (args.matrix_type, args.rbp_p)]
	for cr in check_required:
		if not os.path.isfile("%s/%s" % (args.output, cr)):
			sys.stderr.write("Cannot find %s. Need to run prep.py first.\n" % cr)
			sys.stderr.flush()
			sys.exit(1)
			break

	logger.info("Using existing input binaries...")
	expr = np.load("%s/expr.npy" % args.output)
	Xcen = np.load("%s/Xcen.npy" % args.output)
	genes = np.load("%s/genes.npy" % args.output)
	t_matrix = np.load("%s/t_matrix_%s_%.2f.npy" % (args.output, args.matrix_type, args.rbp_p))	
	ncell = Xcen.shape[0] 

	outdir = "%s/result_fast_sim_5000_%.2f_%.3f" % (args.output, rbp_p, examine_top)
	if matrix_type=="dissim":
		outdir = "%s/result_fast_5000_%.2f_%.3f" % (args.output, rbp_p, examine_top)
	if not os.path.isdir(outdir):
		os.mkdir(outdir)
	if not os.path.isfile("%s/do_gpd.R" % outdir):
		sys.stderr.write("do_gpd.R does not exist.\n")
		sys.stderr.flush()
		sys.exit(1)

	ncell = Xcen.shape[0]
	num_gt = int(ncell*examine_top)
	target = num_gt
	res = spatial_genes.random_pattern(matrix=t_matrix, matrix_type=matrix_type, num_cell = ncell, sizes=[target], trials_per_gene=5000, run_gpd=True, outdir=outdir, logger=logger)

def main():
	parser = argparse.ArgumentParser(description="evaluate.exact.one.2b.py: calculate silhouette score for randomly distributed spatial patterns", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("-r", "--rbp-p", dest="rbp_p", type=float, default=0.95, help="p parameter of RBP")
	parser.add_argument("-e", "--examine-top", dest="examine_top", type=float, default=0.05, help="top proportion of cells per gene to be 1's (expressed)")
	parser.add_argument("-m", "--matrix-type", dest="matrix_type", type=str, choices=["sim", "dissim"], help="whether to calculate similarity matrix or dissimilarity matrix", default="dissim")
	parser.add_argument("-o", "--output-dir", dest="output", type=str, default=".", help="output directory")
	parser.add_argument("-v", "--verbose", dest="verbose", action="store_true", help="print verbose messages to console")

	args = parser.parse_args()
	do_one(args)

if __name__=="__main__":
	main()

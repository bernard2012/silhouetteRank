import math
import sys
import os
import re
import scipy
import scipy.stats
import numpy as np
from operator import itemgetter
import silhouetteRank
import silhouetteRank.prep as prep
import silhouetteRank.evaluate_exact_one_2b as evaluate_exact_one_2b
import silhouetteRank.use_previous_cluster as use_previous_cluster
import silhouetteRank.combine as combine
import logging
import argparse
import subprocess

def read(n):
	f = open(n)
	by_gene = {}
	for l in f:
		l = l.rstrip("\n")
		ll = l.split()
		gene = ll[0]
		pval = float(ll[-2])
		by_gene[gene] = pval
	f.close()
	return by_gene

def do_one(args):
	result = subprocess.call("Rscript --version 2> /dev/null", shell=True)
	if result==127:
		sys.stderr.write("Rscript is not found\n")
		sys.stderr.flush()
		sys.exit(1)	

	check_required = ["expr.npy", "Xcen.npy", "genes.npy"]
	for cr in check_required:
		if not os.path.isfile("%s/%s" % (args.input, cr)):
			sys.stderr.write("%s file does not exist\n" % cr)
			sys.stderr.flush()
			sys.exit(1)

	for rbp_p in args.rbp_ps:
		for examine_top in args.examine_tops:
			random_dir = "%s/result_sim_5000_%.2f_%.3f" % (args.input, rbp_p, examine_top)
			score_file = "%s/silhouette.sim.exact.rbp.%.2f.top.%.3f.txt" % (args.input, rbp_p, examine_top)
			output_score_file = "%s/silhouette.sim.exact.rbp.%.2f.top.%.3f.pval.txt" % (args.input, rbp_p, examine_top)
			if args.matrix_type=="dissim":
				random_dir = "%s/result_5000_%.2f_%.3f" % (args.input, rbp_p, examine_top)
				score_file = "%s/silhouette.exact.rbp.%.2f.top.%.3f.txt" % (args.input, rbp_p, examine_top)
				output_score_file = "%s/silhouette.exact.rbp.%.2f.top.%.3f.pval.txt" % (args.input, rbp_p, examine_top)
			args1 = argparse.Namespace(expr="../expression.txt", centroid="../Xcen.good", examine_top=examine_top, input=score_file, input_random=random_dir, output=output_score_file, outdir=args.input, query_sizes=args.query_sizes, overwrite_input_bin=False, verbose=verbose, log_file="master.pvalue.log")
			use_previous_cluster.do_one(args1)

	combined_file = "%s/silhouette.overall.pval.txt" % args.input
	if args.matrix_type=="sim":
		combined_file = "%s/silhouette.sim.overall.pval.txt" % args.input
	args1 = argparse.Namespace(rbp_ps=args.rbp_ps, examine_tops=args.examine_tops, matrix_type=args.matrix_type, input=args.input, output=combined_file)
	combine.do_one(args1)

if __name__=="__main__":
	parser = argparse.ArgumentParser(description="combine.py: combine spatial scores across parameters", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("-r", "--rbp-ps", dest="rbp_ps", nargs="+", type=float, default=[0.95, 0.99], help="p parameter of RBP")
	parser.add_argument("-e", "--examine-tops", dest="examine_tops", nargs="+", type=float, default=[0.005, 0.010, 0.050, 0.100, 0.300], help="top proportion of cells per gene to be 1's (expressed)")
	parser.add_argument("-m", "--matrix-type", dest="matrix_type", type=str, choices=["sim", "dissim"], help="whether to calculate similarity matrix or dissimilarity matrix", default="dissim")
	parser.add_argument("-i", "--input-dir", dest="input", type=str, default=".", help="input directory containing individual spatial score rankings (to be aggregated)")
	#parser.add_argument("-o", "--output", dest="output", type=str, required=True, help="output file name")
	args = parser.parse_args()
	do_one(args)

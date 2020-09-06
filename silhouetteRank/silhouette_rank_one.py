import sys
import os
import re
import numpy as np
import subprocess
import math
import scipy
import silhouetteRank.spatial_genes as spatial_genes
from shutil import copyfile
from operator import itemgetter
from scipy.spatial.distance import squareform, pdist
from scipy.stats import percentileofscore
from sklearn.metrics import roc_auc_score
import argparse
import silhouetteRank
import silhouetteRank.prep as prep
import silhouetteRank.evaluate_exact_one_2b as evaluate_exact_one_2b
import silhouetteRank.use_previous_cluster as use_previous_cluster
import silhouetteRank.combine as combine
import logging

def main():
	parser = argparse.ArgumentParser(description="silhouette_rank_one.py: calculate silhouette score for randomly distributed spatial patterns", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("-x", "--file-expr", dest="expr", type=str, required=True, help="expression matrix. Will use input binary expr.npy (if exists) to speed up reading.")
	parser.add_argument("-c", "--file-centroid", dest="centroid", type=str, required=True, help="cell coordinate. Will use input binary Xcen.npy (if exists) to speed up reading.")
	parser.add_argument("-w", "--overwrite-input-binary", dest="overwrite_input_bin", action="store_true", help="overwrite input binaries")
	parser.add_argument("-r", "--rbp-ps", dest="rbp_ps", nargs="+", type=float, default=[0.95, 0.99], help="p parameter of RBP")
	parser.add_argument("-e", "--examine-tops", dest="examine_tops", nargs="+", type=float, default=[0.005, 0.010, 0.050, 0.100, 0.300], help="top proportion of cells per gene to be 1's (expressed)")
	parser.add_argument("-m", "--matrix-type", dest="matrix_type", type=str, choices=["sim", "dissim"], help="whether to calculate similarity matrix or dissimilarity matrix", default="dissim")
	parser.add_argument("-p", "--cores", dest="num_core", type=int, default=4, help="number of cores")
	parser.add_argument("-a", "--parallel-path", dest="parallel_path", type=str, default="/usr/bin", help="parallel binary path")
	parser.add_argument("-o", "--output-dir", dest="output", type=str, default=".", help="output directory")
	parser.add_argument("-q", "--query-sizes", dest="query_sizes", type=int, default=10, help="query sizes (advanced user setting)")
	parser.add_argument("-v", "--verbose", dest="verbose", action="store_true", help="print verbose messages to console")
	args = parser.parse_args()

	if not os.path.isdir(args.output):
		os.mkdir(args.output)

	logdir = "%s/logs" % args.output
	if not os.path.isdir(logdir):
		os.mkdir(logdir)

	verbose = args.verbose
	log_file = "%s/master.log" % args.output
	logger = logging.getLogger("master")
	logger.setLevel(logging.DEBUG)

	if not logger.hasHandlers():
		handler = logging.FileHandler(log_file)
		handler.setFormatter(logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s"))
		logger.addHandler(handler)
		if verbose:
			logger.addHandler(logging.StreamHandler())

	args1 = argparse.Namespace(expr=args.expr, centroid=args.centroid, rbp_ps=args.rbp_ps, examine_tops=args.examine_tops, matrix_type=args.matrix_type, output=args.output, query_sizes=args.query_sizes, overwrite_input_bin=args.overwrite_input_bin, verbose=verbose, log_file="master.prep.log")
	prep.do_one(args1)


	fw = open("%s/args" % args.output, "w")
	for rbp_p in args.rbp_ps:
		for examine_top in args.examine_tops:
			freq_file = "%s/result_5000_%.2f_%.3f/gene.freq.good.txt" % (args.output, rbp_p, examine_top)
			if args.matrix_type=="sim":
				freq_file = "%s/result_sim_5000_%.2f_%.3f/gene.freq.good.txt" % (args.output, rbp_p, examine_top)
			uniq_freq = 0
			f = open(freq_file)
			for l in f:
				l = l.rstrip("\n")
				uniq_freq+=1
			f.close()
			num_query_sizes = args.query_sizes
			if uniq_freq<=num_query_sizes:
				num_query_sizes = uniq_freq
			for i in range(num_query_sizes):
				fw.write("%.2f\n" % rbp_p)
				fw.write("%.3f\n" % examine_top)
				fw.write("%d\n" % i)
	fw.close()

	fw = open("%s/args.basic" % args.output, "w")
	for rbp_p in args.rbp_ps:
		for examine_top in args.examine_tops:
			fw.write("%.2f\n" % rbp_p)
			fw.write("%.3f\n" % examine_top)
	fw.close()

	bin_path = os.path.dirname(silhouetteRank.__file__)
	for i in range(4):
		bin_path = os.path.dirname(bin_path)
	bin_path = os.path.join(bin_path, "bin")	

	logger.info("Start calculating silhouette rank, saving logs to log directory (check progress here)...")
	cmd = "cat '%s'/args.basic | '%s'/parallel --jobs %d --max-args=2 \\''%s'\\'''/silhouette_rank_main -x \\''%s'\\''' -c \\''%s'\\''' -r {1} -e {2} -m %s -o \\''%s'\\'''" % (args.output, args.parallel_path, args.num_core, bin_path, args.expr, args.centroid, args.matrix_type, args.output)
	os.system(cmd)

	logger.info("Start randomization, saving logs to log directory (check progress here)...")
	cmd="cat '%s'/args | '%s'/parallel --jobs %d --max-args=3 \\''%s'\\'''/silhouette_rank_random -r {1} -e {2} -m %s -o \\''%s'\\''' -q {3}" % (args.output, args.parallel_path, args.num_core, bin_path, args.matrix_type, args.output)
	os.system(cmd)

	logger.info("Start computing P-values...")
	for rbp_p in args.rbp_ps:
		for examine_top in args.examine_tops:

			random_dir = "%s/result_sim_5000_%.2f_%.3f" % (args.output, rbp_p, examine_top)
			score_file = "%s/silhouette.sim.exact.rbp.%.2f.top.%.3f.txt" % (args.output, rbp_p, examine_top)
			output_score_file = "%s/silhouette.sim.exact.rbp.%.2f.top.%.3f.pval.txt" % (args.output, rbp_p, examine_top)
			if args.matrix_type=="dissim":
				random_dir = "%s/result_5000_%.2f_%.3f" % (args.output, rbp_p, examine_top)
				score_file = "%s/silhouette.exact.rbp.%.2f.top.%.3f.txt" % (args.output, rbp_p, examine_top)
				output_score_file = "%s/silhouette.exact.rbp.%.2f.top.%.3f.pval.txt" % (args.output, rbp_p, examine_top)
				

			args1 = argparse.Namespace(expr=args.expr, centroid=args.centroid, examine_top=examine_top, input=score_file, input_random=random_dir, output=output_score_file, outdir=args.output, query_sizes=args.query_sizes, overwrite_input_bin=args.overwrite_input_bin, verbose=verbose, log_file="master.pvalue.log")
			use_previous_cluster.do_one(args1)

	combined_file = "%s/silhouette.overall.pval.txt" % args.output
	if args.matrix_type=="sim":
		combined_file = "%s/silhouette.sim.overall.pval.txt" % args.output
	args1 = argparse.Namespace(rbp_ps=args.rbp_ps, examine_tops=args.examine_tops, matrix_type=args.matrix_type, input=args.output, output=combined_file)
	combine.do_one(args1)
	
if __name__=="__main__":
	main()

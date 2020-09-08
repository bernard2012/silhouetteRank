import math
import sys
import os
import re
import scipy
import scipy.stats
import numpy as np
from operator import itemgetter
import silhouetteRank
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

	outdir=args.input
	log_file = "%s/%s" % (outdir, args.log_file)
	logger = logging.getLogger("combine_fast")
	logger.setLevel(logging.DEBUG)

	if not logger.hasHandlers():
		handler = logging.FileHandler(log_file)
		handler.setFormatter(logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s"))
		logger.addHandler(handler)
		if args.verbose:
			logger.addHandler(logging.StreamHandler())

	logger.info("Entering %s..." % os.path.basename(args.input))

	logger.info("Using existing input binaries...")
	expr = np.load("%s/expr.npy" % args.input)
	Xcen = np.load("%s/Xcen.npy" % args.input)
	genes = np.load("%s/genes.npy" % args.input)

	ncell = Xcen.shape[0]
	param = {}
	full_list = {}
	for examine_top in args.examine_tops:
		for rbp_p in args.rbp_ps:
			target = int(ncell * examine_top)
			dname = "%s/result_fast_sim_5000_%.2f_%.3f" % (args.input, rbp_p, examine_top)
			if args.matrix_type=="dissim":
				dname = "%s/result_fast_5000_%.2f_%.3f" % (args.input, rbp_p, examine_top)
			f = open("%s/par.%d" % (dname, target))
			n_scale = float(f.readline().rstrip("\n").split("\t")[1])
			n_shape = float(f.readline().rstrip("\n").split("\t")[1])
			f.close()
			scores = []
			f = open("%s/%d" % (dname, target))
			for l in f:
				l = l.rstrip("\n")
				scores.append(float(l))
			f.close()
			scores = np.sort(np.array(scores))
			param[(rbp_p, examine_top, target)] = (n_scale, n_shape, scores[-250])
			full_list[(rbp_p, examine_top, target)] = scores		

	for examine_top in args.examine_tops:
		for rbp_p in args.rbp_ps:	
			for trial in range(args.num_trials):
				fname = "%s/silhouette.sim.fast.rbp.%.2f.top.%.3f.%d.txt" % (args.input, rbp_p, examine_top, trial)
				if args.matrix_type=="dissim":
					fname = "%s/silhouette.fast.rbp.%.2f.top.%.3f.%d.txt" % (args.input, rbp_p, examine_top, trial)	
				target = int(ncell * examine_top)
				f = open(fname)
				entries = []
				by_size = {}
				ids = {}
				ind = 0
				for l in f:
					l = l.rstrip("\n")
					ll = l.split("\t")
					g = ll[1]
					sc = float(ll[2])
					n_scale, n_shape, n_exceed = param[(rbp_p, examine_top, target)]
					if sc>n_exceed:
						xx = math.pow(1.0 - n_shape * (sc - n_exceed) / n_scale, 1.0 / n_shape)
						P = 0.05 * xx
					else:
						P = 0.01 * (100 - scipy.stats.percentileofscore(full_list[(rbp_p, examine_top, target)], sc))
					entries.append([g, target, target, sc, P])
					by_size.setdefault(target, [])
					by_size[target].append(P)
					ids.setdefault(target, [])
					ids[target].append(ind)
					ind+=1
				f.close()

				for targ in ids:
					fw = open("/tmp/1", "w")
					for i in by_size[targ]:
						fw.write(str(i) + "\n")
					fw.close()
					os.system("Rscript %s/qval.R /tmp/1 /tmp/1.qval" % os.path.dirname(silhouetteRank.__file__))
					s_scores = []
					f = open("/tmp/1.qval")
					for l in f:
						l = l.rstrip("\n")
						s_scores.append(float(l))
					f.close()
					for i1, i2 in zip(ids[targ], s_scores):
						entries[i1].append(i2)
				o_name = "%s/silhouette.sim.fast.rbp.%.2f.top.%.3f.%d.pval.txt" % (args.input, rbp_p, examine_top, trial)
				if args.matrix_type=="dissim":
					o_name = "%s/silhouette.fast.rbp.%.2f.top.%.3f.%d.pval.txt" % (args.input, rbp_p, examine_top, trial)
				fw = open(o_name, "w")
				for i1, i2, i3, i4, i5, i6 in entries:
					fw.write("%s %s %s %s %s %s\n" % (str(i1), str(i2), str(i3), str(i4), str(i5), str(i6)))
				fw.close()
				

	by_gene = {}
	for examine_top in args.examine_tops:
		for rbp_p in args.rbp_ps:
			for trial in range(args.num_trials):
				fname = "%s/silhouette.sim.fast.rbp.%.2f.top.%.3f.%d.pval.txt" % (args.input, rbp_p, examine_top, trial)	
				if args.matrix_type=="dissim":
					fname = "%s/silhouette.fast.rbp.%.2f.top.%.3f.%d.pval.txt" % (args.input, rbp_p, examine_top, trial)	
				by_gene[(examine_top, rbp_p, trial)] = read(fname)

	all_genes = list(by_gene[(args.examine_tops[0], args.rbp_ps[0], 0)].keys())
	score = {}
	pval = {}
	for g in all_genes:
		score[g] = 0
		tot_test = 0
		for i in args.examine_tops:
			for j in args.rbp_ps:
				for k in range(args.num_trials):
					score[g] += math.log(by_gene[(i, j, k)][g])
					tot_test+=1
		score[g] *= -2.0
		pval[g] = np.exp(scipy.stats.chi2.logsf(score[g], tot_test*2))

	score_it = list(score.items())
	score_it.sort(key=itemgetter(1), reverse=True)
	fw = open("/tmp/1.pval", "w")
	for i,j in score_it:
		fw.write(str(pval[i]) + "\n")
	fw.close()

	os.system("Rscript %s/qval.R /tmp/1.pval /tmp/1.qval" % os.path.dirname(silhouetteRank.__file__))
	f = open("/tmp/1.qval")
	q_score = []
	for l in f:
		l = l.rstrip("\n")
		q_score.append(float(l))
	f.close()

	fw = open(args.output, "w")
	for (i,j),k in zip(score_it, q_score):
		fw.write("%s %s %s %s\n" % (str(i), str(j), str(pval[i]), str(k)))
	fw.close()

if __name__=="__main__":
	parser = argparse.ArgumentParser(description="combine.py: combine spatial scores across parameters", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("-i", "--input", dest="input", type=str, required=True, help="input directory (should be whatever that contains result_5000)")
	#parser.add_argument("-j", "--input-random", dest="input_random", type=str, required=True, help="input random silhouette scores (for random patterns)")
	parser.add_argument("-o", "--output", dest="output", type=str, required=True, help="output file name")
	#parser.add_argument("-u", "--output-dir", dest="outdir", type=str, help="output directory containing binary files", default=".")
	parser.add_argument("-l", "--log-filename", dest="log_file", type=str, default="master.combine.fast.log", help="log file name (no path), will be generated in same directory as --output")
	parser.add_argument("-v", "--verbose", dest="verbose", action="store_true", help="print verbose messages to console")

	parser.add_argument("-r", "--rbp-ps", dest="rbp_ps", nargs="+", type=float, default=[0.95, 0.99], help="p parameter of RBP")
	parser.add_argument("-e", "--examine-tops", dest="examine_tops", nargs="+", type=float, default=[0.005, 0.010, 0.050, 0.100, 0.300], help="top proportion of cells per gene to be 1's (expressed)")
	parser.add_argument("-m", "--matrix-type", dest="matrix_type", type=str, choices=["sim", "dissim"], help="whether to calculate similarity matrix or dissimilarity matrix", default="dissim")
	parser.add_argument("-t", "--num-trials", dest="num_trials", type=int, default=1, help="number of trials")
	#parser.add_argument("-i", "--input-dir", dest="input", type=str, default=".", help="input directory containing individual spatial score rankings (to be aggregated)")
	#parser.add_argument("-o", "--output", dest="output", type=str, required=True, help="output file name")
	args = parser.parse_args()
	do_one(args)

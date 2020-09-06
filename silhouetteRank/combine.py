import math
import sys
import os
import re
import scipy
import scipy.stats
import numpy as np
from operator import itemgetter
import silhouetteRank

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
	by_gene = {}
	for examine_top in args.examine_tops:
		for rbp in args.rbp_ps:	
			fname = "%s/silhouette.sim.exact.rbp.%.2f.top.%.3f.pval.txt" % (args.input, rbp, examine_top)	
			if args.matrix_type=="dissim":
				fname = "%s/silhouette.exact.rbp.%.2f.top.%.3f.pval.txt" % (args.input, rbp, examine_top)	
			by_gene[(examine_top, rbp)] = read(fname)
	all_genes = list(by_gene[(args.examine_tops[0], args.rbp_ps[0])].keys())
	score = {}
	pval = {}
	for g in all_genes:
		score[g] = 0
		tot_test = 0
		for i in args.examine_tops:
			for j in args.rbp_ps:
				score[g] += math.log(by_gene[(i, j)][g])
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
	parser.add_argument("-r", "--rbp-ps", dest="rbp_ps", nargs="+", type=float, default=[0.95, 0.99], help="p parameter of RBP")
	parser.add_argument("-e", "--examine-tops", dest="examine_tops", nargs="+", type=float, default=[0.005, 0.010, 0.050, 0.100, 0.300], help="top proportion of cells per gene to be 1's (expressed)")
	parser.add_argument("-m", "--matrix-type", dest="matrix_type", type=str, choices=["sim", "dissim"], help="whether to calculate similarity matrix or dissimilarity matrix", default="dissim")
	parser.add_argument("-i", "--input-dir", dest="input", type=str, default=".", help="input directory containing individual spatial score rankings (to be aggregated)")
	parser.add_argument("-o", "--output", dest="output", type=str, required=True, help="output file name")
	args = parser.parse_args()
	do_one(args)

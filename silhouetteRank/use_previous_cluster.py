import sys
import os
import re
import subprocess
import numpy as np
import math
import scipy
from scipy import stats
import argparse
import logging

def read_freq(n):
	freq = []
	f = open(n)
	for l in f:
		l = l.rstrip("\n")
		ll = l.split(" ")
		for i in range(int(ll[0])):
			freq.append(int(ll[1]))
	f.close()
	freq = np.array(freq)
	return freq

def get_pattern_size(expr=None, Xcen=None, genes=None, examine_top=0.05):
	num_cell = Xcen.shape[0]
	ncell = num_cell
	#examine_top = 0.05
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
		pattern_size[g] = t_size
	return pattern_size

def read(f_expr="expression.txt", f_Xcen="Xcen.good", logger=logging):
	#f_expr = "expression.txt"
	#f_Xcen = "Xcen.good"

	#sys.stderr.write("Reading gene expression...\n")
	logger.info("Reading gene expression..")
	#sys.stderr.flush()

	f = open(f_expr)
	h = f.readline().rstrip("\n").split("\t")
	num_cell = len(h[1:])
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
	return expr, Xcen, genes

def do_one(args):
	result = subprocess.call("Rscript --version 2> /dev/null", shell=True)
	if result==127:
		sys.stderr.write("Rscript is not found\n")
		sys.stderr.flush()
		sys.exit(1)	

	outdir=args.input_random

	log_file = "%s/%s" % (os.path.dirname(args.output), args.log_file)
	logger = logging.getLogger("combine")
	logger.setLevel(logging.DEBUG)

	if not logger.hasHandlers():
		handler = logging.FileHandler(log_file)
		handler.setFormatter(logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s"))
		logger.addHandler(handler)
		if args.verbose:
			logger.addHandler(logging.StreamHandler())

	logger.info("Entering %s..." % os.path.basename(args.input))

	uniq_freq = 0
	f = open("%s/gene.freq.good.txt" % outdir)
	for l in f:
		l = l.rstrip("\n")
		uniq_freq+=1
	f.close()
	
	cc = []
	param = {}
	num_query_sizes = args.query_sizes
	if uniq_freq<=num_query_sizes:
		#sys.stderr.write("Setting query size %d instead of %d\n" % (uniq_freq, num_query_sizes))
		#sys.stderr.flush()
		logger.info("Setting query size %d instead of %d" % (uniq_freq, num_query_sizes))
		num_query_sizes = uniq_freq

	alist = ["%d" % x for x in list(range(num_query_sizes))]
	mapto = {}	
	full_list = {}
	for i in alist:
		f = open("%s/good/%s" % (outdir,i))
		line = f.readline().rstrip("\n").split(" ")
		target = int(line[1])
		for ix in line[2].split(";"):
			mapto[int(ix)] = target
		cc.append(np.array([target]))
		f.close()
		f = open("%s/par.%d" % (outdir,target))
		n_scale = float(f.readline().rstrip("\n").split("\t")[1])
		n_shape = float(f.readline().rstrip("\n").split("\t")[1])
		f.close()
		scores = []
		f = open("%s/%d" % (outdir, target))
		for l in f:
			l = l.rstrip("\n")
			scores.append(float(l))
		f.close()
		scores = np.sort(np.array(scores))
		param[target] = (n_scale, n_shape, scores[-250])
		full_list[target] = scores		
		
	cc = np.array(cc)
	t_overwrite = args.overwrite_input_bin
	#if there is no binary input file, create one
	check_required = ["expr.npy", "Xcen.npy", "genes.npy"]
	for cr in check_required:
		if not os.path.isfile("%s/%s" % (args.outdir, cr)):
			t_overwrite = True
			break

	expr, Xcen, genes, pattern_size = None, None, None, None
	if t_overwrite:
		expr, Xcen, genes = read(f_expr=args.expr, f_Xcen=args.centroid)
		pattern_size = get_pattern_size(expr=expr, Xcen=Xcen, genes=genes, examine_top=args.examine_top)
		np.save("%s/expr.npy" % args.outdir, expr)
		np.save("%s/Xcen.npy" % args.outdir, Xcen)
		np.save("%s/genes.npy" % args.outdir, genes)
	else:
		#sys.stderr.write("Using existing input binaries...\n")
		logger.info("Using existing input binaries...")
		sys.stderr.flush()
		expr = np.load("%s/expr.npy" % args.outdir)
		Xcen = np.load("%s/Xcen.npy" % args.outdir)
		genes = np.load("%s/genes.npy" % args.outdir)
		pattern_size = get_pattern_size(expr=expr, Xcen=Xcen, genes=genes, examine_top=args.examine_top)

	p_size = [pattern_size[g] for g in genes]
	pred = [mapto[p] for p in p_size]
	t_size = {}
	t_detect = {}
	for g, p1, p2 in zip(genes, p_size, pred):
		t_size[g] = p2
		t_detect[g] = p1

	f = open(args.input)
	entries = []
	by_size = {}
	ids = {}
	ind = 0

	for l in f:
		l = l.rstrip("\n")
		ll = l.split("\t")
		g = ll[1]
		sc = float(ll[2])
		n_scale, n_shape, n_exceed = param[t_size[g]]
		if sc>n_exceed:
			xx = math.pow(1.0 - n_shape * (sc - n_exceed) / n_scale, 1.0 / n_shape)
			P = 0.05 * xx
			#print(g, t_size[g], t_detect[g], sc, P)
		else:
			P = 0.01 * (100 - stats.percentileofscore(full_list[t_size[g]], sc))
			#print(g, t_size[g], t_detect[g], sc, P)
		entries.append([g, t_size[g], t_detect[g], sc, P])
		targ = t_size[g]
		by_size.setdefault(targ, [])
		by_size[targ].append(P)
		ids.setdefault(targ, [])
		ids[targ].append(ind)
		ind+=1
	f.close()

	for targ in ids:
		if len(by_size[targ])==1:
			s_scores = []
			s_scores.append(by_size[targ][0])
			for i1,i2 in zip(ids[targ], s_scores):
				entries[i1].append(i2)
			continue
		fw = open("/tmp/1", "w")
		for i in by_size[targ]:
			fw.write(str(i) + "\n")
		fw.close()
		os.system("cd '%s' && Rscript qval.R /tmp/1 /tmp/1.qval" % outdir)
		s_scores = []
		f = open("/tmp/1.qval")
		for l in f:
			l = l.rstrip("\n")
			s_scores.append(float(l))
		f.close()
		for i1, i2 in zip(ids[targ], s_scores):
			entries[i1].append(i2)
	fw = open(args.output, "w")
	for i1, i2, i3, i4, i5, i6 in entries:
		fw.write("%s %s %s %s %s %s\n" % (str(i1), str(i2), str(i3), str(i4), str(i5), str(i6)))
	fw.close()

if __name__=="__main__":
	parser = argparse.ArgumentParser(description="use_previous_cluster.py: calculate P-value for silhouette score based on random pattern silhouette distribution", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("-x", "--file-expr", dest="expr", type=str, required=True, help="expression matrix. Will use input binary expr.npy (if exists) to speed up reading.")
	parser.add_argument("-c", "--file-centroid", dest="centroid", type=str, required=True, help="cell coordinate. Will use input binary Xcen.npy (if exists) to speed up reading.")
	parser.add_argument("-w", "--overwrite-input-binary", dest="overwrite_input_bin", action="store_true", help="overwrite input binaries")
	parser.add_argument("-i", "--input", dest="input", type=str, required=True, help="input silhouette scores")
	parser.add_argument("-j", "--input-random", dest="input_random", type=str, required=True, help="input random silhouette scores (for random patterns)")
	parser.add_argument("-o", "--output", dest="output", type=str, required=True, help="output file name")
	parser.add_argument("-u", "--output-dir", dest="outdir", type=str, help="output directory containing binary files", default=".")
	parser.add_argument("-l", "--log-filename", dest="log_file", type=str, default="master.combine.log", help="log file name (no path), will be generated in same directory as --output")
	parser.add_argument("-v", "--verbose", dest="verbose", action="store_true", help="print verbose messages to console")
	parser.add_argument("-e", "--examine-top", dest="examine_top", type=float, default=0.05, help="used by evaluate.2b.py to generate silhouette scores")
	parser.add_argument("-q", "--query-sizes", dest="query_sizes", type=int, default=10, help="query sizes parameter used by evaluate.2b.py (advanced setting)")
	
	args = parser.parse_args()
	do_one(args)

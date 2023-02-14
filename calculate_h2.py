import argparse
import hydra
import os
import sys
from typing import Iterable
import time
from pyspark.sql import SparkSession
from pyspark.sql.types import *
import pyspark.sql.functions as F
from pyspark.sql.window import Window
from pyspark import SparkConf
import pandas as pd
import numpy
import scipy.special as sps


def _read_ld_sumstats(args, log, fh, alleles=False, dropna=True):
    sumstats = _read_sumstats(args, log, fh, alleles=alleles, dropna=dropna)
    ref_ld = _read_ref_ld(args, log)
    n_annot = len(ref_ld.columns) - 1
    M_annot = _read_M(args, log, n_annot)
    M_annot, ref_ld, novar_cols = _check_variance(log, M_annot, ref_ld)
    w_ld = _read_w_ld(args, log)
    sumstats = _merge_and_log(ref_ld, sumstats, 'reference panel LD', log)
    sumstats = _merge_and_log(sumstats, w_ld, 'regression SNP LD', log)
    w_ld_cname = sumstats.columns[-1]
    ref_ld_cnames = ref_ld.columns[1:len(ref_ld.columns)]
    return M_annot, w_ld_cname, ref_ld_cnames, sumstats, novar_cols
    

def estimate_h2(args, log):
    '''Estimate h2 and partitioned h2.'''
    args = copy.deepcopy(args)

    if args.samp_prev is not None and args.pop_prev is not None:
        args.samp_prev, args.pop_prev = map(
            float, [args.samp_prev, args.pop_prev])

    if args.intercept_h2 is not None:
        args.intercept_h2 = float(args.intercept_h2)
    if args.no_intercept:
        args.intercept_h2 = 1
    M_annot, w_ld_cname, ref_ld_cnames, sumstats, novar_cols = _read_ld_sumstats(
        args, log, args.h2)
    ref_ld = np.array(sumstats[ref_ld_cnames])
    _check_ld_condnum(args, log, ref_ld_cnames)
    _warn_length(log, sumstats)
    n_snp = len(sumstats)
    n_blocks = min(n_snp, args.n_blocks)
    n_annot = len(ref_ld_cnames)
    chisq_max = args.chisq_max
    old_weights = False
    if n_annot == 1:
        if args.two_step is None and args.intercept_h2 is None:
            args.two_step = 30
    else:
        old_weights = True
        if args.chisq_max is None:
            chisq_max = max(0.001*sumstats.N.max(), 80)

    s = lambda x: np.array(x).reshape((n_snp, 1))
    chisq = s(sumstats.Z**2)
    if chisq_max is not None:
        ii = np.ravel(chisq < chisq_max)
        sumstats = sumstats.ix[ii, :]
        log.log('Removed {M} SNPs with chi^2 > {C} ({N} SNPs remain)'.format(
                C=chisq_max, N=np.sum(ii), M=n_snp-np.sum(ii)))
        n_snp = np.sum(ii)  # lambdas are late-binding, so this works
        ref_ld = np.array(sumstats[ref_ld_cnames])
        chisq = chisq[ii].reshape((n_snp, 1))

    if args.two_step is not None:
        log.log('Using two-step estimator with cutoff at {M}.'.format(M=args.two_step))

    hsqhat = reg.Hsq(chisq, ref_ld, s(sumstats[w_ld_cname]), s(sumstats.N),
                     M_annot, n_blocks=n_blocks, intercept=args.intercept_h2,
                     twostep=args.two_step, old_weights=old_weights)

    if args.print_cov:
        _print_cov(hsqhat, args.out + '.cov', log)
    if args.print_delete_vals:
        _print_delete_values(hsqhat, args.out + '.delete', log)
        _print_part_delete_values(hsqhat, args.out + '.part_delete', log)

    log.log(hsqhat.summary(ref_ld_cnames, P=args.samp_prev, K=args.pop_prev, overlap = args.overlap_annot))
    if args.overlap_annot:
        overlap_matrix, M_tot = _read_annot(args, log)

        # overlap_matrix = overlap_matrix[np.array(~novar_cols), np.array(~novar_cols)]#np.logical_not
        df_results = hsqhat._overlap_output(ref_ld_cnames, overlap_matrix, M_annot, M_tot, args.print_coefficients)
        df_results.to_csv(args.out+'.results', sep="\t", index=False)
        log.log('Results printed to '+args.out+'.results')

    return hsqhat



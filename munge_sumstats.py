import pandas as pd
from pyspark.sql import DataFrame, SparkSession
import pyspark.sql.functions as f
from pyspark.sql.types import *
from pyspark.sql.window import Window
import numpy as np
import argparse
import gzip 
#from scipy.stats import chi2

""" 
preparation for munge:
	Read in sumstats
	Add Rs_IDs to all SNPs
	retain needed columns
	pass to next section:

ldsc/munge_sumstats.py:
	Read in sumstats
	Remove duplicate rs_IDs
	Convert P-value to Z-score
		p_to_z(p_val, sample_size):
			return np.sqrt(chi2.isf(P, 1))
	Join to w_hm3.snplist
	Write to output .sumstats.gz

	Calculate mean CHISQ stat for QC
		Mean chi^2 < 1.02 = Warning (too small)
	log metrics: (median Chi^2 / 0.4549) ?	| Max chi^2 | N SNPs remaining.
 """

#	L=L.toPandas()

#	gwn=L.at[0,"study_id"]
#	L.to_csv(path2save+gwn+".csv",index=False)
#	L.coalesce(1).write.format

def parse_args():
	parser=argparse.ArgumentParser()
	#parser.add_argument("--input_sumstats", help="input sumstat parquet file", type=str)
	parser.add_argument("--index", help="variant index to join rsID")
	# Variant index: gs://genetics-portal-dev-data/22.09.0/outputs/lut/variant-index
	parser.add_argument("--hm3", help="HM3 SNPs to join")
	# w_hm3.snplist : gs://genetics-portal-dev-analysis/xg1/Configs/w_hm3.snplist
	parser.add_argument("--outdir", help="Output directory for sumstats")
	# output to munge : gs://genetics-portal-dev-analysis/xg1/rsid_sumstats
	args = parser.parse_args()
	return args

def lex_order(str1,str2):
	str1=str1.upper()
	str2=str2.upper()
	if str1<str2:
		out=str1+"_"+str2
	else:
		out=str2+"_"+str1
	return out	

def main():
	args=parse_args()
	spark = SparkSession.builder.getOrCreate()

	# list GWAS sumstats:
	gwas_list=!gsutil ls gs://genetics-portal-dev-sumstats/unfiltered/gwas

	variant_index=spark.read.parquet(args.index)

	lex_orderUDF = f.udf(lambda z1,z2: lex_order(z1,z2),f.StringType())

	variant_index=variant_index.withColumn("lexa1a2",lex_orderUDF(f.col("ref_allele"),f.col("alt_allele")))
	variant_index=variant_index.withColumn("snpid",f.concat_ws("_",f.col("chr_id"),f.col("position"),f.col("lexa1a2")))
	
	VI=variant_index.select(f.col("rs_id"),f.col("snpid"))
	
	HM3_SNPs=spark.read.options(header=True, sep="\t").csv(args.hm3)
	VI=HM3_SNPs.join(VI,["rs_id"]).distinct()

	for gw in gwas_list[1:100]:
		gwas=spark.read.parquet(gw)

		gwas=gwas.withColumn("lexa1a2",lex_orderUDF(f.col("ref"),f.col("alt")))
		gwas=gwas.withColumn("snpid",f.concat_ws("_",f.col("chrom"),f.col("pos"),f.col("lexa1a2")))

		gwas=gwas.join(VI, ["snpid"]).distinct().select(f.col('rs_id').alias('SNP'), 
														f.col('pval').alias('P'), 
														f.col('ref').alias('A1'), 
														f.col('alt').alias('A2'), 
														f.col('n_total').alias('N'),
														f.col('n_cases').alias('N_CASES'),
														f.col('beta').alias('B'),
														f.col('se').alias('SE'),
														f.col('eaf').alias('EAF')).filter(f.col('SNP').isNotNull())
		gwas=gwas.toPandas()
		gwn=gw.replace('gs://genetics-portal-dev-sumstats/unfiltered/gwas/','')
		gwn=gwn.replace('.parquet','')

		gwas.to_csv(args.outdir+"/"+gwn+".sumstats.gz", index=False, compression="gzip", sep="\t")

if __name__ == '__main__':
    
    main()
# python munge_sumstats.py --input_sumstats ${trait}.parquet --index  

# Add rsID
# remove duplicates
# Only keep hapmap3 snps

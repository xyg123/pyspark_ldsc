import pandas as pd
from pyspark.sql import DataFrame, SparkSession
import pyspark.sql.functions as f
from pyspark.sql.types import *
from pyspark.sql.window import Window
import numpy as np
import argparse

spark = SparkSession.builder.getOrCreate()

gwas_list=!gsutil ls gs://genetics-portal-dev-sumstats/unfiltered/gwas

gene_index=spark.read.parquet("gs://genetics-portal-dev-data/22.09.0/outputs/lut/genes-index")
#gene_index=gene_index.toPandas()

variant_index=spark.read.parquet("gs://genetics-portal-dev-data/22.09.0/outputs/lut/variant-index")

HM3_SNPs=spark.read.csv('gs://genetics-portal-dev-analysis/xg1/Configs/listHM3.txt')

thr=5e-8
window_sig=5e5
window_gene=0
path2save="gs://genetics-portal-dev-sumstats/tmp/yt4-filtered-sumstats/"

def lex_order(str1,str2):
	str1=str1.upper()
	str2=str2.upper()
	if str1<str2:
 		out=str1+"_"+str2
	else:
 		out=str2+"_"+str1
	return out	

lex_orderUDF = f.udf(lambda z1,z2: lex_order(z1,z2),f.StringType())

variant_index=variant_index.withColumn("lexa1a2",lex_orderUDF(f.col("ref_allele"),f.col("alt_allele")))
variant_index=variant_index.withColumn("snpid",f.concat_ws("_",f.col("chr_id"),f.col("position"),f.col("lexa1a2")))
VI=variant_index.select(f.col("chr_id_b37"),f.col("position_b37"),f.col("rs_id"),f.col("snpid"))


HM3_SNPs=HM3_SNPs.withColumnRenamed("_c0", "rs_id")
HM3_SNPs=HM3_SNPs.join(VI,["rs_id"]).distinct()

gene_index=gene_index.withColumn("Window_start", f.col("start")-window_gene)
gene_index=gene_index.withColumn("Window_end", f.col("end")+window_gene)
gene_index=gene_index.select("chr", "Window_start", "Window_end")


gw=gwas_list[1]

cts=['rs_id','chrom','pos','ref','alt','beta','se','pval',
'n_total','n_cases','eaf','chr_id_b37','position_b37']

for i,gw in enumerate(gwas_list[1:1001]):
	print("N "+str(i)+": "+gw)
	gwas=spark.read.parquet(gw)
	gwas=gwas.withColumn("lexa1a2",lex_orderUDF(f.col("ref"),f.col("alt")))
	gwas=gwas.withColumn("snpid",f.concat_ws("_",f.col("chrom"),f.col("pos"),f.col("lexa1a2")))
	sig=gwas.filter(gwas.pval<=thr)
	L=gwas.filter("snpid is NULL")
	if sig.count()>0 :
    		sig=sig.withColumnRenamed("chrom", "chr")
			sig=sig.withColumn('Window_start', f.col("pos")-window_sig)
			sig=sig.withColumn('Window_end', f.col("pos")+window_sig)
			sig=sig.select("chr", "Window_start", "Window_end")
			l=gwas.alias('gwas').join(sig.alias('sig_gwas'), ( (f.col('gwas.chrom') == f.col('sig_gwas.chr')) & (f.col('gwas.pos') >= f.col('sig_gwas.Window_start')) & (f.col('gwas.pos') <= f.col('sig_gwas.Window_end')) ), how='inner')
			l=l.drop("chr", "Window_start", "Window_end")
			L=L.union(l).distinct()
	l=gwas.alias('gwas').join(gene_index.alias('gene_index'), ( (f.col('gwas.chrom') == f.col('gene_index.chr')) & (f.col('gwas.pos') >= f.col('gene_index.Window_start')) & (f.col('gwas.pos') <= f.col('gene_index.Window_end')) ), how='inner')
	l=l.drop("chr", "Window_start", "Window_end")
	L=L.union(l).distinct()
	l=gwas.join(HM3_SNPs.select("snpid"), on="snpid", how="inner").distinct()
	l=l.select(L.columns)
	L=L.union(l).distinct()
	L=L.join(VI,["snpid"]).distinct()
	L=L.select(cts)
	gwn=gw.replace('gs://genetics-portal-dev-sumstats/unfiltered/gwas/','')
	gwn=gwn.replace('.parquet','')
	L.write.csv(path2save+gwn)

#	L=L.toPandas()

#	gwn=L.at[0,"study_id"]
#	L.to_csv(path2save+gwn+".csv",index=False)
#	L.coalesce(1).write.format

def parse_args():
    	parser=argparse.ArgumentParser()
		parser.add_argument("--input_sumstats", help="input sumstat parquet file", type=str, required=True)
		parser.add_argument("--index", help="variant index to join rsID")
		parser.add_argument("--HM3_SNPs", help="Hapmap3 SNPs to include in the analysis")
		parser.add_argument("--gene_index")

def main():
    	args=parse_args()

print("Done!")

# python munge_sumstats.py --input_sumstats ${trait}.parquet --index  

# Add rsID
# remove duplicates
# Only keep hapmap3 snps

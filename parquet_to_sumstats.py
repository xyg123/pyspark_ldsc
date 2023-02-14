import pandas as pd
from pyspark.sql import DataFrame, SparkSession
import pyspark.sql.functions as f
import argparse


def parse_args():
    parser=argparse.ArgumentParser()
    parser.add_argument("--input_sumstats", help="input sumstat parquet file", type=str, required=True)
    parser.add_argument("--index", help="variant index to join rsID")
    parser.add_argument("--outdir", help="Output directory for sumstats")
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
    spark = SparkSession.builder.getOrCreate()

    args=parse_args()

    variant_index=spark.read.parquet(args.index)
    
    lex_orderUDF = f.udf(lambda z1,z2: lex_order(z1,z2),f.StringType())

    variant_index=variant_index.withColumn("lexa1a2",lex_orderUDF(f.col("ref_allele"),f.col("alt_allele")))
    variant_index=variant_index.withColumn("snpid",f.concat_ws("_",f.col("chr_id"),f.col("position"),f.col("lexa1a2")))
    VI=variant_index.select(f.col("chr_id_b37"),f.col("position_b37"),f.col("rs_id"),f.col("snpid"))


    gwas=spark.read.parquet(args.input_sumstats)
    cts=['rs_id','chrom','pos','ref','alt','beta','se','pval',
        'n_total','n_cases','eaf']

    gwas=gwas.withColumn("lexa1a2",lex_orderUDF(f.col("ref"),f.col("alt")))
    gwas=gwas.withColumn("snpid",f.concat_ws("_",f.col("chrom"),f.col("pos"),f.col("lexa1a2")))

    gwas=gwas.join(VI, ["snpid"]).distinct().select(cts)
    gwas=gwas.toPandas()
    gwas.to_csv(args.outdir+args.input_sumstats+".sumstats", index=False, sep="\t")

if __name__ == '__main__':
    
    main()

# python ~/pyspark_ldsc/parquet_to_sumstats.py --input_sumstats FINNGEN_R6_FG_CVD.parquet --index ~/ldsc/variant_index --outdir ../../input_sumstats/
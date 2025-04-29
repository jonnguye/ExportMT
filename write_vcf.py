import hail as hl
import argparse
import os

def write_vcf(inputs):
    #LOAD TABLES AND FIND SUBSET
    mt = hl.read_matrix_table(inputs['matrix_table'])
    samples_table = hl.import_table(inputs['samples_list'], no_header=True)
    #coerce first column to 's' in case the above fails
    print(samples_table.describe(), flush=True)
    samples_table = samples_table.rename({'0': 's'})
    
    mt = mt.filter_cols(hl.is_defined(samples_table[mt.s]))
    print(f"Filtering to {mt.count_cols()} samples")

    #SELECT WHICH CHROMOSOME TO FILTER BY

    if inputs['chr'] is None or inputs['chr'].upper() == 'ALL':
        print("NO CHR FILTER APPLIED")
        print(f"CHR value provided: {inputs['chr']}")
    else:
        print(f"Filtering on {inputs['chr']}")
        mt = mt.filter_rows( (mt.locus.contig == inputs['chr']) )
    
    #ONLY CONTAINS PASS IN FT
    #IF FT is not pass, set to 0,0
    mt = mt.annotate_entries(GT = hl.if_else(hl.is_defined(mt.FT) & (mt.FT == "PASS"),mt.GT,hl.call(0, 0)))

    #Save total pop data
    mt = mt.annotate_rows(
        total = mt.info.annotate(
            ALL_AF = hl.min(mt.info.AF),
            ALL_AN = mt.info.AN,
            ALL_AC = hl.min(mt.info.AC),
            ALL_p_value_hwe = mt.variant_qc.p_value_hwe,
            ALL_p_value_excess_het = mt.variant_qc.p_value_excess_het
        )
    )
    
    #OVERWRITE TOTAL POPULATION INFO WITH SUBPOPULATION INFO
    mt = mt.annotate_rows( info = hl.agg.call_stats(mt.GT, mt.alleles) )
    
    #95% of alleles called in the population
    mt = mt.filter_rows(mt.info.AN >= 0.95 * mt.count_cols() * 2)

    #recalculate per subpopulation
    mt = hl.variant_qc(mt)

    #only show minor allele data
    mt = mt.annotate_rows(
        info = mt.info.annotate(
            ALL_AF = mt.total.ALL_AF,
            ALL_AC = mt.total.ALL_AC,
            ALL_AN = mt.total.ALL_AN,
            ALL_p_value_hwe = mt.total.ALL_p_value_hwe,
            ALL_p_value_excess_het = mt.total.ALL_p_value_excess_het,
            AF = hl.min(mt.info.AF),
            AC = hl.min(mt.info.AC),
            AN = mt.info.AN,
            p_value_hwe = mt.variant_qc.p_value_hwe,
            p_value_excess_het = mt.variant_qc.p_value_excess_het
        )
    )

    #FILTER OUT MONOMORPHIC ALLELES
    threshold = 0.001
    mt = mt.filter_rows(mt.info.AF > threshold)
    mt = mt.filter_rows(mt.info.AF < 1 - threshold)

    #FILTER BY MIN AC
    mt = mt.filter_rows(mt.info.AC >= inputs['MinimumAC_inclusive'])
    
    hl.export_vcf(mt, f"{inputs['output_prefix']}.vcf.bgz")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix_table", required=True)
    parser.add_argument("--samples_list", required=True)
    parser.add_argument("--chr", required=True)
    parser.add_argument("--MinimumAC_inclusive", type=int, required=True)
    parser.add_argument("--output_prefix", required=True)

    args = parser.parse_args()

    inputs = {
        'matrix_table': args.matrix_table,
        'samples_list': args.samples_list,
        'chr': args.chr,
        'MinimumAC_inclusive': args.MinimumAC_inclusive,
        'output_prefix': args.output_prefix,
    }

    hl.init(
        app_name='hail_job',
        master='local[*]',
        tmp_dir='gs://fc-secure-b8771cfd-5455-4292-a720-8533eb501a93/hail-tmp/',
        spark_conf={
            'spark.executor.instances': '4',
            'spark.executor.cores': '8',
            'spark.executor.memory': '25g',
            'spark.driver.memory': '30g',
            'spark.local.dir': '/cromwell_root',
            'spark.sql.shuffle.partitions': '100',
            'spark.default.parallelism': '100',
            'spark.memory.fraction': '0.8',
            'spark.memory.storageFraction': '0.2',
        },
        default_reference='GRCh38'
    )
    
    print("Spark local directories:", os.getenv("SPARK_LOCAL_DIRS"), flush=True)
    print("Disk usage:", flush=True)
    os.system("df -h")
    
    write_vcf(inputs)

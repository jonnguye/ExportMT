version 1.0

workflow WriteVCFWorkflow {
    input {
        String matrix_table
        String ancestry_table
        String ancestry
        String chr
        Int MinimumAC_inclusive
        String output_path
        Int new_id_max_allele_len
        String output_prefix
        File genotype_rscript
    }

    call WriteVCFTask {
        input:
            matrix_table = matrix_table,
            ancestry_table = ancestry_table,
            ancestry = ancestry,
            chr = chr,
            MinimumAC_inclusive = MinimumAC_inclusive,
            output_path = output_path,
    }

    call plink2 {
        input:
            vcf_file = WriteVCFTask.output_vcf,
            output_prefix = output_prefix,
            new_id_max_allele_len = new_id_max_allele_len
    }

    call ComputeGenotypePCS {
        input:
            vcf_file = WriteVCFTask.output_vcf,
            output_prefix = output_prefix,
            genotype_rscript = genotype_rscript
    }

    output {
        File output_vcf = WriteVCFTask.output_vcf
        Array[File] plink_outputs = plink2.plink_outputs
        File genotype_pcs = ComputeGenotypePCS.output_tsv
    }
}

task WriteVCFTask {
    input {
        String matrix_table
        String ancestry_table
        String ancestry
        String chr
        Int MinimumAC_inclusive
        String output_path
    }

    command <<<
        set -e

        export SPARK_LOCAL_DIRS=/cromwell_root

        echo "Checking disk mounts and usage:"
        df -h
        echo "Checking Spark local directory:"
        echo $SPARK_LOCAL_DIRS
        echo "Checking /cromwell_root directory:"
        ls -lah /cromwell_root

        curl -O https://raw.githubusercontent.com/jonnguye/ExportMT/main/write_vcf.py

        python3 write_vcf.py \
            --matrix_table "~{matrix_table}" \
            --ancestry_table "~{ancestry_table}" \
            --ancestry "~{ancestry}" \
            --chr "~{chr}" \
            --MinimumAC_inclusive "~{MinimumAC_inclusive}" \
            --output_path "~{output_path}" \
    >>>

    runtime {
        docker: "quay.io/jonnguye/hail:latest"
        memory: "256G"
        cpu: 64
        disks: "local-disk 1000 SSD"
    }

    output {
        File output_vcf = "~{output_path}"
    }
}

task plink2 {
    input {
        File vcf_file
        String output_prefix
        Int new_id_max_allele_len
    }

    command <<<
        set -e

        mkdir -p plink_output

        plink2 --vcf "~{vcf_file}" \
        --make-pgen \
        --out plink_output/"~{output_prefix}" \
        --set-all-var-ids @:#\$r_\$a \
        --new-id-max-allele-len "~{new_id_max_allele_len}" \
        --output-chr chrM \
        --chr 1-22
    >>>

    runtime {
        docker: "quay.io/biocontainers/plink2:2.0.0a.6.9--h9948957_0"
        memory: "16G"
        cpu: 4
        disks: "local-disk 100 SSD"
    }

    output {
        Array[File] plink_outputs = glob("plink_output/*")
    }
}

task ComputeGenotypePCS {
    input {
        File vcf_file
        String output_prefix
        File genotype_rscript
    }

    command <<<
        set -e

        Rscript "~{genotype_rscript}" \
            --vcf_path "~{vcf_file}" \
            --prefix "~{output_prefix}"
        >>>
    
        runtime {
            docker: "quay.io/jonnguye/genotype_pcs:micromamba"
            memory: "8G"
            cpu: 2
        }
    
        output {
            File output_tsv = "~{output_prefix}_genetic_PCs.tsv"
        }
}

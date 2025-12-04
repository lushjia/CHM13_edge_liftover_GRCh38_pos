#!/usr/bin/env nextflow

/*
 * find edge corresponding CHM13 variants 
 */

nextflow.enable.dsl = 2

// 1. Define flexible parameters
params.edge_file = "data/chr{chr}.edge_list.txt"
params.pangenie_decom_bi_vcf = "data/chr{chr}.var.vcf"
params.liftover_vcf = "data/chr{chr}.liftover.vcf"

params.script_edge_graphvar = "scripts/edge_graphvar.py"
params.script_edge_b38pos = "scripts/edge_b38pos.py"
params.outdir = "results"

// 2. Process 1: find CHM13 edge corresponding CHM13 variants
process find_edge_var {
    tag "Chr ${chr}" // print processing chr in log 

    input:
    path py_script
    tuple val(chr), path(edge), path(pangenie_decom_vcf), path(liftover)

    output:
    // We pass the 'chr' variable along with the output file
    tuple val(chr), path("chr${chr}.edge_var.txt"), path("chr${chr}.var_index.txt"), path(liftover), emit: intermediate_data

    script:
    """
    python ${py_script} -edge ${edge} -var ${pangenie_decom_vcf} -o chr${chr}
    """
}

// 3. Process 2: find b38 position corresponding CHM13 edge
process find_edge_b38pos {
    tag "Chr ${chr}" // print processing chr in log 
    publishDir params.outdir, mode: 'copy'

    input:
    path py_script
    tuple val(chr), path(edge_var), path(var_index), path(liftover_vcf)

    output:
    tuple path("chr${chr}_edge_b38pos.txt"), path("chr${chr}_edge_b38pos_leftmost.txt")

    script:
    """
    python ${py_script} -edge_var ${edge_var} -var_index ${var_index} -liftover ${liftover_vcf} -o chr${chr}
    """
}

// 4. Workflow
workflow {
    chr_ch = channel.of(1..22)

    // creates a channel of tuples: 
    inputs_ch = chr_ch.map { chr_num ->
        def edge_file = file(params.edge_file.replace("{chr}", "${chr_num}"))
        def pangenie_decom_bi_vcf = file(params.pangenie_decom_bi_vcf.replace("{chr}", "${chr_num}"))
        def liftover_vcf = file(params.liftover_vcf.replace("{chr}", "${chr_num}"))
        // Return the tuple structure Step 1 expects
        return [chr_num, edge_file, pangenie_decom_bi_vcf, liftover_vcf]
    }

    // Notice inputs_ch contains the data tuples, so we pass it as the 2nd argument
    step1_results = find_edge_var(file(params.script_edge_graphvar), inputs_ch)

    find_edge_b38pos(file(params.script_edge_b38pos), step1_results.intermediate_data)
}


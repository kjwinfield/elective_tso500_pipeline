#!/usr/bin/env nextflow
params.outdir = "./results"
params.input = "*.fastq.gz"
params.bwaindex = "*"
params.fasta = "*.fa"
params.fasta_fai = "*.fa.fai"
params.bam = "*.bam"
params.bam_bai = "*.bam.bai"

process FASTQC {
    tag "$meta.id"
    // label 'process_medium'
    publishDir "$params.outdir/fastqc", mode:'copy'

    conda "bioconda::fastqc=0.11.9"
    container "quay.io/biocontainers/fastqc:0.11.9--0"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip") , emit: zip
    path  "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Make list of old name and new name pairs to use for renaming in the bash while loop
    def old_new_pairs = reads instanceof Path || reads.size() == 1 ? [[ reads, "${prefix}.${reads.extension}" ]] : reads.withIndex().collect { entry, index -> [ entry, "${prefix}_${index + 1}.${entry.extension}" ] }
    def rename_to = old_new_pairs*.join(' ').join(' ')
    def renamed_files = old_new_pairs.collect{ old_name, new_name -> new_name }.join(' ')
    """
    printf "%s %s\\n" $rename_to | while read old_name new_name; do
        [ -f "\${new_name}" ] || ln -s \$old_name \$new_name
    done
    fastqc $args --threads $task.cpus $renamed_files
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$( fastqc --version | sed -e "s/FastQC v//g" )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.html
    touch ${prefix}.zip
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$( fastqc --version | sed -e "s/FastQC v//g" )
    END_VERSIONS
    """
}

process BWAMEM2_INDEX {
    tag "$fasta"
    label 'process_single'

    conda "bioconda::bwa-mem2=2.2.1"
    container 'quay.io/biocontainers/bwa-mem2:2.2.1--he513fc3_0'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("bwamem2"), emit: index
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir bwamem2
    bwa-mem2 \\
        index \\
        $args \\
        $fasta -p bwamem2/${fasta}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwamem2: \$(echo \$(bwa-mem2 version 2>&1) | sed 's/.* //')
    END_VERSIONS
    """

    stub:
    """
    mkdir bwamem2
    touch bwamem2/${fasta}.0123
    touch bwamem2/${fasta}.ann
    touch bwamem2/${fasta}.pac
    touch bwamem2/${fasta}.amb
    touch bwamem2/${fasta}.bwt.2bit.64
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwamem2: \$(echo \$(bwa-mem2 version 2>&1) | sed 's/.* //')
    END_VERSIONS
    """
}

process BWAMEM2_MEM {
    tag "$meta.id"
    // label 'process_high'
    publishDir "$params.outdir/bwa-mem", mode:'copy'

    conda "bioconda::bwa-mem2=2.2.1 bioconda::samtools=1.16.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-e5d375990341c5aef3c9aff74f96f66f65375ef6:2cdf6bf1e92acbeb9b2834b1c58754167173a410-0' :
        'quay.io/biocontainers/mulled-v2-e5d375990341c5aef3c9aff74f96f66f65375ef6:2cdf6bf1e92acbeb9b2834b1c58754167173a410-0' }"

    input:
    tuple val(meta), path(reads)
    path(index)
    val   sort_bam

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def samtools_command = sort_bam ? 'sort' : 'view'
    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`
    bwa-mem2 \\
        mem \\
        $args \\
        -t $task.cpus \\
        \$INDEX \\
        $reads \\
        | samtools $samtools_command $args2 -@ $task.cpus -o ${prefix}.bam -
        | samtools index -
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwamem2: \$(echo \$(bwa-mem2 version 2>&1) | sed 's/.* //')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwamem2: \$(echo \$(bwa-mem2 version 2>&1) | sed 's/.* //')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}

process FREEBAYES {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::freebayes=1.3.6"
    container 'quay.io/biocontainers/freebayes:1.3.6--hbfe0e7f_2'

    input:
    tuple val(meta), path(input_1), path(input_1_index), path(input_2), path(input_2_index), path(target_bed)
    path fasta
    path fasta_fai
    // path samples
    // path populations
    // path cnv

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    path  "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input            = input_2        ? "${input_1} ${input_2}"        : "${input_1}"
    def targets_file     = target_bed     ? "--target ${target_bed}"       : ""
    // def samples_file     = samples        ? "--samples ${samples}"         : ""
    // def populations_file = populations    ? "--populations ${populations}" : ""
    // def cnv_file         = cnv            ? "--cnv-map ${cnv}"             : ""

    """
    freebayes \\
        -f $fasta \\
        $targets_file \\
      #  $samples_file \\
       # $populations_file \\
      #  $cnv_file \\
        $args \\
        $input > ${prefix}.vcf
    bgzip ${prefix}.vcf
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        freebayes: \$(echo \$(freebayes --version 2>&1) | sed 's/version:\s*v//g' )
    END_VERSIONS
    """
}

process ENSEMBLVEP_VEP {
    tag "$meta.id"
    // label 'process_medium'
    publishDir "$params.outdir/vep", mode:'copy'


    conda "bioconda::ensembl-vep=108.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ensembl-vep:108.2--pl5321h4a94de4_0' :
        'quay.io/biocontainers/ensembl-vep:108.2--pl5321h4a94de4_0' }"

    input:
    tuple val(meta), path(vcf), path(custom_extra_files)
    val   genome
    val   species
    val   cache_version
    path  cache
    tuple val(meta2), path(fasta)
    path  extra_files

    output:
    tuple val(meta), path("*.vcf.gz")  , optional:true, emit: vcf
    tuple val(meta), path("*.tab.gz")  , optional:true, emit: tab
    tuple val(meta), path("*.json.gz") , optional:true, emit: json
    path "*.summary.html"              , emit: report
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def file_extension = args.contains("--vcf") ? 'vcf' : args.contains("--json")? 'json' : args.contains("--tab")? 'tab' : 'vcf'
    def compress_cmd = args.contains("--compress_output") ? '' : '--compress_output bgzip'
    def prefix = task.ext.prefix ?: "${meta.id}"
    def dir_cache = cache ? "\${PWD}/${cache}" : "/.vep"
    def reference = fasta ? "--fasta $fasta" : ""
    """
    vep \\
        -i $vcf \\
        -o ${prefix}.${file_extension}.gz \\
        $args \\
        $compress_cmd \\
        $reference \\
        --assembly $genome \\
        --species $species \\
        --cache \\
        --cache_version $cache_version \\
        --dir_cache $dir_cache \\
        --fork $task.cpus \\
        --stats_file ${prefix}.summary.html \\
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ensemblvep: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : //;s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf.gz
    touch ${prefix}.tab.gz
    touch ${prefix}.json.gz
    touch ${prefix}.summary.html
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ensemblvep: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : //;s/ .*\$//')
    END_VERSIONS
    """
}

// define workflow
workflow {
    ch_fastq = Channel.fromFilePairs(params.input, size: -1)
        .map {
            meta, fastq ->
            def fmeta = [:]
            // Set meta.id
            fmeta.id = meta
            // Set meta.single_end
            if (fastq.size() == 1) {
                fmeta.single_end = true
            } else {
                fmeta.single_end = false
            }
            [ fmeta, fastq ]
	}
    ch_fastq.view()
   // FASTQC ( ch_fastq )
   // BWAMEM2_INDEX([ch_fastq, params.fasta])
   // BWAMEM2_MEM(ch_fastq, params.bwaindex, true)
     FREEBAYES([ch_fastq, params.bam, params.bam_bai], params.fasta, params.fasta_fai)
    // | ENSEMBLVEP_VEP
}
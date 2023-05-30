#!/usr/bin/env nextflow
params.outdir = "./results"
params.input = "*_R{1,2}.fastq.gz"
params.bwaindex = "*"
params.fasta = "*.fa"
params.fasta_fai = "*.fa.fai"
params.dict = "*.dict"
params.index = false

include { ALIGN_BAM       as ALIGN_RAW_BAM               } from './modules/alignbam/main'
include { ALIGN_BAM       as ALIGN_CONSENSUS_BAM         } from './modules/alignbam/main'



process FASTQC {
    tag "$meta.id"
    label 'process_medium'
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

process FGBIO_FASTQTOBAM {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::fgbio=2.0.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fgbio:2.0.2--hdfd78af_0' :
        'quay.io/biocontainers/fgbio:2.0.2--hdfd78af_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.bam") , emit: bam , optional: true
    tuple val(meta), path("*.cram"), emit: cram, optional: true
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sample_name = args.contains("--sample") ? "" : "--sample ${prefix}"
    def library_name = args.contains("--library") ? "" : "--library ${prefix}"
    def output = prefix =~ /\.(bam|cram)$/ ? prefix : "${prefix}.bam"
    """

    fgbio \\
        --tmp-dir=. \\
        FastqToBam \\
        ${args} \\
        --input ${reads} \\
        --output ${output} \\
        ${sample_name} \\
        ${library_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
    END_VERSIONS
    """
}

process SAMTOOLS_INDEX {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::samtools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.17--h00cdaf9_0' :
        'quay.io/biocontainers/samtools:1.17--h00cdaf9_0' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*.bai") , optional:true, emit: bai
    tuple val(meta), path("*.csi") , optional:true, emit: csi
    tuple val(meta), path("*.crai"), optional:true, emit: crai
    path  "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    samtools \\
        index \\
        -@ ${task.cpus-1} \\
        $args \\
        $input
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch ${input}.bai
    touch ${input}.crai
    touch ${input}.csi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}

process FREEBAYES {
    tag "$meta.id"
    label 'process_medium'
    publishDir "$params.outdir/freebayes", mode:'copy'

    conda "bioconda::freebayes=1.3.6"
    container 'quay.io/biocontainers/freebayes:1.3.6--hbfe0e7f_2'

    input:
    tuple val(meta), path(input_1), path(input_1_index), path(input_2), path(input_2_index), path(target_bed)
    path fasta
    path fasta_fai
    path samples
    path populations
    path cnv

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
    def samples_file     = samples        ? "--samples ${samples}"         : ""
    def populations_file = populations    ? "--populations ${populations}" : ""
    def cnv_file         = cnv            ? "--cnv-map ${cnv}"             : ""

    """
    freebayes \\
        -f $fasta \\
        $targets_file \\
        $samples_file \\
        $populations_file \\
        $cnv_file \\
        $args \\
        $input > ${prefix}.vcf
    bgzip ${prefix}.vcf
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        freebayes: \$(echo \$(freebayes --version 2>&1) | sed 's/version:\s*v//g' )
    END_VERSIONS
    """
}

process FGBIO_SETMATEINFORMATION{
    label 'process_medium'
    container "davelabhub/fgbio:2.0.2--hdfd78af_0"

    input:
    tuple val(meta), path(bamfile)

    output:
    tuple val(meta), path("*_w_mate_info.bam"), emit: bam

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    fgbio \\
        SetMateInformation \\
        -i $bamfile \\
        -o ${prefix}_w_mate_info.bam
        

    """
}

process FGBIO_GROUPREADSBYUMI {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::fgbio=2.0.2"
    container 'quay.io/biocontainers/fgbio:2.0.2--hdfd78af_0'

    input:
    tuple val(meta), path(taggedbam)
    val(strategy)

    output:
    tuple val(meta), path("*_umi-grouped.bam")  , emit: bam
    tuple val(meta), path("*_umi_histogram.txt"), emit: histogram
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """

    fgbio \\
        --tmp-dir=. \\
        GroupReadsByUmi \\
        -s $strategy \\
        $args \\
        -i $taggedbam \\
        -o ${prefix}_umi-grouped.bam \\
        -f ${prefix}_umi_histogram.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
    END_VERSIONS
    """
}

process FGBIO_CALLDUPLEXCONSENSUSREADS {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::fgbio=2.0.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fgbio:2.0.2--hdfd78af_0' :
        'quay.io/biocontainers/fgbio:2.0.2--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam)
    // please note:
    // --min-reads is a required argument with no default
    // --min-input-base-quality is a required argument with no default
    // make sure they are specified via ext.args in your config

    output:
    tuple val(meta), path("${prefix}.bam"), emit: bam
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_consensus"

    def mem_gb = 8
    if (!task.memory) {
        log.info '[fgbio CallDuplexConsensusReads] Available memory not known - defaulting to 8GB. Specify process memory requirements to change this.'
    } else if (mem_gb > task.memory.giga) {
        if (task.memory.giga < 2) {
            mem_gb = 1
        } else {
            mem_gb = task.memory.giga - 1
        }
    }

    """
    fgbio \\
        -Xmx${mem_gb}g \\
        --tmp-dir=. \\
        --async-io=true \\
        --compression=1 \\
        CallDuplexConsensusReads \\
        --input $bam \\
        --output ${prefix}.bam \\
        --threads ${task.cpus} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
    END_VERSIONS
    """
}

process ENSEMBLVEP_DOWNLOAD {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::ensembl-vep=108.2"
    container "quay.io/biocontainers/ensembl-vep:108.2--pl5321h4a94de4_0"

    input:
    tuple val(meta), val(assembly), val(species), val(cache_version)

    output:
    tuple val(meta), path("vep_cache"), emit: cache
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    vep_install \\
        --CACHEDIR vep_cache \\
        --SPECIES $species \\
        --ASSEMBLY $assembly \\
        --CACHE_VERSION $cache_version \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ensemblvep: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : //;s/ .*\$//')
    END_VERSIONS
    """

    stub:
    """
    mkdir vep_cache

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ensemblvep: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : //;s/ .*\$//')
    END_VERSIONS
    """
}

process ENSEMBLVEP_VEP {
    tag "$meta.id"
    label 'process_medium'
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
  //  def reference = fasta ? "--fasta $fasta" : ""
    """
    vep \\
        -i $vcf \\
        -o ${prefix}.${file_extension}.gz \\
        $args \\
        $compress_cmd \\
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

    // run fastqc to assess read quality
    FASTQC ( ch_fastq )
    
    //fix RX tags 
    FGBIO_FASTQTOBAM(ch_fastq)

    // initial alignment
    ALIGN_RAW_BAM(FGBIO_FASTQTOBAM.out.bam, params.fasta, params.dict, true)
    
    // handle UMIs and collapse consensus reads based on UMI tags and alignment
    FGBIO_SETMATEINFORMATION(ALIGN_RAW_BAM.out.bam)
    FGBIO_GROUPREADSBYUMI(FGBIO_SETMATEINFORMATION.out.bam, "Edit")
    FGBIO_CALLDUPLEXCONSENSUSREADS(FGBIO_GROUPREADSBYUMI.out.bam)
    
    // realign (using only consensus reads)
    ALIGN_CONSENSUS_BAM(FGBIO_CALLDUPLEXCONSENSUSREADS.out.bam, params.fasta, params.dict, false)
    
    // index bam to create bam.bai
    SAMTOOLS_INDEX(ALIGN_CONSENSUS_BAM.out.bam)

    // create channel with empty tuples for input to freebayes
    freebayes_ch = ALIGN_CONSENSUS_BAM.out.bam.join(SAMTOOLS_INDEX.out.bai)
    freebayes_ch.view()
    blank_ch = Channel.of([[],[],[]])
    blank_ch.view()
    freebayes_concat_ch = freebayes_ch.concat(blank_ch).toList()
    freebayes_concat_ch.view()
    freebayes_mapped = freebayes_concat_ch.map { it -> [it[0][0],it[0][1], it[0][2], it[1][0], it[1][1], it[1][2]]}
    freebayes_mapped.view()
   
    // run freebayes to create vcf
    FREEBAYES(freebayes_mapped, params.fasta, params.fasta_fai, [], [], [])

    // implement ENSEMBL VEP vcf annotation in future 
  // ENSEMBLVEP_DOWNLOAD([['id': 'test', single_end: false],"GRCh37","homo_sapiens","108"])
  //  ENSEMBLVEP_VEP([['id': 'test', single_end: false], '/home/dnanexus/test.vcf.gz', [],], "GRCh37", "homo_sapiens", "108", '/home/dnanexus/vep_cache')
}
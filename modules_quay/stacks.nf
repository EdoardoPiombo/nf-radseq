process process_radtags_pe {

  container "quay.io/biocontainers/stacks:2.65--hdcf5f25_0"
  executor params.executor
  publishDir "analysis/stacks"
  publishDir "report/logs/", pattern: "*.log"
  cpus 1

  input:
    tuple path(read1), path(read2), val(name)
    val args
    val meta

  output:
    tuple path("${name}_trimmed{.1,.rem.1,.2,.rem.2}.fq.gz"), val(name), emit: data
    path "${name}_process_radtags.log", emit: log
    val meta, emit: meta

  script:
    _args = args.join(" ")

    """
    process_radtags -1 ${read1} -2 ${read2} \
      -o . \
      --paired \
      ${_args} &> ${name}_process_radtags.log
    mv ${name}_trimmed_1.1.fq.gz ${name}_trimmed.1.fq.gz
    mv ${name}_trimmed_1.rem.1.fq.gz ${name}_trimmed.rem.1.fq.gz
    mv ${name}_trimmed_2.2.fq.gz ${name}_trimmed.2.fq.gz
    mv ${name}_trimmed_2.rem.2.fq.gz ${name}_trimmed.rem.2.fq.gz
    """
}

process ustacks {

  container "quay.io/biocontainers/stacks:2.65--hdcf5f25_0"
  executor params.executor
  cpus params.ustacks_cpus
  publishDir "analysis/stacks"
  publishDir "report/logs/", pattern: "*.log"

  input:
    tuple path(reads), val(name)
    val args
    val meta

  output:
    tuple path("_tmpout/*"), val(name), emit: data
    path "${name}_ustacks.log", emit: log
    val meta, emit: meta

  script:
    _args = args.join(" ")
    _id = meta.NumId
    """
    mkdir _tmpout
    ustacks -f ${reads[0]} \
      -i ${_id} \
      -o _tmpout \
      -p ${params.ustacks_cpus} \
      ${_args} &> ${name}_ustacks.log
    mv _tmpout/${name}_trimmed.1.alleles.tsv.gz _tmpout/${name}_trimmed.alleles.tsv.gz
    mv _tmpout/${name}_trimmed.1.snps.tsv.gz _tmpout/${name}_trimmed.snps.tsv.gz
    mv _tmpout/${name}_trimmed.1.tags.tsv.gz _tmpout/${name}_trimmed.tags.tsv.gz 
    """
}

process cstacks {

  container "quay.io/biocontainers/stacks:2.65--hdcf5f25_0"
  executor params.executor
  cpus params.cstacks_cpus
  publishDir "analysis/stacks"
  publishDir "report/logs/", pattern: "*.log"

  input:
    path stacks
    path popmap
    val args
    val meta

  output:
    path "*.{gz,tsv}", emit: data
    path "cstacks.log", emit: log
    val meta, emit: meta

  script:
    _args = args.join(" ")
    """
    cstacks -P . \
      -p ${params.cstacks_cpus} \
      -M ${popmap} \
      ${_args} &> cstacks.log
    """
}

process sstacks {

  container "quay.io/biocontainers/stacks:2.65--hdcf5f25_0"
  executor params.executor
  cpus params.sstacks_cpus
  publishDir "analysis/stacks"
  publishDir "report/logs/", pattern: "*.log"

  input:
    path catalog
    path stacks
    path popmap
    val args
    val meta

  output:
    path "*.{gz,tsv}", emit: data
    path "sstacks.log", emit: log
    val meta, emit: meta

  script:
    _args = args.join(" ")
    """
    sstacks -P . \
      -p ${params.sstacks_cpus} \
      -M ${popmap} \
      ${_args} &> sstacks.log
    """
}

process tsv2bam {

  container "quay.io/biocontainers/stacks:2.65--hdcf5f25_0"
  executor params.executor
  cpus params.tsv2bam_cpus
  publishDir "analysis/stacks"
  publishDir "report/logs/", pattern: "*.log"

  input:
    path catalog
    path stacks
    path maches
    path pe_reads
    path popmap
    val args
    val meta

  output:
    path "*.bam", emit: data
    path "tsv2bam.log", emit: log
    val meta, emit: meta

  script:
    _args = args.join(" ")
    """
    tsv2bam -P . \
      -t ${params.tsv2bam_cpus} \
      -M ${popmap} \
      --pe-reads-dir . \
      ${_args} &> tsv2bam.log
    """
}

process gstacks {

  container "quay.io/biocontainers/stacks:2.65--hdcf5f25_0"
  executor params.executor
  cpus params.gstacks_cpus
  publishDir "analysis/stacks"
  publishDir "report/logs/", pattern: "*.log"
  publishDir "report/logs/", pattern: "*.distribs"

  input:
    path bam_files
    path popmap
    val args
    val meta

  output:
    path "*.{gz,calls}", emit: data
    path "gstacks.{log,log.distribs}", emit: log
    val meta, emit: meta

  script:
    _args = args.join(" ")
    """
    gstacks -P . \
      -M ${popmap} \
      -t ${params.gstacks_cpus} \
      ${_args} &> gstacks.log
    """
}

process populations {

  container "quay.io/biocontainers/stacks:2.65--hdcf5f25_0"
  executor params.executor
  cpus params.populations_cpus
  publishDir "analysis/stacks"
  publishDir "report/logs/", pattern: "*.log"
  publishDir "report/logs/", pattern: "*.distribs"
  publishDir "report/logs/", pattern: "*sumstats_summary.tsv"

  input:
    path catalog
    path stacks
    path maches
    path pe_reads
    path bam_files
    path calls
    path popmap
    val args
    val meta

  output:
    path "*.{tsv,fa,gz,vcf,gtf}", emit: data
    tuple path("populations.{log,log.distribs}"), path("*.sumstats_summary.tsv"), emit: log
    val meta, emit: meta

  script:
    _args = args.join(" ")
    """
    populations -P . \
      -M ${popmap} \
      -t ${params.populations_cpus} \
      --vcf \
      --fasta-loci \
      --fasta-samples \
      --verbose \
      ${_args} &> populations.log
    """
}

def bibtex() {
  """
  @article{catchen2013stacks,
    title={Stacks: an analysis tool set for population genomics},
    author={Catchen, Julian and Hohenlohe, Paul A and Bassham, Susan and Amores, Angel and Cresko, William A},
    journal={Molecular ecology},
    volume={22},
    number={11},
    pages={3124--3140},
    year={2013},
    publisher={Wiley Online Library}
  }
  """
}

def citekey() {
  "@catchen2013stacks"
}
nextflow.enable.dsl = 2

include { fastqc_pe as fq_raw } from './modules/fastqc.nf'
include { citekey as fq_citekey } from './modules/fastqc.nf'
include { bibtex as fq_bib } from './modules/fastqc.nf'
include { map_input_pe as fqc_input } from './modules/fastqc.nf'
include { fastqc_pe as fq_tri } from './modules/fastqc.nf'
include { fastqc_pe as fq_cut } from './modules/fastqc.nf'

include { trimmomatic_pe as trimmomatic } from './modules/trimmomatic.nf'
include { map_input_pe as tri_input } from './modules/trimmomatic.nf'
include { citekey as tr_citekey } from './modules/trimmomatic.nf'
include { bibtex as tr_bib } from './modules/trimmomatic.nf'

include { process_radtags_pe as process_radtags } from './modules/stacks_cutadapt.nf'
include { ustacks } from './modules/stacks_cutadapt.nf'
include { cstacks } from './modules/stacks_cutadapt.nf'
include { sstacks } from './modules/stacks_cutadapt.nf'
include { tsv2bam } from './modules/stacks_cutadapt.nf'
include { gstacks } from './modules/stacks_cutadapt.nf'
include { populations } from './modules/stacks_cutadapt.nf'
include { citekey as st_citekey } from './modules/stacks_cutadapt.nf'
include { bibtex as st_bib } from './modules/stacks_cutadapt.nf'

include { multiqc } from './modules/multiqc.nf'
include { citekey as mq_citekey } from './modules/multiqc.nf'
include { bibtex as mq_bib } from './modules/multiqc.nf'

/*
  This is a quick helper process to generate a population file for Stacks
  from the CSV metadata.
*/
process make_popfile {
  container "bschiffthaler/stacks:" + params.stacks_version
  executor params.executor
  cpus 1
  publishDir "analysis/populations"

  input:
  val data

  output:
  path "populations.tsv"

  script:

  _tmp = data

  """
  echo -e "${_tmp}" > populations.tsv
  """

}


process cutadapt {
  container "quay.io/biocontainers/cutadapt:4.9--py39hff71179_1"
  publishDir "analysis/cutadapt", pattern: "*.fastq.gz"
  publishDir "report/logs/", pattern: "*.log"
  cpus 4

  input:
   tuple path("${name}_trimmed_1.fastq.gz"), path("${name}_trimmed_2.fastq.gz"), val(name)
    val meta

  output:
    tuple path("${name}_cutadapt_1.fastq.gz"), path("${name}_cutadapt_2.fastq.gz"), val(name), emit: data
    tuple path("${name}_cutadapt_1.log"), path("${name}_cutadapt_2.log"), emit: log
    val meta, emit: meta

  script:
    """
    cutadapt -a file:${baseDir}/data/overrepresented.fasta --cores=2 -m 140 -o ${name}_cutadapt_1.fastq.gz ${name}_trimmed_1.fastq.gz &> ${name}_cutadapt_1.log
    cutadapt -a file:${baseDir}/data/overrepresented.fasta --cores=2 -m 140 -o ${name}_cutadapt_2.fastq.gz ${name}_trimmed_2.fastq.gz &> ${name}_cutadapt_2.log
    """
}



workflow {

  // Main data channel reads linewise from the csv
  data = Channel
  .fromPath("${baseDir}/data/metadata.csv")
  .splitCsv(header: true)

  // Make a population description based on the CSV
  make_popfile(
    data.map(it -> {
      [it.Id + "_cutadapt", it.Population].join("\t")
    })
    .reduce {
      a, b -> [a, b].join("\n")
    }
  )
  
  // Raw data FastQC
  fq_raw("raw", fqc_input(data), data)
  // Trim and crop
  trimmomatic(tri_input(data), params.trimmomatic_trimmers, data)
  // Post-trimming FastQC
  fq_tri("trimmomatic", trimmomatic.out.data, trimmomatic.out.meta)
  // Cutadapt step
  cutadapt(trimmomatic.out.data, trimmomatic.out.meta)
  // Post-cutadapt FastQC
  fq_cut("cutadapt", cutadapt.out.data, cutadapt.out.meta)
  // Stacks 1: process radtags
  process_radtags(
    cutadapt.out.data, params.process_radtags_params, cutadapt.out.meta
  )
  // Stacks 2: ustacks
  ustacks(process_radtags.out.data, params.ustacks_params, cutadapt.out.meta)
  // Stacks 3: cstacks
  cstacks(
    ustacks.out.data.map(it -> it[0]).collect(),
    make_popfile.out,
    params.cstacks_params,
    ustacks.out.meta.collect()
  )
  // Stacks 4: sstacks
  sstacks(
    cstacks.out.data,
    ustacks.out.data.map(it -> it[0]).collect(),
    make_popfile.out,
    params.sstacks_params,
    cstacks.out.meta
  )
  // Stacks 5: tsv2bam
  tsv2bam(
    cstacks.out.data,
    ustacks.out.data.map(it -> it[0]).collect(),
    sstacks.out.data,
    process_radtags.out.data.map(it -> it[0]).collect(),
    make_popfile.out,
    params.tsv2bam_params,
    sstacks.out.meta
  )
  // Stacks 6: gstacks
  gstacks(
    tsv2bam.out.data,
    make_popfile.out,
    params.gstacks_params,
    tsv2bam.out.meta
  )
  // Stacks 7: populations
  populations(
    cstacks.out.data,
    ustacks.out.data.map(it -> it[0]).collect(),
    sstacks.out.data,
    process_radtags.out.data.map(it -> it[0]).collect(),
    tsv2bam.out.data,
    gstacks.out.data,
    make_popfile.out,
    params.populations_params,
    gstacks.out.meta
  )
  // Generate final report
  multiqc(
    fq_tri.out.data
    .concat(fq_cut.out.data)
    .concat(trimmomatic.out.log)
    .concat(cutadapt.out.log)
    .concat(gstacks.out.log)
    .concat(populations.out.log)
    .flatten()
    .collect(),
    params.multiqc_params,
    populations.out.meta
  )
}
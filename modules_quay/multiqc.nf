process multiqc {

  container "quay.io/biocontainers/multiqc:1.24.1--pyhdfd78af_0"
  publishDir "report/multiqc"
  cpus 1

  input:
    path location
    val args
    val meta

  output:
    path "multiqc_*", emit: data
    path "multiqc.log", emit: log
    val meta, emit: meta

  script:
    _args = args.join(" ")
    """
    multiqc . \
      ${_args} &> multiqc.log
    """
}

def bibtex() {
  """
  @article{ewels2016multiqc,
    title={MultiQC: summarize analysis results for multiple tools and samples in a single report},
    author={Ewels, Philip and Magnusson, M{\\aa}ns and Lundin, Sverker and K{\\"a}ller, Max},
    journal={Bioinformatics},
    volume={32},
    number={19},
    pages={3047--3048},
    year={2016},
    publisher={Oxford University Press}
  }
  """
}

def citekey() {
  "@ewels2016multiqc"
}
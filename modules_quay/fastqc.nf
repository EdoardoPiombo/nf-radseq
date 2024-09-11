process fastqc_pe {

  container "quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"
  publishDir "report/qc/${stage}", pattern: "*.{html,zip}"
  publishDir "report/logs/", pattern: "*.log"
  executor params.executor
  cpus 2

  input:
    val stage
    tuple path(read1), path(read2), val(name)
    val meta

  output:
    path "*.{zip,html}", emit: data
    path "*.log", emit: log
    val meta, emit: meta

  """
  fastqc -t 2 -o . ${read1} ${read2} &> ${name}_fastqc_${stage}.log
  """
}

/*
Function to process raw CSV input into expected 
*/
def map_input_pe(Ch) {
  Ch.map(row -> { 
    ["${baseDir}/" + row.RF, "${baseDir}/" + row.RS, row.Id] 
  })
}

def bibtex() {
  """
  @misc{andrews2010fastqc,
    title={FastQC: a quality control tool for high throughput sequence data},
    author={Andrews, Simon and others},
    year={2010},
    publisher={Babraham Bioinformatics, Babraham Institute, Cambridge, United Kingdom}
  }
  """
}

def citekey() {
  "@andrews2010fastqc"
}
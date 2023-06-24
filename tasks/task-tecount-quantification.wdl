version 1.0

# Family-level quantification of transposable elements.

task te_count {
  input {
    File bam
    File gtf_rmsk
    File gtf_gene
    String mode= "multi"
    String stranded
    String output_prefix
    String docker_image = "docker.io/polumechanos/telocal-count:latest"
  }
  command<<<
    TEcount -b ~{bam} --sortByPos --stranded ~{stranded} --mode ~{mode} --TE ~{gtf_rmsk} --GTF ~{gtf_gene} --project ~{output_prefix} --format BAM
  >>>
  output {
    File count_table= output_prefix + ".cntTable"
  }
  runtime {
    cpu: 4
    memory: "32G"
    docker: docker_image
    disks: "local-disk 500 SSD"
  }
}
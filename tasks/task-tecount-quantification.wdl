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
    if [[ '~{gtf_rmsk}' == *.gz ]]; then
      echo '------ Decompressing the repeats GTF ------' 1>&2
      gzip -dc ~{gtf_rmsk} > repeats.gtf
    else
      echo '------ No decompression needed for the repeats GTF ------' 1>&2
      cat ~{gtf_rmsk} > repeats.gtf
    fi
    if [[ '~{gtf_gene}' == *.gz ]]; then
      echo '------ Decompressing the genes GTF ------' 1>&2
      gzip -dc ~{gtf_gene} > genes.gtf
    else
      echo '------ No decompression needed for the genes GTF ------' 1>&2
      cat ~{gtf_gene} > genes.gtf
    fi
    
    TEcount -b ~{bam} --sortByPos --stranded ~{stranded} --mode ~{mode} --TE repeats.gtf --GTF genes.gtf --project ~{output_prefix} --format BAM
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
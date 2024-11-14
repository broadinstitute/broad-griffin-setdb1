version 1.0

# Loci-level quantification of transposable elements.

task te_local {
  input {
    File bam
    File te_locIndex
    File gtf_gene
    String mode= "multi"
    String stranded
    String output_prefix
    String docker_image = "docker.io/polumechanos/telocal-count:latest"
  }
  command <<<
    if [[ '~{te_locIndex}' == *.gz ]]; then
      echo '------ Decompressing the repeats GTF ------' 1>&2
      gzip -dc ~{te_locIndex} > repeats.gtf.locInd
    else
      echo '------ No decompression needed for the repeats GTF ------' 1>&2
      cat ~{te_locIndex} > repeats.gtf.locInd
    fi
    if [[ '~{gtf_gene}' == *.gz ]]; then
      echo '------ Decompressing the genes GTF ------' 1>&2
      gzip -dc ~{gtf_gene} > genes.gtf
    else
      echo '------ No decompression needed for the genes GTF ------' 1>&2
      cat ~{gtf_gene} > genes.gtf
    fi
    
    TElocal -b ~{bam} --sortByPos --stranded ~{stranded} --mode ~{mode} --TE repeats.gtf.locInd --GTF genes.gtf --project ~{output_prefix}
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
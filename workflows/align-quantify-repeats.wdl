version 1.0

import "../tasks/task-star-alignment.wdl" as star
import "../tasks/task-telocal-quantification.wdl" as te_local
import "../tasks/task-tecount-quantification.wdl" as te_count

workflow wf {
	input{
    File genome_index_tar
    Array[File] fastq1
    Array[File]? fastq2
    File gene_gtf
    File te_gtf
    File te_loci_gtf
    File sjdbGTFfile
    Int sjdbOverhang
    String twopassMode = "Basic"
    String prefix
    String stranded
    }

  call star.star_align as align {
    input: 
      fastq1=fastq1,
      fastq2=fastq2,
      genome_index_tar=genome_index_tar,
      outFileNamePrefix=prefix,
      sjdbGTFfile=sjdbGTFfile,
      sjdbOverhang=sjdbOverhang,
      twopassMode= twopassMode
  }

  call te_local.te_local as loci_uniq{
    input: 
      bam=align.bamFile,
      stranded=stranded,
      mode="uniq",
      gtf_rmsk=te_loci_gtf,
      gtf_gene=gene_gtf,
      output_prefix="${prefix}.loci_level.uniq_only"
            
  }

  call te_local.te_local as loci_multi{
    input: 
      bam=align.bamFile,
      stranded=stranded,
      mode="multi",
      gtf_rmsk=te_loci_gtf,
      gtf_gene=gene_gtf,
      output_prefix="${prefix}.loci_level.uniq_multi"
  }
  
  call te_count.te_count as family_uniq {
    input: 
      bam=align.bamFile,
      stranded=stranded,
      mode="uniq",
      gtf_rmsk=te_gtf,
      gtf_gene=gene_gtf,
      output_prefix="${prefix}.family_level.uniq_only"
  }
  
  call te_count.te_count as family_multi {
    input: 
      bam=align.bamFile,
      stranded=stranded,
      mode="multi",
      gtf_rmsk=te_gtf,
      gtf_gene=gene_gtf,
      output_prefix="${prefix}.family_level.uniq_multi"
  }  

  output{
    File align_log=align.logFinalOut
    File bam_repeats_optimized=align.bamFile
    File family_level_counts_uniq_mappers=family_uniq.count_table
  	File family_level_counts_multi_mappers=family_multi.count_table
    File loci_level_counts_uniq_mappers=loci_uniq.count_table
  	File loci_level_counts_multi_mappers=loci_multi.count_table
  }

}
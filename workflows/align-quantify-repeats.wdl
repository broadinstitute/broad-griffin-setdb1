version 1.0

import "../tasks/task-star-alignment.wdl" as star
import "../tasks/task-telocal-quantification.wdl" as te_local
import "../tasks/task-tecount-quantification.wdl" as te_count
import "../tasks/task-generate-rna-tracks.wdl" as generate_tracks

workflow wf {
	input{
    File genome_index_tar
    File chromosome_sizes_file
    Array[File] fastq1
    Array[File]? fastq2
    File te_gtf
    File te_locIndex
    File genes_gtf
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
      outFileNamePrefix="${prefix}.",
      genes_gtf=genes_gtf,
      sjdbOverhang=sjdbOverhang,
      twopassMode= twopassMode
  }

  call te_local.te_local as loci_uniq{
    input: 
      bam=align.bamFile,
      stranded=stranded,
      mode="uniq",
      te_locIndex=te_locIndex,
      gtf_gene=genes_gtf,
      output_prefix="${prefix}.loci_level.uniq_only"
            
  }

  call te_local.te_local as loci_multi{
    input: 
      bam=align.bamFile,
      stranded=stranded,
      mode="multi",
      te_locIndex=te_locIndex,
      gtf_gene=genes_gtf,
      output_prefix="${prefix}.loci_level.uniq_multi"
  }
  
  call te_count.te_count as family_uniq {
    input: 
      bam=align.bamFile,
      stranded=stranded,
      mode="uniq",
      gtf_rmsk=te_gtf,
      gtf_gene=genes_gtf,
      output_prefix="${prefix}.family_level.uniq_only"
  }
  
  call te_count.te_count as family_multi {
    input: 
      bam=align.bamFile,
      stranded=stranded,
      mode="multi",
      gtf_rmsk=te_gtf,
      gtf_gene=genes_gtf,
      output_prefix="${prefix}.family_level.uniq_multi"
  }

  call generate_tracks.generate_tracks as tracks {
    input: 
      bam=align.bamFile,
      chromosome_sizes_file=chromosome_sizes_file,
      library_size=align.rna_aligned_reads,
      output_prefix=prefix
  }

  output{
    File align_log=align.logFinalOut
    File bam_repeats_optimized=align.bamFile
    Int rna_input_reads=align.rna_input_reads
    Int rna_aligned_reads=align.rna_aligned_reads
    Int rna_aligned_uniquely=align.rna_aligned_uniquely
    Int rna_aligned_multimap=align.rna_aligned_multimap
    Int rna_unaligned_reads=align.rna_unaligned_reads
    File family_level_counts_uniq_mappers=family_uniq.count_table
  	File family_level_counts_multi_mappers=family_multi.count_table
    File loci_level_counts_uniq_mappers=loci_uniq.count_table
  	File loci_level_counts_multi_mappers=loci_multi.count_table
    File rna_fwd_strand_bw_track=tracks.fwd_bw
    File rna_rev_strand_bw_track=tracks.rev_bw
    File rna_bw_tracks=tracks.bw
    File rna_unique_bw_tracks=tracks.unique_bw
  }

}

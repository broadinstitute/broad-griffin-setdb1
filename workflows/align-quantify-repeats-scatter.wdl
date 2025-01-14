version 1.0

import "../workflows/align-quantify-repeats.wdl" as align_quantify_repeats

workflow wf {
	input{
    File genome_index_tar
    File chromosome_sizes_file
    Array[Array[File]] fastq1
    #Array[Array[File]]? fastq2
    File te_gtf
    File te_locIndex
    File genes_gtf
    Int sjdbOverhang
    String twopassMode = "Basic"
    String prefix
    String stranded
    }

  scatter(index in range(length(fastq1))) {
      call align_quantify_repeats.wf as wf_align_quantify_repeats {
            input: 
              fastq1=fastq1[index],
              #fastq2=fastq2[index],
              genome_index_tar=genome_index_tar,
              genes_gtf=genes_gtf,
              sjdbOverhang=sjdbOverhang,
              twopassMode=twopassMode
        }
  }

  output{
    Array[File] align_log=wf_align_quantify_repeats.align_log
    Array[File] bam_repeats_optimized=wf_align_quantify_repeats.bam_repeats_optimized
    Array[Int] rna_input_reads=wf_align_quantify_repeats.rna_input_reads
    Array[Int] rna_aligned_reads=wf_align_quantify_repeats.rna_aligned_reads
    Array[Int] rna_aligned_uniquely=wf_align_quantify_repeats.rna_aligned_uniquely
    Array[Int] rna_aligned_multimap=wf_align_quantify_repeats.rna_aligned_multimap
    Array[Int] rna_unaligned_reads=wf_align_quantify_repeats.rna_unaligned_reads
    Array[File] family_level_counts_uniq_mappers=wf_align_quantify_repeats.family_level_counts_uniq_mappers
  	Array[File] family_level_counts_multi_mappers=wf_align_quantify_repeats.family_level_counts_multi_mappers
    Array[File] loci_level_counts_uniq_mappers=wf_align_quantify_repeats.loci_level_counts_uniq_mappers
  	Array[File] loci_level_counts_multi_mappers=wf_align_quantify_repeats.loci_level_counts_multi_mappers
    Array[File] rna_fwd_strand_bw_track=wf_align_quantify_repeats.rna_fwd_strand_bw_track
    Array[File] rna_rev_strand_bw_track=wf_align_quantify_repeats.rna_rev_strand_bw_track
    Array[File] rna_bw_tracks=wf_align_quantify_repeats.rna_bw_tracks
    Array[File] rna_unique_bw_tracks=wf_align_quantify_repeats.rna_unique_bw_tracks
  }
}

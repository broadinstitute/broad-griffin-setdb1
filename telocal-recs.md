Recommendations for TElocal input files

TElocal can perform transposable element quantification from alignment results (e.g. BAM files) generated from a variety of programs. Given the variety of experimental systems, we could not provide an optimal alignment strategy for every approach. Therefore, we recommend that users identify the optimal parameters for their particular genome and alignment program in order to get the best results.

When optimizing the alignment parameters, we recommend taking these points into consideration:

Allowing sufficient number of multi-mappers during alignment

Most alignment programs provide only 1 alignment per read by default. We recommend reporting multiple alignments per read. We have found that reporting a maximum of 100 alignments per read provides an optimal compromise between the size of the alignment file and recovery of multi-mappers in many genome builds. However, we highly suggest that users optimize this parameter for their particular experiment, as this could significantly improve the quality of transposable element quantification.

Paired end sequencing input

For paired-end libraries, it is recommended that only alignments from properly paired reads are present in the input BAM file. I.e., each read 1 alignment should only have a single read 2 alignment. For example, if read 1 matched 3 genomic locations (A, B, C), then if read 2 also match 3 genomic locations (A', B', C'), then all three pairs of alignments could be used (and should be in the BAM file). However, if alignment C of read 1 was matched with more than one alignment of read 2 (e.g. C' and C*), then alignment C should be discarded (as there are unmatched alignments between read 1 and read 2). STAR only outputs properly paired alignments by default, while Bowtie2 requires the --no-mixed parameter to be used.

Specific recommendations when using STAR

STAR utilizes two parameters for optimal identification of multi-mappers --outFilterMultimapNmax and --outAnchorMultimapNmax. The author of STAR recommends that --winAnchorMultimapNmax should be set at twice the value used in --outFilterMultimapNmax, but no less than 50. In our study, we used the same number for both parameters (100), and found negligible differences in identifying multi-mappers. Upon further discussion with the author of STAR, we recommend that setting the same value for --winAnchorMultimapNmax and --outFilterMultimapNmax, though we highly suggest users test multiple values of --winAnchorMultimapNmax to identify the optimal value for their experiment.
version 1.0

# Generate RNA tracks for visualization in the UCSC Genome Browser.
# https://www.biostars.org/p/92935/
# This task assumes you are usiing the Illumina strand specific protocol that produces R2-R1 read orientation.

task generate_tracks{
    input{
        File bam
        File chromosome_sizes_file
        Int library_size
        String output_prefix
        Int cpus = 2
        Int memory_gb = 16
        Int disk_gb = 500
        String docker_image = "polumechanos/generate-rna-tracks:repeats-pipeline"
    }

    # Determining disk type base on the size of disk.
    String disk_type = if disk_gb > 375 then "SSD" else "LOCAL"

    command <<<
        scale_factor_merged=$(bc <<< "scale=6;1000000/~{library_size}")
        # Dividing by 2 because we are splitting by strand
        scale_factor_strand=$(bc <<< "scale=6;1000000/~{library_size}/2")

        # To get the file for transcripts that originated from the forward strand:

        # include reads that are 2nd in a pair (128);
        # exclude reads that are mapped to the reverse strand (16)
        samtools view -b -f 128 -F 16 ~{bam} > a.fwd1.bam

        # exclude reads that are mapped to the reverse strand (16) and
        # first in a pair (64): 64 + 16 = 80
        samtools view -b -f 80 ~{bam} > a.fwd2.bam

        # combine the temporary files
        samtools merge -f ~{output_prefix}_fwd.bam a.fwd1.bam a.fwd2.bam

        # index the filtered BAM file
        samtools index ~{output_prefix}_fwd.bam

        # run bedtools genomecov
        bedtools genomecov -g ~{chromosome_sizes_file} -ibam ~{output_prefix}_fwd.bam -bg -scale $scale_factor_strand > ~{output_prefix}_cpm_fwd.bedgraph
        bedGraphToBigWig ~{output_prefix}_cpm_fwd.bedgraph ~{chromosome_sizes_file} ~{output_prefix}_cpm_fwd.bw

        # remove the temporary files
        rm a.fwd*.bam ~{output_prefix}_cpm_fwd.bedgraph
        
        # To get the file for transcripts that originated from the reverse strand:

        # include reads that map to the reverse strand (128)
        # and are second in a pair (16): 128 + 16 = 144
        samtools view -b -f 144 ~{bam} > a.rev1.bam

        # include reads that are first in a pair (64), but
        # exclude those ones that map to the reverse strand (16)
        samtools view -b -f 64 -F 16 ~{bam} > a.rev2.bam

        # merge the temporary files
        samtools merge -f ~{output_prefix}_rev.bam a.rev1.bam a.rev2.bam

        # index the merged, filtered BAM file
        samtools index ~{output_prefix}_rev.bam

        # run bedtools genomecov
        bedtools genomecov -g ~{chromosome_sizes_file} -ibam ~{output_prefix}_rev.bam -bg -scale $scale_factor_strand > ~{output_prefix}_cpm_rev.bedgraph
        bedGraphToBigWig ~{output_prefix}_cpm_rev.bedgraph ~{chromosome_sizes_file} ~{output_prefix}_cpm_rev.bw

        # remove temporary files
        rm a.rev*.bam ~{output_prefix}_cpm_rev.bedgraph

        #  RNA bigwig
        bedtools genomecov -g ~{chromosome_sizes_file} -ibam ~{bam} -bg -scale $scale_factor_merged > ~{output_prefix}_cpm.bedgraph
        bedGraphToBigWig ~{output_prefix}_cpm.bedgraph ~{chromosome_sizes_file} ~{output_prefix}_cpm.bw
    >>>

    output{
        File fwd_bw = output_prefix + "_cpm_fwd.bw"
        File rev_bw = output_prefix + "_cpm_rev.bw"
        File bw = output_prefix + "_cpm.bw"
    }
    runtime{
        cpu: cpus
        memory: "~{memory_gb}G"
        docker: docker_image
        disks: "local-disk ~{disk_gb} ~{disk_type}"
    }

    parameter_meta{
        bam: "The input BAM file"
        chromosome_sizes_file: "The file containing the chromosome sizes"
        library_size: "The library size"
        output_prefix: "The prefix for the output files"
        cpus: "The number of CPUs to use"
        memory_gb: "The amount of memory to use"
        disk_gb: "The amount of disk space to use"
        docker_image: "The docker image to use for the task"
    }
}
















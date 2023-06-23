version 1.0

# Align reads to the genome using repeats-aware parameters using bowtie2 for chromatin and HiSat2 for RNA

task smartmap {
    input {
        Array[File] fastq1
        Array[File]? fastq2
        #File index_tar # For the future. don't think is important now.
        Array[File]+ indexFiles
        String outFileNamePrefix
        String outSAMtype = "BAM SortedByCoordinate"
        String readFilesCommand = "zcat"
        Int outBAMcompression = 1
        String? twopassMode
        String outSAMunmapped = "Within KeepPairs"
      
        File? sjdbGTFfile
        Int? sjdbOverhang
        Int chimSegmentMin= 0
        Int outFilterMultimapNmax = 100
        Int winAnchorMultimapNmax = 100
        String alignEndsType = "EndToEnd"
        String alignEndsProtrude = "100 DiscordantPair"
        Float outFilterScoreMinOverLread = 0.4
        Float outFilterMatchNminOverLread = 0.4
        String outSAMattributes = "All"
        String outSAMattrIHstart = "0"
        #String outSAMstrandField = "intronMotif" # This is a placeholder for the future. Unclear if we need it.


        Int cpus = 16
        String? memory = "40G"
        # 20 minute initialization + time reading in index (1 minute per G) + time aligning data.
        #Int timeMinutes = 20 + ceil(size(indexFiles, "G")) + ceil(size(flatten([inputR1, inputR2]), "G") * 300 / runThreadN)
        String docker_image = "docker.io/polumechanos/starmap"
    }

    # Use a margin of 30% index size. Real memory usage is ~30 GiB for a 27 GiB index. 
    Int memoryGb = 40
    # For some reason doing above calculation inside a string does not work.
    # So we solve it with an optional memory string and using select_first
    # in the runtime section.

    #TODO: Could be extended for all possible output extensions.
    Map[String, String] samOutputNames = {"BAM SortedByCoordinate": "sortedByCoord.out.bam"}

    command <<<
        set -e
        # Create a repeats-aware bam file for chromatin data.
        SmartMapPrep -k 51 -I 100 -L 2000 -p ~{cpus} -x [Bowtie2 index] -o [output prefix] -1 [R1 fastq] -2 [R2 fastq]

        # Create a repeats-aware bam file for transcriptomic data.
        SmartMapRNAPrep -k 51 -I 100 -L 2000 -p ~{cpus} -x [Bowtie2 index] -o [output prefix] -1 [R1 fastq] -2 [R2 fastq]

        # After prepping we can run SmartMap
        # -c : Flag for continuous output bedgraphs. Default off.
        # -S : Flag for strand-specific mode. Default off.
        # -r : Flag for read output mode with weights. Default off.
        SmartMap [options] -m 50 -s 0 -i 1 -v 1 -l 1 -g [genome length file] -o [output prefix] [BED or BED.gz file input(s)]
    >>>

    output {
        File bamFile = outFileNamePrefix + "Aligned." +  samOutputNames[outSAMtype]
        File logFinalOut = outFileNamePrefix + "Log.final.out"
    }

    runtime {
        cpu: runThreadN
        memory: select_first([memory, "${memoryGb}G"])
        docker: docker_image
        disks: "local-disk 500 SSD"
    }

    parameter_meta {
        # inputs
        indexFiles: {description: "The star index files.", category: "required"}
        outFileNamePrefix: {description: "The prefix for the output files. May include directories.", category: "required"}
        outSAMtype: {description: "The type of alignment file to be produced. Currently only `BAM SortedByCoordinate` is supported.", category: "advanced"}
        readFilesCommand: {description: "Equivalent to star's `--readFilesCommand` option.", category: "advanced"}
        outBAMcompression: {description: "The compression level of the output BAM.", category: "advanced"}
        outFilterScoreMinOverLread: {description: "Equivalent to star's `--outFilterScoreMinOverLread` option.", category: "advanced"}
        outFilterMatchNminOverLread: {description: "Equivalent to star's `--outFilterMatchNminOverLread` option.", category: "advanced"}
        twopassMode: {description: "Equivalent to star's `--twopassMode` option.", category: "advanced"}
        outSAMunmapped: {description: "Equivalent to star's `--outSAMunmapped` option.", category: "advanced"}
        runThreadN: {description: "The number of threads to use.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        docker_image: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}
    }
}


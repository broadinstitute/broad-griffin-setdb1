version 1.0

# Copyright (c) 2017 Leiden University Medical Center
#
# Modified by Eugenio Mattei

# Align reads to the genome using repeats-aware parameters. The parameters are chose using the 
# recommendations from Squire and TElocal.

# when passing the gtf file and the overhang int length, the twopassmode input variable needs to be set to 'basic'

task Star {
    input {
        File fastq1
        File? fastq2
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


        Int runThreadN = 4
        String? memory = "40G"
        # 20 minute initialization + time reading in index (1 minute per G) + time aligning data.
        #Int timeMinutes = 20 + ceil(size(indexFiles, "G")) + ceil(size(flatten([inputR1, inputR2]), "G") * 300 / runThreadN)
        String docker_image = "docker.io/polumechanos/star:2.10.b"
    }

    # Use a margin of 30% index size. Real memory usage is ~30 GiB for a 27 GiB index. 
    Int memoryGb = 40
    # For some reason doing above calculation inside a string does not work.
    # So we solve it with an optional memory string and using select_first
    # in the runtime section.

    #TODO: Could be extended for all possible output extensions.
    Map[String, String] samOutputNames = {"BAM SortedByCoordinate": "sortedByCoord.out.bam"}

    command {
        set -e
        mkdir -p `dirname ${outFileNamePrefix}`
        STAR \
        --readFilesIn ${fastq1} ${fastq2} \
        --outFileNamePrefix ${outFileNamePrefix} \
        --genomeDir ${sub(indexFiles[0], basename(indexFiles[0]), "")} \
        --outSAMtype ${outSAMtype} \
        --outBAMcompression ${outBAMcompression} \
        ${"--readFilesCommand " + readFilesCommand} \
        ${"--outFilterScoreMinOverLread " + outFilterScoreMinOverLread} \
        ${"--outFilterMatchNminOverLread " + outFilterMatchNminOverLread} \
        ${"--outSAMunmapped " + outSAMunmapped} \
        ${"--runThreadN " + runThreadN} \
        ${"--twopassMode " + twopassMode} \
        ${"--sjdbGTFfile " + sjdbGTFfile} \
        ${"--sjdbOverhang " + sjdbOverhang} \
        ${"--chimSegmentMin " + chimSegmentMin} \
        ${"--outFilterMultimapNmax " + outFilterMultimapNmax} \
        ${"--winAnchorMultimapNmax " + winAnchorMultimapNmax} \
        ${"--alignEndsType " + alignEndsType} \
        ${"--alignEndsProtrude " + alignEndsProtrude} \
        ${"--outSAMattributes " + outSAMattributes} \
        ${"--outSAMattrIHstart " + outSAMattrIHstart}
    }

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

        # outputs
        bamFile: {description: "Alignment file."}
        logFinalOut: {description: "Log information file."}
    }
}


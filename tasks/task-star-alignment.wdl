version 1.0

# Copyright (c) 2017 Leiden University Medical Center
#
# Modified by Eugenio Mattei

# Align reads to the genome using repeats-aware parameters. The parameters are chose using the 
# recommendations from Squire and TElocal.

# when passing the gtf file and the overhang int length, the twopassmode input variable needs to be set to 'basic'

task star_align {
    input {
        Array[File] fastq1
        Array[File]? fastq2
        #File index_tar # For the future. don't think is important now.
        File genome_index_tar
        String outFileNamePrefix
        String outSAMtype = "BAM SortedByCoordinate"
        String readFilesCommand = "zcat"
        Int outBAMcompression = 1
        String twopassMode = "Basic"
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


        Int cpus = 4
        Int memory_gb = 40
        Int disk_size_gb = 500

        String docker_image = "docker.io/polumechanos/star:2.10.b"
    }



    #TODO: Could be extended for all possible output extensions.
    Map[String, String] samOutputNames = {"BAM SortedByCoordinate": "sortedByCoord.out.bam"}

    command {
        set -e

        tar -xvzf ${genome_index_tar} --no-same-owner -C ./

        if [[ '${sjdbGTFfile}' == *.gz ]]; then
            echo '------ Decompressing the GTF ------' 1>&2
            gunzip -c ${sjdbGTFfile} > genes.gtf
        else
            echo '------ No decompression needed for the GTF ------' 1>&2
            cat ${sjdbGTFfile} > genes.gtf
        fi

        mkdir -p `dirname ${outFileNamePrefix}`
        STAR \
        --readFilesIn ${sep="," fastq1} ${sep=","fastq2} \
        --outFileNamePrefix ${outFileNamePrefix} \
        --genomeDir ./ \
        --outSAMtype ${outSAMtype} \
        --outBAMcompression ${outBAMcompression} \
        ${"--readFilesCommand " + readFilesCommand} \
        ${"--outFilterScoreMinOverLread " + outFilterScoreMinOverLread} \
        ${"--outFilterMatchNminOverLread " + outFilterMatchNminOverLread} \
        ${"--outSAMunmapped " + outSAMunmapped} \
        ${"--runThreadN " + cpus} \
        ${"--twopassMode " + twopassMode} \
        ${if defined(sjdbGTFfile) then "--sjdbGTFfile genes.gtf" else ""} \
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
        cpu: "${cpus}"
        memory: "${memory_gb}G"
        docker: docker_image
        disks: "local-disk ${disk_size_gb} SSD"
    }

    parameter_meta {
        # inputs
        fastq1: {description: "The first read file.", type: "array[file]", category: "required"}
        fastq2: {description: "The second read file.", type: "array[file]", category: "optional"}
        genome_index_tar: {description: "The tarball containing the genome index.", type: "file", category: "required"}
        outFileNamePrefix: {description: "The prefix for the output files. May include directories.", category: "required"}
        outSAMtype: {description: "The type of alignment file to be produced. Currently only `BAM SortedByCoordinate` is supported.", category: "advanced"}
        readFilesCommand: {description: "Equivalent to star's `--readFilesCommand` option.", category: "advanced"}
        outBAMcompression: {description: "The compression level of the output BAM.", category: "advanced"}
        outFilterScoreMinOverLread: {description: "Equivalent to star's `--outFilterScoreMinOverLread` option.", category: "advanced"}
        outFilterMatchNminOverLread: {description: "Equivalent to star's `--outFilterMatchNminOverLread` option.", category: "advanced"}
        twopassMode: {description: "Equivalent to star's `--twopassMode` option.", category: "advanced"}
        outSAMunmapped: {description: "Equivalent to star's `--outSAMunmapped` option.", category: "advanced"}
        cpus: {description: "The number of threads to use.", category: "advanced"}
        memory_gb: {description: "The amount of memory this job will use.", category: "advanced"}
        disk_size_gb: {description: "The amount of disk space this job will use.", category: "advanced"}
        docker_image: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        bamFile: {description: "Alignment file."}
        logFinalOut: {description: "Log information file."}
    }
}


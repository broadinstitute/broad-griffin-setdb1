version 1.0

workflow wf {
    call GenomeGenerate
}


task GenomeGenerate {
    input {
        String genomeDir = "STAR_index"
        File referenceFasta

        File? referenceGtf
        Int? sjdbOverhang

        Int threads = 4
        String memory = "32G"
        Int timeMinutes = ceil(size(referenceFasta, "G") * 240 / threads)
        String dockerImage = "quay.io/biocontainers/star:2.7.3a--0"
    }

    command {
        set -e
        mkdir -p ~{genomeDir}
        STAR \
        --runMode genomeGenerate \
        --runThreadN ~{threads} \
        --genomeDir ~{genomeDir} \
        --genomeFastaFiles ~{referenceFasta} \
        ~{"--sjdbGTFfile " + referenceGtf} \
        ~{"--sjdbOverhang " + sjdbOverhang}
    }

    output {
        File chrLength = "~{genomeDir}/chrLength.txt"
        File chrNameLength = "~{genomeDir}/chrNameLength.txt"
        File chrName = "~{genomeDir}/chrName.txt"
        File chrStart = "~{genomeDir}/chrStart.txt"
        File genome = "~{genomeDir}/Genome"
        File genomeParameters = "~{genomeDir}/genomeParameters.txt"
        File sa = "~{genomeDir}/SA"
        File saIndex = "~{genomeDir}/SAindex"
        File? exonGeTrInfo = "~{genomeDir}/exonGeTrInfo.tab"
        File? exonInfo = "~{genomeDir}/exonInfo.tab"
        File? geneInfo = "~{genomeDir}/geneInfo.tab"
        File? sjdbInfo = "~{genomeDir}/sjdbInfo.txt"
        File? sjdbListFromGtfOut = "~{genomeDir}/sjdbList.fromGTF.out.tab"
        File? sjdbListOut = "~{genomeDir}/sjdbList.out.tab"
        File? transcriptInfo = "~{genomeDir}/transcriptInfo.tab"
        Array[File] starIndex = select_all([chrLength, chrNameLength, chrName,
                                            chrStart, genome, genomeParameters,
                                            sa, saIndex, exonGeTrInfo, exonInfo,
                                            geneInfo, sjdbInfo, sjdbListFromGtfOut,
                                            sjdbListOut, transcriptInfo])
    }

    runtime {
        cpu: threads
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
        disks: "local-disk 500 SSD"
    }

    parameter_meta {
        # inputs
        genomeDir: {description:"The directory the STAR index should be written to.", categroy: "common"}
        referenceFasta: {description: "The reference Fasta file.", category: "required"}
        referenceGtf: {description: "The reference GTF file.", category: "common"}
        sjdbOverhang: {description: "Equivalent to STAR's `--sjdbOverhang` option.", category: "advanced"}
        threads: {description: "The number of threads to use.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        chrLength: {description: "Text chromosome lengths file."}
        chrNameLength: {description: "Text chromosome name lengths file."}
        chrName: {description: "Text chromosome names file."}
        chrStart: {description: "Chromosome start sites file."}
        genome: {description: "Binary genome sequence file."}
        genomeParameters: {description: "Genome parameters file."}
        sa: {description: "Suffix arrays file."}
        saIndex: {description: "Index file of suffix arrays."}
        exonGeTrInfo: {description: "Exon, gene and transcript information file."}
        exonInfo: {description: "Exon information file."}
        geneInfo: {description: "Gene information file."}
        sjdbInfo: {description: "Splice junctions coordinates file."}
        sjdbListFromGtfOut: {description: "Splice junctions from input GTF file."}
        sjdbListOut: {description: "Splice junction list file."}
        transcriptInfo: {description: "Transcripts information file."}
        starIndex: {description: "A collection of all STAR index files."}
    }
}
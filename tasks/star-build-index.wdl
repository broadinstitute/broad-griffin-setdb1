version 1.0

task GenomeGenerate {
    input {
        String genome_dir = "STAR_index"
        File reference_fasta

        Int? sjdbOverhang
        File? reference_gtf
        String? prefix_output= "star-index"

        Int threads = 16
        String memory = "32G"
        Int time_minutes = ceil(size(reference_fasta, "G") * 240 / threads)
        String docker_image = "docker.io/polumehcanos/star:2.10.b"
    }

    command {
        set -e

        mkdir -p ~{genome_dir}

        STAR \
        --runMode genomeGenerate \
        --runThreadN ~{threads} \
        --genome_dir ~{genome_dir} \
        --genomeFastaFiles ~{reference_fasta} \
        ~{"--sjdbGTFfile " + reference_gtf} \
        ~{if defined(sjdbOverhang) then "--sjdbOverhang ~{sjdbOverhang}" else ""}

        cd ~{genome_dir}
        tar -cvzf ../~{prefix_output}.tar.gz .
    }

    output {
        File chrLength = "~{genome_dir}/chrLength.txt"
        File chrNameLength = "~{genome_dir}/chrNameLength.txt"
        File chrName = "~{genome_dir}/chrName.txt"
        File chrStart = "~{genome_dir}/chrStart.txt"
        File genome = "~{genome_dir}/Genome"
        File genomeParameters = "~{genome_dir}/genomeParameters.txt"
        File sa = "~{genome_dir}/SA"
        File saIndex = "~{genome_dir}/SAindex"
        File index_tarred = "~{prefix_output}.tar.gz"
        File? exonGeTrInfo = "~{genome_dir}/exonGeTrInfo.tab"
        File? exonInfo = "~{genome_dir}/exonInfo.tab"
        File? geneInfo = "~{genome_dir}/geneInfo.tab"
        File? sjdbInfo = "~{genome_dir}/sjdbInfo.txt"
        File? sjdbListFromGtfOut = "~{genome_dir}/sjdbList.fromGTF.out.tab"
        File? sjdbListOut = "~{genome_dir}/sjdbList.out.tab"
        File? transcriptInfo = "~{genome_dir}/transcriptInfo.tab"
        Array[File] starIndex = select_all([chrLength, chrNameLength, chrName,
                                            chrStart, genome, genomeParameters,
                                            sa, saIndex, exonGeTrInfo, exonInfo,
                                            geneInfo, sjdbInfo, sjdbListFromGtfOut,
                                            sjdbListOut, transcriptInfo])
    }

    runtime {
        cpu: threads
        memory: memory
        time_minutes: time_minutes
        docker: docker_image
        disks: "local-disk 500 SSD"
    }

    parameter_meta {
        # inputs
        genome_dir: {description:"The directory the STAR index should be written to.", categroy: "common"}
        reference_fasta: {description: "The reference Fasta file.", category: "required"}
        reference_gtf: {description: "The reference GTF file.", category: "common"}
        sjdbOverhang: {description: "Equivalent to STAR's `--sjdbOverhang` option.", category: "advanced"}
        threads: {description: "The number of threads to use.", category: "advanced"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        time_minutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        docker_image: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

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
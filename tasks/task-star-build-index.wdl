version 1.0

task GenomeGenerate {
    input {
        File reference_fasta
        File reference_gtf

        Int read_length=100
        
        String prefix_output="star-index"

        Int threads = 16
        String memory = "32G"
        Int time_minutes = ceil(size(reference_fasta, "G") * 240 / threads)
        String docker_image = "docker.io/polumechanos/star:2.10.b"
    }

    command <<<
        set -e

        mkdir -p temp_index

        overhang=$(bc <<< "~{read_length}-1")

        STAR \
        --runMode genomeGenerate \
        --runThreadN ~{threads} \
        --genomeDir temp_index \
        --genomeFastaFiles ~{reference_fasta} \
        --sjdbGTFfile ~{reference_gtf} \
        --sjdbOverhang "$overhang"

        cd temp_index
        tar -cvzf ../~{prefix_output}.tar.gz .
    >>>

    output {
        File genome_index_tar = "~{prefix_output}.tar.gz"
    }

    runtime {
        cpu: threads
        memory: memory
        time_minutes: time_minutes
        docker: docker_image
        disks: "local-disk 500 SSD"
    }

    parameter_meta {
        reference_fasta: "Reference genome in fasta format"
        reference_gtf: "Reference annotation in GTF format"
        read_length: "Read length"
        prefix_output: "Prefix for output files"
    }
}
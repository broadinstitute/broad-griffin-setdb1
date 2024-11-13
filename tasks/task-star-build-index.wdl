version 1.0

task GenomeGenerate {
    input {
        File reference_fasta
        File reference_gtf

        Int read_length=100
        
        String prefix_output="star-index"

        Int threads = 16
        Int memory = 32
        String docker_image = "docker.io/polumechanos/star:2.10.b"
    }

    command <<<
        set -e

        mkdir -p temp_index

        if [[ '~{reference_fasta}' == *.gz ]]; then
            echo '------ Decompressing the genome ------' 1>&2
            gunzip -c ~{reference_fasta} > genome.fa
        else
            echo '------ No decompression needed for the genome ------' 1>&2
            cat ~{reference_fasta} > genome.fa
        fi

        if [[ '~{reference_gtf}' == *.gz ]]; then
            echo '------ Decompressing the GTF ------' 1>&2
            gunzip -c ~{reference_gtf} > genes.gtf
        else
            echo '------ No decompression needed for the GTF ------' 1>&2
            cat ~{reference_gtf} > genes.gtf
        fi


        overhang=$(bc <<< "~{read_length}-1")

        STAR \
        --runMode genomeGenerate \
        --runThreadN ~{threads} \
        --genomeDir temp_index \
        --genomeFastaFiles genome.fa \
        --sjdbGTFfile genes.gtf \
        --sjdbOverhang "$overhang"

        cd temp_index
        tar -cvzf ../~{prefix_output}.tar.gz .
    >>>

    output {
        File genome_index_tar = "~{prefix_output}.tar.gz"
    }

    runtime {
        cpu: "~{threads}"
        memory: "~{memory}G"
        docker: "~{docker_image}"
        disks: "local-disk 500 SSD"
    }

    parameter_meta {
        reference_fasta: "Reference genome in fasta format"
        reference_gtf: "Reference annotation in GTF format"
        read_length: "Read length"
        prefix_output: "Prefix for output files"
    }
}
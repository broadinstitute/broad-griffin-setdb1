set -e

mkdir -p mm10/genome mm10/genes mm10/repeats

# Download the chromosome sizes from UCSC
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/p6/mm10.p6.chrom.sizes -O mm10/genome/mm10.p6.chrom.sizes

# Download the chromosome FASTAs from UCSC
mkdir tmp_mm10; cd tmp_mm10
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/p6/mm10.p6.chromFa.tar.gz
tar xvzf mm10.p6.chromFa.tar.gz

# Remove chromosome M and create a unique FASTA file.
# Chromosome will be ordered following the chromosome size ordering.
cut -f1 ../mm10/genome/mm10.p6.chrom.sizes | awk '{print "chroms/"$1".fa"}' | xargs cat | gzip -c > ../mm10/genome/mm10.p6.ucsc.fa.gz

# Clean up
cd ..; rm -r tmp_mm10

# Download the GENCODE gene annotations
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz -O mm10/genes/gencode.vM25.mm10.p6.annotation.gtf.gz

# Download TEtranscript pre-indexed files. The repeat masker is filtered and contains less entries than the UCSC one. 
# tRNA are removed and other simple and short entries. Exact details are unclear.
# Note that the script will remove most simple repetitive sequences and short non-coding RNA (e.g. tRNA)
# All these files need to be uncompressed using the `gzip -d` command. TElocal does not accept compressed files.
wget https://labshare.cshl.edu/shares/mhammelllab/www-data/TEtranscripts/TE_GTF/mm10_rmsk_TE.gtf.gz -O mm10/repeats/mm10_rmsk_TE.gtf.gz
wget https://labshare.cshl.edu/shares/mhammelllab/www-data/TElocal/annotation_tables/mm10_rmsk_TE.gtf.locInd.locations.gz -O mm10/repeats/mm10_rmsk_TE.gtf.locInd.locations.gz

# Generate md5 checksums for all downloaded files
md5sum mm10/genome/mm10.p6.ucsc.fa.gz > mm10/genome/mm10.p6.ucsc.fa.gz.md5
md5sum mm10/genome/mm10.p6.chrom.sizes > mm10/genome/mm10.p6.chrom.sizes.md5
md5sum mm10/genes/gencode.vM25.mm10.p6.annotation.gtf.gz > mm10/genes/gencode.vM25.mm10.p6.annotation.gtf.gz.md5
md5sum mm10/repeats/mm10_rmsk_TE.gtf.gz > mm10/repeats/mm10_rmsk_TE.gtf.gz.md5
md5sum mm10/repeats/mm10_rmsk_TE.gtf.locInd.locations.gz > mm10/repeats/mm10_rmsk_TE.gtf.locInd.locations.gz.md5
md5sum mm10/repeats/mm10_rmsk_TE.gtf.locInd.gz > mm10/repeats/mm10_rmsk_TE.gtf.locInd.gz.md5
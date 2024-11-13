set -e

mkdir -p hg38/genome hg38/genes hg38/repeats

# Download the chromosome sizes from UCSC
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/p14/hg38.p14.chrom.sizes -O hg38/genome/hg38.p14.chrom.sizes

# Download the chromosome FASTAs from UCSC
mkdir tmp; cd tmp
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/p14/hg38.p14.chromFa.tar.gz
tar xvzf hg38.p14.chromFa.tar.gz

# Remove chromosome M and create a unique FASTA file.
# Chromosome will be ordered following the chromosome size ordering.
cut -f1 ../hg38/genome/hg38.p14.chrom.sizes | awk '{print "chroms/"$1".fa"}' | xargs cat  | gzip -c > ../hg38/genome/hg38.p14.ucsc.fa.gz

# Clean up
cd ..; rm -r tmp

# Download the GENCODE complete gene annotations.
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.chr_patch_hapl_scaff.annotation.gtf.gz -O hg38/genes/gencode.v47.hg38.p14.chr_patch_hapl_scaff.annotation.gtf.gz

# Download TEtranscript pre-indexed files. The repeat masker is filtered and contains less entries than the UCSC one. 
# tRNA are removed and other simple and short entries. Exact details are unclear.
# Note that the script will remove most simple repetitive sequences and short non-coding RNA (e.g. tRNA)
# All these files need to be uncompressed using the `gzip -d` command. TElocal does not accept compressed files.
wget https://labshare.cshl.edu/shares/mhammelllab/www-data/TEtranscripts/TE_GTF/hg38_rmsk_TE.gtf.gz -O hg38/repeats/hg38.p14_rmsk_TE.gtf.gz
wget https://labshare.cshl.edu/shares/mhammelllab/www-data/TElocal/annotation_tables/hg38_rmsk_TE.gtf.locInd.locations.gz -O hg38/repeats/hg38.p14_rmsk_TE.gtf.locInd.locations.gz
wget https://labshare.cshl.edu/shares/mhammelllab/www-data/TElocal/prebuilt_indices/hg38_rmsk_TE.gtf.locInd.gz -O hg38/repeats/hg38.p14_rmsk_TE.gtf.locInd.gz

# Generate md5 checksums for all downloaded and processed files
md5sum-lite hg38/genome/hg38.p14.ucsc.fa.gz > hg38/genome/hg38.p14.ucsc.fa.gz.md5
md5sum-lite hg38/genome/hg38.p14.chrom.sizes > hg38/genome/hg38.p14.chrom.sizes.md5
md5sum-lite hg38/genes/gencode.v47.hg38.p14.chr_patch_hapl_scaff.annotation.gtf.gz > hg38/genes/gencode.v47.hg38.p14.chr_patch_hapl_scaff.annotation.gtf.gz.md5
md5sum-lite hg38/repeats/hg38.p14_rmsk_TE.gtf.gz > hg38/repeats/hg38.p14_rmsk_TE.gtf.gz.md5
md5sum-lite hg38/repeats/hg38.p14_rmsk_TE.gtf.locInd.locations.gz > hg38/repeats/hg38.p14_rmsk_TE.gtf.locInd.locations.gz.md5
md5sum-lite hg38/repeats/hg38.p14_rmsk_TE.gtf.locInd.gz > hg38/repeats/hg38.p14_rmsk_TE.gtf.locInd.gz.md5
set -e

mkdir -p mm39/genome mm39/genes mm39/repeats

# Download the chromosome sizes from UCSC
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.chrom.sizes -O mm39/genome/mm39.chrom.sizes

# Download the chromosome FASTAs from UCSC
mkdir tmp_mm39; cd tmp_mm39
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.chromFa.tar.gz
tar xvzf mm39.chromFa.tar.gz

# Remove chromosome M and create a unique FASTA file.
# Chromosome will be ordered following the chromosome size ordering.
cut -f1 ../mm39/genome/mm39.chrom.sizes | awk '{print $1".fa"}' | xargs cat | gzip -c > ../mm39/genome/mm39.ucsc.fa.gz

# Clean up
cd ..; rm -r tmp_mm39

# Download the GENCODE complete gene annotations.
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M36/gencode.vM36.chr_patch_hapl_scaff.annotation.gtf.gz -O mm39/genes/gencode.vM36.mm39.p6.chr_patch_hapl_scaff.annotation.gtf.gz

# Download TEtranscript pre-indexed files. The repeat masker is filtered and contains less entries than the UCSC one. 
# tRNA are removed and other simple and short entries. Exact details are unclear.
# Note that the script will remove most simple repetitive sequences and short non-coding RNA (e.g. tRNA)
# All these files need to be uncompressed using the `gzip -d` command. TElocal does not accept compressed files.
wget https://labshare.cshl.edu/shares/mhammelllab/www-data/TEtranscripts/TE_GTF/mm39_rmsk_TE.gtf.gz -O mm39/repeats/mm39_rmsk_TE.gtf.gz
wget https://labshare.cshl.edu/shares/mhammelllab/www-data/TElocal/annotation_tables/mm39_rmsk_TE.gtf.locInd.locations.gz -O mm39/repeats/mm39_rmsk_TE.gtf.locInd.locations.gz
wget https://labshare.cshl.edu/shares/mhammelllab/www-data/TElocal/prebuilt_indices/mm39_rmsk_TE.gtf.locInd.gz -O mm39/repeats/mm39_rmsk_TE.gtf.locInd.gz

# Generate md5 checksums for all downloaded files
md5sum-lite mm39/genome/mm39.ucsc.fa.gz > mm39/genome/mm39.ucsc.fa.gz.md5
md5sum-lite mm39/genome/mm39.chrom.sizes > mm39/genome/mm39.chrom.sizes.md5
md5sum-lite mm39/genes/gencode.vM36.mm39.p6.chr_patch_hapl_scaff.annotation.gtf.gz > mm39/genes/gencode.vM36.mm39.p6.chr_patch_hapl_scaff.annotation.gtf.gz.md5
md5sum-lite mm39/repeats/mm39_rmsk_TE.gtf.gz > mm39/repeats/mm39_rmsk_TE.gtf.gz.md5
md5sum-lite mm39/repeats/mm39_rmsk_TE.gtf.locInd.locations.gz > mm39/repeats/mm39_rmsk_TE.gtf.locInd.locations.gz.md5
md5sum-lite mm39/repeats/mm39_rmsk_TE.gtf.locInd.gz > mm39/repeats/mm39_rmsk_TE.gtf.locInd.gz.md5
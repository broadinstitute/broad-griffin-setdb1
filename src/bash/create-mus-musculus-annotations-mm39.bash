set -e

mkdir -p mm39/genome mm39/genes mm39/repeats

# Download the chromosome sizes from UCSC
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.chrom.sizes -O mm39/genome/mm39.chrom.sizes

# Download the chromosome FASTAs from UCSC
mkdir tmp;cd tmp
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.chromFa.tar.gz
tar xvzf mm39.chromFa.tar.gz

# Remove chromosome M and create a unique FASTA file.
# Chromosome will be ordered following the chromosome size ordering.
cut -f1 ../mm39/genome/mm39.chrom.sizes | awk '{print $1".fa"}' | xargs cat > ../mm39/genome/mm39.ucsc.fa

# Clean up
cd ..; rm -r tmp

# Download the GENCODE complete gene annotations.
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M36/gencode.vM36.chr_patch_hapl_scaff.annotation.gtf.gz -O - | gzip -dc > mm39/genes/gencode.vM36.mm39.p6.chr_patch_hapl_scaff.annotation.gtf

# Download TEtranscript pre-indexed files. The repeat masker is filtered and contains less entries than the UCSC one. 
# tRNA are removed and other simple and short entries. Exact details are unclear.
# Note that the script will remove most simple repetitive sequences and short non-coding RNA (e.g. tRNA)
# All these files need to be uncompressed using the `gzip -d` command. TElocal does not accept compressed files.
wget https://labshare.cshl.edu/shares/mhammelllab/www-data/TEtranscripts/TE_GTF/mm39_rmsk_TE.gtf.gz -O - | gzip -dc > mm39/repeats/mm39_rmsk_TE.gtf
wget https://labshare.cshl.edu/shares/mhammelllab/www-data/TElocal/annotation_tables/mm39_rmsk_TE.gtf.locInd.locations.gz -O - | gzip -dc > mm39/repeats/mm39_rmsk_TE.gtf.locInd.locations
wget https://labshare.cshl.edu/shares/mhammelllab/www-data/TElocal/prebuilt_indices/mm39_rmsk_TE.gtf.locInd.gz -O - | gzip -dc > mm39/repeats/mm39_rmsk_TE.gtf.locInd

# Generate md5 checksums for all downloaded files
md5sum-lite mm39/genome/mm39.ucsc.fa > mm39/genome/mm39.ucsc.fa.md5
md5sum-lite mm39/genome/mm39.chrom.sizes > mm39/genome/mm39.chrom.sizes.md5
md5sum-lite mm39/genes/gencode.vM36.mm39.p6.chr_patch_hapl_scaff.annotation.gtf > mm39/genes/gencode.vM36.mm39.p6.chr_patch_hapl_scaff.annotation.gtf.md5
md5sum-lite mm39/repeats/mm39_rmsk_TE.gtf > mm39/repeats/mm39_rmsk_TE.gtf.md5
md5sum-lite mm39/repeats/mm39_rmsk_TE.gtf.locInd.locations > mm39/repeats/mm39_rmsk_TE.gtf.locInd.locations.md5
md5sum-lite mm39/repeats/mm39_rmsk_TE.gtf.locInd > mm39/repeats/mm39_rmsk_TE.gtf.locInd.md5
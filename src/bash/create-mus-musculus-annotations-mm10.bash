set -e

mkdir -p mm10/genome mm10/genes mm10/repeats

# Download the chromosome sizes from UCSC
wget http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes -O mm10/genome/mm10.ucsc.chrom.sizes

# Download the chromosome FASTAs from UCSC
mkdir tmp; cd tmp
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/chromFa.tar.gz
tar xvzf chromFA.tar.gz

# Remove chromosome M and create a unique FASTA file.
# Chromosome will be ordered following the chromosome size ordering.
cut -f1 ../mm10/genome/mm10.ucsc.chrom.sizes | awk '{print $1".fa"}' | xargs cat > ../mm10/genome/mm10.ucsc.fa

# Clean up
cd ..; rm -r tmp

# Download the GENCODE gene annotations
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M1/gencode.vM1.annotation.gtf.gz -O - | gzip -dc > mm10/genes/gencode.vM1.mm10.annotation.gtf

# Download TEtranscript pre-indexed files. The repeat masker is filtered and contains less entries than the UCSC one. 
# tRNA are removed and other simple and short entries. Exact details are unclear.
# Note that the script will remove most simple repetitive sequences and short non-coding RNA (e.g. tRNA)
# All these files need to be uncompressed using the `gzip -d` command. TElocal does not accept compressed files.
wget https://labshare.cshl.edu/shares/mhammelllab/www-data/TEtranscripts/TE_GTF/mm10_rmsk_TE.gtf.gz -O - | gzip -dc > mm10/repeats/mm10_rmsk_TE.gtf
wget https://labshare.cshl.edu/shares/mhammelllab/www-data/TElocal/annotation_tables/mm10_rmsk_TE.gtf.locInd.locations.gz -O - | gzip -dc > mm10/repeats/mm10_rmsk_TE.gtf.locInd.locations
wget https://labshare.cshl.edu/shares/mhammelllab/www-data/TElocal/prebuilt_indices/mm10_rmsk_TE.gtf.locInd.gz -O - | gzip -dc > mm10/repeats/mm10_rmsk_TE.gtf.locInd

# Generate md5 checksums for all downloaded files
md5sum-lite mm10/genome/mm10.ucsc.fa > mm10/genome/mm10.ucsc.fa.md5
md5sum-lite mm10/genome/mm10.ucsc.chrom.sizes > mm10/genome/mm10.ucsc.chrom.sizes.md5
md5sum-lite mm10/genes/gencode.vM1.mm10.annotation.gtf > mm10/genes/gencode.vM1.mm10.annotation.gtf.md5
md5sum-lite mm10/repeats/mm10_rmsk_TE.gtf > mm10/repeats/mm10_rmsk_TE.gtf.md5
md5sum-lite mm10/repeats/mm10_rmsk_TE.gtf.locInd.locations > mm10/repeats/mm10_rmsk_TE.gtf.locInd.locations.md5
md5sum-lite mm10/repeats/mm10_rmsk_TE.gtf.locInd > mm10/repeats/mm10_rmsk_TE.gtf.locInd.md5
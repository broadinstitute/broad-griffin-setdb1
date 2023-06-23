set -e


mkdir -p mm10/genome mm10/genes mm10/repeats
# Download the chromosome FASTAs from UCSC
mkdir tmp
cd tmp
wget http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/chromFa.tar.gz
tar xvzf chromFA.tar.gz

# Download the chromosome sizes from UCSC
wget http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes

# Remove chromosome M and create a unique FASTA file.
# Chromosome will be ordered following the chromosome size ordering.
grep -v chrM mm10.chrom.sizes > mm10.ucsc.chrom.sizes.x_mito
cut -f1 mm10.ucsc.chrom.sizes.x_mito | awk '{print $1".fa"}' | xargs cat > mm10.ucsc.x_mito.fa

# Move files in final location
mv mm10.ucsc.x_mito.fa mm10.ucsc.chrom.sizes.x_mito ../mm10/genome

# Clean up
cd ..
rm -r tmp

mkdir tmp; cd tmp
# Download the refseq gene annotations from UCSC and convert them to GTF.
# Following these instructions
# https://genomewiki.ucsc.edu/index.php/Genes_in_gtf_or_gff_format
wget http://hgdownload.cse.ucsc.edu/goldenPath/mm10/database/refGene.txt.gz
gzip -d refGene.txt.gz
cut -f 2- refGene.txt > refGene.input
../bins/genePredToGtf file refGene.input mm10.refGene.tmp
sed 's/refGene.input/mm10_ucsc_refGene/g' mm10.refGene.tmp > mm10.refGene.tmp.fixed
sort -k1,1 -k4,4n mm10.refGene.tmp.fixed > mm10.refGene.gtf

mv mm10.refGene.gtf ../mm10/genes
cd ..
rm -r tmp

# Download the reeat masker file from UCSC
#wget http://hgdownload.cse.ucsc.edu/goldenPath/mm10/database/rmsk.txt.gz
#mv rmsk.txt.gz mm10/repeats/mm10.ucsc.rmsk.txt.gz


# Download TEtranscript pre-indexed files. The repeat masker is filtered and contains less entries than the UCSC one. 
# tRNA are removed and other simple and short entries. Exact details are unclear.
# Note that the script will remove most simple repetitive sequences and short non-coding RNA (e.g. tRNA)
# All these files need to be uncompressed using the `gzip -d` command. TElocal does not accept compressed files.
 wget https://labshare.cshl.edu/shares/mhammelllab/www-data/TElocal/annotation_tables/mm10_rmsk_TE.gtf.locInd.locations.gz -O mm10/repeats/mm10_rmsk_TE.gtf.gz
 #mv mm10_rmsk_TE.gtf.locInd.locations.gz ../mm10/repeats
 wget https://labshare.cshl.edu/shares/mhammelllab/www-data/TEtranscripts/TE_GTF/mm10_rmsk_TE.gtf.gz -O mm10/repeats/mm10_rmsk_TE.gtf.gz
 #mv  mm10_rmsk_TE.gtf.gz ../mm10/repeats
 wget https://labshare.cshl.edu/shares/mhammelllab/www-data/TEtranscripts/prebuilt_indices/mm10_rmsk_TE.gtf.ind.gz -O mm10/repeats/mm10_rmsk_TE.gtf.ind.gz
 #mv mm10_rmsk_TE.gtf.ind.gz ../mm10/repeats
 wget https://labshare.cshl.edu/shares/mhammelllab/www-data/TElocal/prebuilt_indices/mm10_rmsk_TE.gtf.locInd.gz -O mm10/repeats/mm10_rmsk_TE.gtf.locInd.gz
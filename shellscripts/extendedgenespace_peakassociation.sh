#!/bin/bash
# This script takes the gasAcu1 gtf we used for the analysis, generates a bed file from it and extends gene intervals to include  the intergenic space around them (so that neighbouring intervals will overlap)
# It then uses bedtools intersect to find which ATAC-seq peaks are within each extended gene interval

# define filepaths
# note: resulting files will be saved to current working directory
gtf=/path/to/projectfolder/Gasterosteus_aculeatus.BROADS1.104.gasAcu1.lifted.gtf
base_name=${gtf%.*}
chromsize_file=/path/to/projectfolder/gasAcu1.sizes

# extract relevant columns from gtf for each gene entry
awk '{if ($3 == "gene") {print $1 "\t" $4 "\t" $5 "\t" $10} }' ${gtf} > ${basename}.genes.bed
#awk '{if ($3 < $2){print $0}}' ${basename}.genes.bed

# remove special characters carried over from the gtf
#Note: sed options specified like this to run on MACOS!
sed -i '' 's/"//g' ${basename}.genes.bed 
sed -i '' 's/;//g' ${basename}.genes.bed

sort -k1,1 -k2,2n ${basename}.genes.bed > ${basename}.genes.sorted.bed

# for each gene, find closest upstream and downstream gene
bedtools closest -io -id -D ref -t last -a ${basename}.genes.sorted.bed -b ${basename}.genes.sorted.bed > upstream.bed
bedtools closest -io -iu -D ref -t first -a ${basename}.genes.sorted.bed -b ${basename}.genes.sorted.bed > downstream.bed

# extend each gene's interval so that start is end of previous gene + 1 or start of chromosome, and stop is start of next gene - 1 or tag it
paste ${basename}.genes.sorted.bed upstream.bed downstream.bed | awk 'BEGIN {OFS="\t"} {
    if ($9 == "."){
        $2 = 1
    } else {
        $2 = ($11 + 1)
    }; #set start of gene to stop of previous gene + 1, or 1 if no gene upstream
    if ($18 == "."){$3 = $3 "***" # tag genes without downstream neighbour
    } else {
        $3 = ($19 - 1);
    }
    print $1, $2, $3, $4;
}' > extended.bed

# extend genes with tagged stop to length of chromosome
awk 'BEGIN {OFS="\t"}
    NR==FNR {sizes[$1]=$2; next}
    {
        if ($3 ~ /\*\*\*$/) {
            gsub(/\*\*\*$/, "", $3)
            if ($1 in sizes) $3 = sizes[$1]
        }
        print
    }' ${chromsize_file} extended.bed > ${basename}.genes_extended.bed



# do some checks
diff -s <(bedtools closest -N -fd -D ref -a ${basename}.genes.sorted.bed -b ${basename}.genes.sorted.bed) <(bedtools closest -io -fd -D ref -a ${basename}.genes.sorted.bed -b ${basename}.genes.sorted.bed)

# intersect peaks with extended gene intervals
bedtools intersect -wa -wb -a ${basename}.genes_extended.bed -b peaks_filtered.bed > genes_peaks.intersect
awk '{print $4, $8}' genes_peaks.intersect > genes_peaks.intersect.list
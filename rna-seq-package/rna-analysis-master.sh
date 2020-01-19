#!/usr/bin/env bash

#https://github.com/OpenGene/fastp
#http://labshare.cshl.edu/shares/gingeraslab/www-data/dobin/STAR/STAR.posix/doc/STARmanual.pdf

#NOTE: create a new directory for each instance this script is run to avoid overwriting

#setting variables to hold input reference and annotation, both must be unzipped
ref=$(grep "ref.celegan>" ./SAMPLE_SHEET.csv | cut -d "," -f 2)
gtf=$(grep "gtf.celegan>" ./SAMPLE_SHEET.csv | cut -d "," -f 2)

#making directories to hold STAR genome during STAR genome generation
mkdir ./STAR_genomeDIR

#setting variables to hold genome location
STAR_genomeDIR="./STAR_genomeDIR"

#fetching in directory specified in SAMPLE_SHEET.csv
in_dir=$(grep "IN.DIR>" ./SAMPLE_SHEET.csv | cut -d "," -f 2)
#fetching output directory specified in SAMPLE_SHEET.csv
out_dir=$(grep "OUT.DIR>" ./SAMPLE_SHEET.csv | cut -d "," -f 2)
#isolated files from SAMPLE_SHEET.txt
files=$(awk '/[0-9]+.r[1-2]>,[A-Za-z0-9_-]+*/{print $0}' ./SAMPLE_SHEET.csv | cut -d "," -f 2)
#isolated aliases from SAMPLE_SHEET.txt
grab_aliases=$(awk '/[0-9]+.alias>,[A-Za-z0-9_-]+*/{print $0}' ./SAMPLE_SHEET.csv | cut -d "," -f 2)

#file extension to be added to outputted filtered/trimmed fastq files
ext="QC-OUT.fastq.gz"

#appending the aliases into a list to index later in for loop
aliases=()
for a in $grab_aliases; do
  aliases+=("$a")
done

echo -e "Directory containing input FASTQ file(s):" $in_dir "\n"

echo "Files submitted to workflow:"
for id in $files; do
  echo $id
done

echo -e "\nCorresponding aliases submitted to workflow:"
for ali in $grab_aliases; do
  echo $ali
done

count=0
alias_count=0
for file in $files; do
  count=$((count+1))
  if [[ $((count%2)) == 1 ]]; then
    in_r1=$file
  else
    [[ $((count%2)) == 0 ]]
    echo -e "\n---------Start of analysis for next read pair---------\n"
    in_r2=$file
    out_r1=$(echo $in_r1 | cut -d'.' -f 1)_$ext
    out_r2=$(echo $in_r2 | cut -d'.' -f 1)_$ext
    id=${aliases[$alias_count]}
    count_file=$id"_"$count"_counts.tsv"
    html_id=$id"_REPORT.html"
    json_id=$id"_REPORT.json"
    echo -e "Beginning QC, quality filtering, and adapter trimming for samples:"
    echo $in_r1
    echo -e $in_r2"\n"
    fastp -i $in_dir$in_r1 -I $in_dir$in_r2 -o $out_dir$out_r1 -O $out_dir$out_r2 \
    --html $html_id --json $json_id --report_title "Nemametrix_QC_Report_"$id \
    --average_qual 30 --length_required 20 \
    --detect_adapter_for_pe \
    --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
    --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    --trim_poly_g \
    --thread 7
    #--dont_overwrite
    echo -e "\nBeginning alignment and feature counting for samples:"
    echo -e "\nAlias id:" $id"\n"
    echo $out_r1
    echo -e $out_r2"\n"
    STAR --runThreadN 7 --runMode genomeGenerate \
    --genomeDir $STAR_genomeDIR \
    --genomeFastaFiles $ref
    STAR --runThreadN 7 --runMode alignReads \
    --readFilesCommand zcat \
    --readFilesIn $out_dir$out_r1 $out_dir$out_r2 \
    --genomeDir $STAR_genomeDIR \
    --outFileNamePrefix $id"_"$count"_" \
    --sjdbGTFfile $gtf \
    --outStd SAM \
    | samtools sort - -o - \
    | samtools view - -o - \
    | htseq-count --stranded no --type exon --idattr gene_id --mode union --order pos \
    - $gtf > $count_file
    echo -e "feature\t"$id > tmp.txt
    cat tmp.txt $count_file > rename.txt && mv rename.txt $count_file
    alias_count=$((alias_count+1))
  fi
done

#organizing .json and .html fastp qc reports into a single folder
mkdir QC_REPORTS
mv *_REPORT.* ./QC_REPORTS/

#removing remaining tmp.txt to avoid confusion as its utility was only needed in loop
rm tmp.txt

#combining all individual counts files ending with counts.tsv, ensure no other files with that ending are in the same dir
paste *counts.tsv | awk 'BEGIN {OFS="\t"; FS="\t"}; {j=$1; for (i=2;i<=NF;i+=2) {j=j FS $i} print j}' > "merged_counts.tsv"

echo -e "\n------------QC, mapping, and feature counting complete.------------"

echo -e "\n------------Running differential gene expression analysis in R and knitting results to html file------------"

#creates html report
Rscript -e 'library(rmarkdown); rmarkdown::render("./dge.Rmd", "html_document")'

#opens R shiny application
R -e "shiny::runApp('shiny', launch.browser=TRUE)"

echo -e "\n*** Analysis Complete ***"

echo -e "\nTo relaunch the R shiny dashboard, run the command: R -e 'shiny::runApp('shiny', launch.browser=TRUE)'"

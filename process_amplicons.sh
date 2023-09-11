#!/bin/bash

#read in directories
bioawk_dir=/usr/local/bin/bioawk
read_dir=/PATH/TO/FASTQ
out_dir=/PATH/FOR/PROCESSED/READS
blast_dir=/PATH/FOR/BEST/BLAST

#read in BLAST database
blast_db=/PATH/TO/macro.noura_blast_db.fasta

#run each file at a time
for f in /Users/katjakasimatis/Documents/UofT/experiments/Rebecca_ReproductiveInterference/RS_amplicon_round3/00_fastq/*.fastq
do
	base=`basename $f _001.fastq`
	
	#create table-delimited table of read names and lengths
	${bioawk_dir}/bioawk -cfastx '{print $name, length($seq)}' $f > reads/${base}_LENGTHS.txt
	
	#filter reads by mean quality score (Phred Score > 30) and length (>50bp)
	#${bioawk_dir}/bioawk -c fastx '{ if(meanqual($qual) > 30 && length($seq) > 50) { print "@"$name; print $seq; print "+"; print $qual; }}' $f > ${out_dir}/${base}_FILTERED.fastq

	#convert FASTQ file to FASTA file format
	#${bioawk_dir}/bioawk -c fastx '{print ">"$name; print $seq}' ${out_dir}/${base}_FILTERED.fastq > ${out_dir}/${base}_FILTERED.fasta
	
	#BLAST to amplicon sequences
	#blastn -subject ${blast_db}/macro.noura_blast_db.fasta -query ${out_dir}/${base}_FILTERED.fasta -outfmt 6 -out ${blast_dir}/${base}.txt
done

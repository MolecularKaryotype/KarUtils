wget -O - "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/gencode.v40.annotation.gtf.gz" |\
gunzip -c | grep 'transcript_type "protein_coding"' |\
awk '($3=="exon") {printf("%s\t%s\t%s\n",$1,int($4)-1,$5);}' |\
sort -T . -t $'\t' -k1,1 -k2,2n | bedtools merge
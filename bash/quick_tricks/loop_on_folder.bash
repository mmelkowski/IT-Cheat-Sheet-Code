# To loop on file on a folder:
ls *.fastq.gz | while read fq; do echo $fq; done

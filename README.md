# Count-Reads-aligned-to-non-overlapping-bins-in-the-genome
This python code will count reads aligned to non-overlapping bins in the genome. It will retrive genome information from the input bam file. 
It's output will be a bedgraph file which can be converted to a bigwig file using the bedGraphToBigWig utitlity tool from UCSC (http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/). The bigwig file can be used for visulalization in tools like IGV or uploaded to UCSC genome browser. It will also print a genome wide line plot of the read counts for each chromosome per bin.
Suce a program is useful when working with NGS data. For example, when working with ChIP-Seq data you might want to see if the randomly sonicated chromatin control (Input) contain any large copy number alterations which might cause bias during peak calling.

The program as an inbuild argument parser which will print the help massage if no arguments are passed to it or when the -h option is used. 

Example Usage:
1. Mandatory fields, this will print the raw read counts per bin 
  python countReadsPerNonOverlappingBins.py --bam in.bam --binSize 1000 --outFile out.bw 
2. With Reads Per Million normalization
  python countReadsPerNonOverlappingBins.py --bam in.bam --binSize 1000 --outFile out.bw --RPM True
3. Plot normalized genome wide line graph
   python countReadsPerNonOverlappingBins.py --bam in.bam --binSize 1000 --outFile out.bw --RPM True --makePlot True

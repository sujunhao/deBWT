# deBWT

Implementation of [deBWT](http://bioinformatics.oxfordjournals.org/content/32/12/i174.long) (de Bruijn branch-based BWT constructor ), which is a parallelizable [BWT](http://www.cs.jhu.edu/~langmea/resources/bwt_fm.pdf) constructing algorithm

deBWT uses the de Bruijn graph to reduce the time consumed in the repeating regions when constructing BWT

step:
1. deBWT counts the kmers existing on the target string. 
2. uses kmers to construct the de Bruijn graph implicitly.
3. uses the graph to construct BWT.


> usage

```
#should specify the target string and kmer length
deBWT.sh file.fa kmer_len


#example
>./deBWT.sh E.coli.fa 20
Kmer counting time use  : 7 's
g++ deBWT.cpp -o deBWT
BWT construct time use  : 7 's


#the output bwt file  is *.bwt in current dir
```

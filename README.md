# deBWT

[deBWT](http://www.ncbi.nlm.nih.gov/pubmed/27307614) (de Bruijn branch-based BWT constructor ) is a parallelizable [BWT](http://www.cs.jhu.edu/~langmea/resources/bwt_fm.pdf) constructing algorithm

deBWT uses de Bruijn graph to reduce the time consumed in the repeating regions when constructing BWT

Firstly, deBWT counts the kmers existed on target string. 

Then it uses kmers to construct the de Bruijn implicitly.

Then it uses the graph to construct BWT

usage

```
#should specify the target string and kmer length
deBWT.sh file.fa kmer_len


#example
>./deBWT.sh E.coli.fa 20
mer_counts time use   : 8 's
g++ deBWT.cpp -o deBWT
bwt construct time use: 7 's

#the output bwt file  is *.bwt in current dir
```

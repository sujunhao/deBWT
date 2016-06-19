# deBWT

[deBWT](http://www.ncbi.nlm.nih.gov/pubmed/27307614) (de Bruijn branch-based BWT constructor ) is a parallelizable [BWT](http://www.cs.jhu.edu/~langmea/resources/bwt_fm.pdf) constructing algorithm

deBWT use de Bruijn graph to reduce the time consumed in repeat region when construct BWT
Firstly, deBWT count the kmers exist in target string. Then it use kmers to construct the de Bruijn implicitly.
Then it uses the graph to construct BWT

usage

```
#should specific the target string and kmer length
deBWT.sh file.fa kmer_len


#example
>./deBWT.sh E.coli.fa 20
mer_counts time use   : 8 's
g++ deBWT.cpp -o deBWT
bwt construct time use: 7 's

#the output bwt file  is *.bwt in current dir
```
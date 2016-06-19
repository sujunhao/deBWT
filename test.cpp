#include <cstdio>
#include <iostream>
#include <cstring>
#include <string>
#include <algorithm>
#include <cstdlib>
#include <bitset>
#include <fstream>
#include <unistd.h>
#include <stdint.h>
#define PX(X) std::cout << (X) << std::endl
using namespace std;

// #define PRINTLOG
size_t KMER=20;
size_t KMER_N=4864385;

int main(int argc, char** argv) 
{

    int c;
    opterr = 0;
    char i_f[100]="E.coli.fa";
    char j_f[100]="mer_counts_dumps.fa";
    //set argument,p for PointerListLen,w for WindowListLen,f for dna_read file(if the -pwf is provided, argument is required)
    while ((c = getopt (argc, argv, "k:d:m:n:h")) != -1) 
    {
        switch (c)
        {
            case 'k':
                if (optarg)
                    KMER = atol(optarg);
                break;
            case 'd':
                if (optarg)
                    strcpy(i_f, optarg);
                break;
            case 'm':
                if (optarg)
                    strcpy(j_f, optarg);
                break;
            case 'n':
                if (optarg)
                    KMER_N = atol(optarg);
                break;
            case 'h':
                printf("\tuse dna file to create a RHT table to file out_RHT\n");
                printf("\tcreate_RHT [option]\n");
                printf("\t-k \t  kmer len            [%lu]\n", (unsigned long)KMER);
                printf("\t-d \t  dna file            [%s]\n", i_f);
                printf("\t-m \t  mer_counts_file     [%s]\n", i_f);
                printf("\t-n \t  mer_counts_file len [%lu]\n", (unsigned long)KMER_N);
                printf("\t-h \t  help\n");
                return 0;
                break;
            default:
                abort ();
        }
    }



    std::cout << "kmer: " << KMER << endl
              << "kmer_n: " << KMER_N << endl
              << "| dna file: " << string(i_f) << endl;


    

    return 0;
    
}

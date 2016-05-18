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

uint64_t get_c[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                    0, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                    3};


size_t KMER=20;
size_t KMER_N=4864385;
uint64_t K_mask = 0xffffffff;

int main()
{
	char i_f[100]="E.coli.fa";
	char o_f[100]="out";
    char j_f[100]="mer_counts_dumps.fa";
    ifstream inRef;
    ofstream out;
    inRef.open(i_f);
    out.open(o_f);

    FILE *inJF;
    inJF = fopen(j_f, "r");



    //-----------------------------------get dna_ref string
    string dna_name, dna_s, dna_f;
    getline(inRef, dna_name);
    while (inRef >> dna_s) dna_f.append(dna_s);
    size_t dna_z=dna_f.size();


    //-----------------------------------k_mer counting get from jellyfish
    uint64_t *K2c;
    size_t kmer=KMER, kmer2=KMER_N;

    K2c = new uint64_t[kmer2];
    memset(K2c, 0, sizeof(K2c));

    int tn;
    // printf("asd\n");
    char tkmer[100];
    for (size_t i=0; i<KMER_N; ++i)
    {

        // fgets(tkmer, 100, inJF);
        fscanf(inJF, ">%d%s\n", &tn, tkmer);
        {printf(">%d\n", tn); puts(tkmer);}

    }

 //    size_t kmer=KMER, kmer2=((KMER+2)<<1);
 //    PX(kmer2);
 //    uint64_t *K2, *K2c;


 //    K2 = new uint64_t[(size_t)1<<kmer2];
 //    size_t K2_index=0;

	// uint64_t tmp = 0, tar;
 //    uint64_t MASK = -1;
 //    PX(std::bitset<64>(MASK));
 //    //init mask
 //    MASK = MASK >> (64 - kmer2);

 //    PX(kmer2);
 //    PX(std::bitset<64>(MASK));
 //    //test k_mer max hitted time
 //    size_t mx=0;
 //    for (size_t i=0; i<dna_z; ++i)
 //    {
 //    	tmp = (tmp << 2) | get_c[dna_f[i]];
 //        tar = tmp & MASK;
 //        if (i>=kmer+1)
 //        {
 //        	// PX(std::bitset<64>(tar));
 //        	if (K2c[tar] == 0) K2[K2_index++] = tar;
 //        	++K2c[tar];
 //            mx = max(mx, (size_t)K2c[tar]);
 //        }
 //    }
 //    PX(1<<kmer2);
 //    cout << K2_index << endl;
 //    PX(mx);

    // //-----------------------------------sort k_mer 
    // // PX(std::bitset<32>(MASK));
    // K_mask = MASK >> 2;
    // sort(K2, K2+K2_index, cmp);
    // for (size_t i=0; i<100; ++i) cout << std::bitset<32>(K2[i]) << endl;

    //-----------------------------------for each k_mer, mark multip in or out if it is
    // bool in_, out_;

    // uint64_t k_mid, k_left, k_right;
    // uint64_t t_mask = 0xffffffff, t1_mask = t_mask, t2_mask = t_mask;
    // t_mask = t_mask >> (64 - (kmer<<1));
    // t_mask = t_mask << 2;
    // t1_mask = t1_mask >> (64 - 2);
    // t1_mask = t1_mask << ((kmer << 1) + 2);
    // t2_mask = t2_mask >> (64 - 2);

    // uint64_t *K;
    // K = new uint64_t[1<<(kmer<<1)];
    // memset(K, 0, sizeof(K));
    // size_t BC_index = 0, count_bc;
    // int out_sp = 63, in_sp = 32;
    // for (size_t i=0, j; i<K2_index; )
    // {
    //     k_mid = K2[i] & t_mask;
    //     k_left = K2[i] & t1_mask;
    //     k_right = K2[i] & t2_mask;
    //     in_ = 0;
    //     out_ = 0;
    //     count_bc = 0;
    //     count_bc += K2c[K2[i]];
    //     j = i+1;
    //     // if (i<5000) {PX(std::bitset<32>(K2[i])); PX(K2c[K2[i]]);}
    //     while ((K2[j]&t_mask) == k_mid && j<K2_index)
    //     {
    //         // if (i<5000) {PX(std::bitset<32>(K2[j])); PX(K2c[K2[j]]);}
    //         count_bc += K2c[K2[j]];
    //         if (in_ && out_) ++j;
    //         else
    //         {
    //             if (k_left != (K2[j]&t1_mask)) in_ = 1;
    //             if (k_right != (K2[j]&t2_mask)) out_ = 1;
    //             ++j;
    //         }
    //     }
    //     //check multip in 
    //     //check multip out
    //     //if multip in 
    //     //write multip in first index
    //     k_mid = k_mid >> 2;
    //     if (out_)
    //     {
    //         //if multip out mark it
    //         K[k_mid] |= ((uint64_t)1<<out_sp);
            
    //     }
    //     if (in_) 
    //     {
    //         //if multip in write multip in first index and multip in number
    //         K[k_mid] |= ((uint64_t)(count_bc)<<in_sp);
    //         K[k_mid] |= (BC_index);
    //         BC_index += count_bc;
    //     }

    //     // PX("---");
    //     // if (i<5000) PX(std::bitset<64>(K[k_mid]));
    //     i = j;
    // }


    //-----------------------------------for each k_mer in dna string, use binary search to find k_mer in K
    //-----------------------------------if is multip out, construct branch code (index++)
    //-----------------------------------if is multip in, store the branch code index(**** index should +1 ****)


  	//-----------------------------------use K2 and K2c to construct FM-index
    //-----------------------------------handel the last k-1 k_mer (ATG$, TG$, G$)
    //-----------------------------------insert in BWT, concern as case 1 

	// delete [] K2;
    delete [] K2c;
    inRef.close();
    out.close();
    fclose(inJF);

    printf("Time used = %.2f\n",  (double)clock() / CLOCKS_PER_SEC);
    return 0;
}

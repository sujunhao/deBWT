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

typedef struct kmernote {
    uint64_t km;
    uint32_t in;
    uint32_t index;
}KMER_NODE;


bool cmp(const uint64_t a, const uint64_t b)
{
    // return a<b;
    return (a<<2) < (b<<2);
}
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
    inJF = fopen(j_f, "rb");



    //-----------------------------------get dna_ref string
    string dna_name, dna_s, dna_f;
    getline(inRef, dna_name);
    while (inRef >> dna_s) dna_f.append(dna_s);
    size_t dna_z=dna_f.size();


    //-----------------------------------k_mer counting get from jellyfish
    uint64_t *K2;
    size_t kmer=KMER, kmer2=KMER+2, kmer2len=KMER_N;

    K2 = new uint64_t[kmer2len];

    int tn;
    char tkmer[100];

    uint64_t tmp = 0, tar;
    uint64_t MASK2 = -1, k2clen=10;
    MASK2 = MASK2 >> (64 - kmer2*2);

    for (size_t i=0; i<kmer2len; ++i)
    {

        fscanf(inJF, ">%d\n%s\n", &tn, tkmer);
        tmp=0;
        for (size_t j=0; j<kmer+2; ++j) tmp = (tmp << 2) | get_c[tkmer[j]];
        tar = tmp << (64 - kmer2*2);
        tar = tar|((uint64_t)tn);

        K2[i]=tar;
        // if (i<100){
        //     printf("%d %s ", tn, tkmer);
            // cout << std::bitset<64>(tar) << endl;
            
        // }
    }

    //-----------------------------------sort kmer+2
    sort(K2, K2+kmer2len, cmp);
    


    //-----------------------------------use center kmer to mark multip in & out
    //-----------------------------------mark multip [1]out muitip in number[31] and mutip in index[32]


    uint64_t *io_info;
    io_info = new uint64_t[kmer2len];

    uint64_t thekmer;
    //the center kmer is ((thekmer << 2) >> (64 - kmer*2));
    //first 2 is (thekmer >> 62)
    //the last 2 is (thekmer << 2*(kmer+1) >> 62)
    //the number is (thekmer << 2*(kmer+2) >> 64-2*(kmer+2))

    //for splite the center, left ,right char and number;
    uint64_t mask_c=-1, mask_l=-1, mask_r=-1, mask_n=-1;
    mask_c = mask_c << 2 >> 2 >> (64 - 2*kmer2 + 2) << (64 - kmer2*2 + 2);
    mask_l = mask_l >> 62 << 62;
    mask_r = mask_r >> (64 - 2*kmer2) << (64 - 2*kmer2) << (kmer+1)*2 >> (kmer+1)*2;
    mask_n = mask_n << (2*kmer2) >> (2*kmer2);

    uint64_t mask_out=-1, mask_in=-1, mask_index=-1;
    mask_out = mask_out >> 63 << 63;
    mask_in = mask_in >> 32 << 32 << 1 >> 1;
    mask_index = mask_index << 32 >> 32;

    bool is_in=false, is_out=false;
    // cout << std::bitset<64>(mask_out) << endl << std::bitset<64>(mask_in) << endl << 
    // std::bitset<64>(mask_index)<< endl;

    size_t kmer_index = 0, theindex = 0;
    uint64_t tmp_mask, tmp_num;
    for (size_t i=0, j, k; i<kmer2len;)
    {
        is_in=false; is_out=false; tmp_mask=0;
        j = i+1;
        while (j<kmer2len && (K2[j]&mask_c == K2[i]&mask_c))
        {
            if (K2[j]&mask_l != K2[i]&mask_l) is_out = true;
            if (K2[j]&mask_r != K2[i]&mask_r) is_in = true;
            if (is_in || is_out) printf("asd\n");
            ++j;
        }
        if (is_in || is_out)
        {
            cout << "Bin\n";
            if (is_out) tmp_mask = tmp_mask|mask_out;
            if (is_in)
            {
                tmp_num=0;
                for (k = i; k<j; ++k)
                {
                    tmp_num += (K2[k]&mask_n);
                }
                tmp_mask = tmp_mask|(tmp_num<<32);
                tmp_mask = tmp_mask|theindex;
                theindex += tmp_num;
            }
            io_info[kmer_index++] = tmp_mask;
            if (kmer_index<100){
                cout << "Asd";

                cout << std::bitset<64>(io_info[kmer_index-1]) << endl << std::bitset<64>(io_info[kmer_index-1]<<1>>1>>32) << endl
                << (io_info[kmer_index-1]<<1>>1>>32) << endl;
                
            }
        }
        i = j;

    }

    //-----------------------------------mark multip in and out

    // for (size_t i=0; i<kmer2len; ++i)
    // {
    //     //-----------------------------------splite the center kmer
    //     tkmer = K2[i];
    //     tkmer << 2;
    //     tkmer >> (64 - kmer*2);

    //     //-----------------------------------check multip out
    //     //-----------------------------------mark the multip in number
    //     //-----------------------------------assign the index

    // }

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

	delete [] K2;
    inRef.close();
    out.close();
    fclose(inJF);

    printf("Time used = %.2f\n",  (double)clock() / CLOCKS_PER_SEC);
    return 0;
}

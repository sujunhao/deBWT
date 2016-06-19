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


string dna_name, dna_s, dna_f;

size_t KMER=20;
size_t KMER_N=4864385;
//ie. kmer len 20, kmer2 22, kmer_l, 20*2 kmer2_l 22*2, k2len the K2 len, klen, the K len
size_t kmer, kmer2, kmer_l, kmer2_l, k2len, klen;
uint64_t *K;
string bc;
size_t bc_index;

bool cmp(const uint64_t a, const uint64_t b)
{
    // return a<b;
    return ((a<<2) < (b<<2));
}

bool bwt_cmp(const uint64_t a, const uint64_t b)
{
    uint64_t ta=(a>>3), tb=(b>>3);
    while (ta<bc_index && tb<bc_index && bc[ta] == bc[tb])
    {
        ++ta;
        ++tb;
    }
    if (ta==bc_index) return true;
    if (tb==bc_index) return false;
    return (bc[ta] < bc[tb]);
}

size_t search_k(uint64_t n)
{
    //binary search
    size_t l=0, r=klen-1, m;
    while (l < r)
    {
        m = (l+r)/2;
        if (n<= K[m]) r = m;
        else l = m + 1;
    }
    if (K[r] == n) return r;
    return -1;
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
    getline(inRef, dna_name);
    while (inRef >> dna_s) dna_f.append(dna_s);
    size_t dna_z=dna_f.size();


    //-----------------------------------k_mer counting get from jellyfish

    kmer=KMER;
    kmer2=KMER+2;
    kmer_l=kmer*2;
    kmer2_l=kmer2*2;
    k2len=KMER_N;

    uint64_t *K2;
    K2 = new uint64_t[k2len];

    int tn;
    char tkmer[100];

    uint64_t tmp = 0, tar;

    for (size_t i=0; i<k2len; ++i)
    {

        fscanf(inJF, ">%d\n%s\n", &tn, tkmer);
        tmp=0;
        for (size_t j=0; j<kmer2; ++j) tmp = (tmp << 2) | get_c[tkmer[j]];
        tar = tmp << (64 - kmer2_l);
        tar = tar|((uint64_t)tn);

        K2[i]=tar;
        // if (i<3000){
        //     printf("%d %s ", tn, tkmer);
        //     cout << std::bitset<64>(tar) << endl;
            
        // }
    }

    //-----------------------------------sort kmer+2
    sort(K2, K2+k2len, cmp);
    
    // for (size_t i = 0; i < 50000; ++i)
    // {
    //     cout << std::bitset<64>(K2[i]) << endl;
    //     /* code */
    // }

    // return 0; 

    //-----------------------------------use center kmer to mark multip in & out info.
    //-----------------------------------mark multip [1]out muitip in number[31] and mutip in index[32]


    uint64_t *io_info;
    io_info = new uint64_t[k2len];
    K = new uint64_t[k2len];

    // uint64_t thekmer;
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

    //to write in io_info
    uint64_t mask_out=-1, mask_in=-1, mask_index=-1;
    mask_out = mask_out >> 63 << 63;
    mask_in = mask_in >> 32 << 32 << 1 >> 1;
    mask_index = mask_index << 32 >> 32;

    bool is_in=false, is_out=false;
    // cout << std::bitset<64>(mask_out) << endl << std::bitset<64>(mask_in) << endl;
    // std::bitset<64>(mask_index)<< endl;

    size_t theindex = 0;    //for the BCN store index
    klen = 0;               //for the K[] and io_info
    uint64_t tmp_mask, tmp_num;

    //get the first k_mer
    for (size_t i=0; i<kmer; ++i)
    {
        tmp = (tmp << 2) | get_c[dna_f[i]];
    }
    uint64_t begin_kmer = tmp << (64 - kmer_l) >> 2, last_kmer=0, now_kmer=0;

    // cout << std::bitset<64>(begin_kmer) << endl;
    for (size_t i=0, j, k; i<k2len;)
    {
        // if (i>100000) break;
        is_in=false; is_out=false; tmp_mask=0;
        tmp_num=0;

        j = i+1;
        now_kmer = K2[i]&mask_c;

        //insert the first kmer into proper site (begin_kmer>last_kmer for begin kmer may = 0)
        if (now_kmer >= begin_kmer && begin_kmer >= last_kmer) 
        {
            if (now_kmer == begin_kmer)
            {
                is_in = true;
                if ((K2[i]&mask_r)>>(64-kmer2_l) != get_c[dna_f[kmer]]) is_out = true;
                ++tmp_num;
            }
            //else if the begin_kmer no appear in K2, no need to create the io_info of begin_kmer
        }
        last_kmer = now_kmer;
        //if kmer exit in Kmer+2 then check is_out & is_in
        while (j<k2len && ((K2[j]&mask_c) == (now_kmer)))
        {
            if ((K2[j]&mask_l) != (K2[i]&mask_l)) is_in = true;
            if ((K2[j]&mask_r) != (K2[i]&mask_r)) is_out = true;
            // if (is_in || is_out) printf("asd\n");
            ++j;
        }
        // cout << "----\n";
        // for (k = i; k<j; ++k)
        // {
        //     cout << std::bitset<64>(K2[k]) << endl;
        // }
        if (is_in || is_out)
        {
            // cout << "Bin\n";
            if (is_out) tmp_mask = tmp_mask|mask_out;
            if (is_in)
            {
                for (k = i; k<j; ++k)
                {
                    tmp_num += (K2[k]&mask_n);
                }
                tmp_mask = tmp_mask|theindex;
                tmp_mask = tmp_mask|(tmp_num<<32);
                theindex += tmp_num;
            }
            io_info[klen] = tmp_mask;
            K[klen] = (K2[i]&mask_c) >> (64-kmer2_l+2);
            ++klen;
            // {
                // cout << "*****\n";

                // cout << std::bitset<64>(K[klen-1]) << endl << "i: " << std::bitset<64>(io_info[klen-1]) << endl;
                // << (io_info[klen-1]<<1>>1>>32) << endl;
                
            // }
        }
        // cout << "--\n";
        i = j;
    }
    //-----------------------------------for each k_mer in dna string, use binary search to find k_mer in K
    //-----------------------------------if is multip out, construct branch code (index++)
    //-----------------------------------if is multip in, store the branch code index(**** index should +1 ****)
    uint64_t *BCN;
    size_t BCN_size = (theindex);
    BCN = new uint64_t[BCN_size];
    for (size_t i=0; i<BCN_size; ++i) BCN[i]=0;

    // tmp=0;
    uint64_t mask_k = -1;
    mask_k = mask_k << (64-kmer_l) >> (64-kmer_l);
    // cout << std::bitset<64>(mask_k) << endl;
    bc = dna_f;
    bc_index=1;
    uint64_t BCN_index, BCN_len, tk, tmp_i;
    uint64_t sum = 0;
    for (size_t i=0, index; i<dna_z-1; ++i)
    {
        tmp = (tmp << 2) | get_c[dna_f[i]];
        if (i>=kmer-1)
        {
            tar = tmp&mask_k;
            index = search_k(tar);
            // if (i<100)
            // {
            //     cout << std::bitset<64>(tar) << endl;
            //     // cout << std::bitset<64>(K[index]) << endl;
            //     cout << index << endl;
            //     cout << "~~~~\n";
            // }
            if (index!=-1)
            {
                // if (i<100000)
                // {
                //     cout << std::bitset<64>(tar) << endl;
                //     cout << std::bitset<64>(K[index]) << endl;
                //     cout << "~~~~\n";
                // }


                tmp_i = io_info[index];

                //if is multip in
                if (tmp_i&mask_in)
                {
                    // cout << std::bitset<64>(tar) << endl;
                    // cout << std::bitset<64>(io_info[index]) << endl;
                    BCN_index = (tmp_i&mask_index);
                    BCN_len = ((tmp_i&mask_in) >> 32);
                    tk = BCN_index;
                    //the BC index is start from 1
                    while (tk<BCN_size && BCN[tk])
                        ++tk;
                    if (tk < (BCN_index + BCN_len))
                    {
                        // for (size_t u=BCN_index; u<BCN_index+BCN_len; ++u)
                        //     cout << BCN[u] << endl;
                        if (i>=kmer)
                            BCN[tk] = (bc_index << 3)|get_c[dna_f[i-kmer]];
                        else 
                            BCN[tk] = (bc_index << 3)|4;
                        // if (tk + 1 == BCN_index + BCN_len)
                        // {
                        //     cout << std::bitset<64>(tar) << endl;
                        //     cout << index << endl;
                        //     for (size_t u=BCN_index; u<BCN_index+BCN_len; ++u)
                        //     cout << std::bitset<64>(BCN[u]) << " " << (BCN[u]>>3) << endl;
                        //     cout << "---\n";
                        //     sum += BCN_len;
                        // }
                        // cout << "---" << tar << endl
                        //      << BCN_index << " " << BCN_len << endl
                        //      << std::bitset<64>(io_info[index]) << endl;
                    }
                    else
                    {
                        cout << "something wrong! attention! build BCN out of range\n";
                        // cout << tar << endl
                        //      << std::bitset<64>(io_info[index]) << endl;
                        // for (size_t u=BCN_index; u<BCN_index+BCN_len; ++u)
                        //     cout << BCN[u] << endl;
                    }

                }
                //if is multip out write dna_f[i+1] in bc
                if (tmp_i&mask_out)
                {
                    // cout << std::bitset<64>(tar) << endl;
                    // cout << std::bitset<64>(io_info[index]) << endl;
                    bc[bc_index++] = dna_f[i+1];
                }
            }
        }
    }
    bc.erase(bc.begin()+bc_index, bc.end());
    // cout << bc << endl;
    // cout << "the init dna string len is: " << dna_z << endl
    // << "the branch string len is: " << bc_index << endl
    // << "bc/dna_z: " << (double)bc_index/dna_z << endl;

    // if (sum != BCN_size) cout << "wrong\n" << endl;
    // for (size_t i=0; i<BCN_size; ++i)
    //     cout << i << " " << BCN[i] << endl;
  	//-----------------------------------use K2 and K2c to construct FM-index
    //-----------------------------------handel the last k-1 k_mer (ATG$, TG$, G$)
    //-----------------------------------insert in BWT, concern as case 1 
    uint64_t mask_code = -1;
    mask_code = mask_code << 61 >> 61;
    last_kmer = 0;

    uint8_t *BWT, tmp_code;
    BWT = new uint8_t[dna_z];
    size_t bwt_index = 0;

    //construct the last kmer
    uint64_t *last_string;
    last_string = new uint64_t[kmer+1];
    for (size_t i=kmer+1; i>=1; --i)
    {
        for (size_t j=dna_z-i; j<dna_z; ++j) tmp = (tmp << 2) | get_c[dna_f[j]];
        last_string[kmer+1-i] = tmp<<(64-(i)*2);
        // cout << last_string[kmer+1-i] << endl;
    }
    sort(last_string, last_string+kmer+1, cmp);
    // for (size_t i=0; i<kmer+1; ++i)
    // {
    //     cout << std::bitset<64>(last_string[i]) << endl;
    // }

    char cc[6]={'A', 'C', 'G', 'T', '$', 'X'};
    for (size_t i=0, j, tmp_index=0, l_index=0, index; i<k2len;)
    {
        //find the io_info and check if mulip in
        // if (i>100000) break;
        is_in=false; is_out=false; tmp_mask=0;
        j = i+1;
        now_kmer = K2[i]&mask_c;
        while (l_index < kmer+1 && (last_string[l_index]&mask_c) <= now_kmer)
        {
            // cout << std::bitset<64>(last_string[l_index]) << endl;
            // cout << cc[((last_string[l_index++]&mask_l)>>62)];

            BWT[bwt_index++] = ((last_string[l_index++]&mask_l)>>62);
            // BWT[bwt_index-1] = 4;
        }
        if (now_kmer >= begin_kmer && begin_kmer >= last_kmer) 
        {
            if (now_kmer == begin_kmer)
            {
                is_in = true;
            }
            else
            {
                // cout << cc[4];
                // cout << std::bitset<64>(begin_kmer) << endl;

                BWT[bwt_index++] = 4;
            }
            //else if the begin_kmer no appear in K2, no need to create the io_info of begin_kmer
        }
        last_kmer = now_kmer;
        while (j<k2len && ((K2[j]&mask_c) == now_kmer))
        {
            if ((K2[j]&mask_l) != (K2[i]&mask_l)) is_in = true;
            if ((K2[j]&mask_r) !=( K2[i]&mask_r)) is_out = true;
            // if (is_in || is_out) printf("asd\n");
            ++j;
        }
        if (is_in || is_out)
        {
            // if (K[tmp_index] != (now_kmer >> (64-kmer2_l+2)))
            // {
            //     cout << std::bitset<64>(K[tmp_index]) << endl;
            //     cout << std::bitset<64>((now_kmer >> (64-kmer2_l+2))) << endl;
            //     PX("asd");
            // }
            tmp_mask = io_info[tmp_index++];
        }
        if (is_in)
        {
            BCN_index = (tmp_mask&mask_index);
            BCN_len = ((tmp_mask&mask_in) >> 32);
            // cout << std::bitset<64>(tmp_mask) << endl;
            // cout << BCN_index <<endl;
            // for (size_t o=BCN_index; o<BCN_index+BCN_len; ++o)
            // {
            //     cout << bc[o] << endl;
            // }
            sort(BCN+BCN_index, BCN+BCN_index+BCN_len, bwt_cmp);
                // PX("asd");
            // cout << "----\n";
            // for (size_t k = i; k<j; ++k)
            // {
            //     cout << std::bitset<64>(K2[k]) << endl;
            // }
            // cout << "--\n";
            for (size_t k=BCN_index; k<BCN_index+BCN_len; ++k)
            {
                BWT[bwt_index++] = BCN[k]&mask_code;
                // cout << std::bitset<64>(K2[k]) << endl;
                // cout << cc[BWT[bwt_index-1]] << endl;
            }
            // cout << "    --------------\n";
        }
        else
        {
            tmp_num=0;  
            // cout << "*****\n"; 
            for (size_t k = i; k<j; ++k)
            {
                // cout << std::bitset<64>(K2[k]) << endl;
                tmp_num += (K2[k]&mask_n);
            }
            tmp_code = ((K2[i]&mask_l)>>62);
            while (tmp_num--)
            {
                BWT[bwt_index++] = tmp_code;
                // cout << cc[BWT[bwt_index-1]];
            }
            // cout << "     *******\n";
        }
        
        
        i = j;
    }
    for (size_t i=0; i<=dna_z; ++i)
    {
        printf("%c", cc[BWT[i]]);
        if ((i+1)%80==0) cout << "\n";
    }
	delete [] K2;
    delete [] K;
    delete [] io_info;
    delete [] BCN;
    delete [] BWT;
    inRef.close();
    out.close();
    fclose(inJF);

    // printf("Time used = %.2f\n",  (double)clock() / CLOCKS_PER_SEC);
    return 0;
}

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
//ie. kmer len 20, kmer2 22, kmer_l, 20*2 kmer2_l 22*2, k2len the K2 len, klen, the K len
size_t kmer, kmer2, kmer_l, kmer2_l, k2len, klen;
uint64_t *K;

bool cmp(const uint64_t a, const uint64_t b)
{
    // return a<b;
    return ((a<<2) < (b<<2));
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
    string dna_name, dna_s, dna_f;
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
        // if (i<100){
        //     printf("%d %s ", tn, tkmer);
        //     cout << std::bitset<64>(tar) << endl;
            
        // }
    }

    //-----------------------------------sort kmer+2
    sort(K2, K2+k2len, cmp);
    


    //-----------------------------------use center kmer to mark multip in & out
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

    bool is_in=false, is_out=false, open=false;
    // cout << std::bitset<64>(mask_out) << endl << std::bitset<64>(mask_in) << endl;
    // std::bitset<64>(mask_index)<< endl;

    size_t theindex = 0;
    klen = 0;
    uint64_t tmp_mask, tmp_num;

    for (size_t i=0; i<kmer; ++i)
    {
        tmp = (tmp << 2) | get_c[dna_f[i]];
    }
    uint64_t first_kmer = tmp << (64 - kmer_l) >> 2;
    size_t first_p=-1;

    for (size_t i=0, j, k; i<k2len;)
    {
        is_in=false; is_out=false; tmp_mask=0; open=false;
        j = i+1;
        if (K2[i]&mask_c == first_kmer) 
        {
            open=true;
            is_in = true;
            if ((K2[i]&mask_r)>>(64-kmer2_l) != get_c[dna_f[kmer]]) is_out = true;
        }
        while (j<k2len && ((K2[j]&mask_c) == (K2[i]&mask_c)))
        {
            if ((K2[j]&mask_l) != (K2[i]&mask_l)) is_out = true;
            if ((K2[j]&mask_r) != (K2[i]&mask_r)) is_in = true;
            // if (is_in || is_out) printf("asd\n");
            ++j;
        }
        else if (is_in || is_out)
        {
            // cout << "Bin\n";
            if (is_out) tmp_mask = tmp_mask|mask_out;
            if (is_in)
            {
                tmp_num=0;
                for (k = i; k<j; ++k)
                {
                    tmp_num += (K2[k]&mask_n);
                }
                tmp_mask = tmp_mask|theindex;
                if (open) 
                {
                    ++tmp_num;
                    first_p = theindex;
                }
                tmp_mask = tmp_mask|(tmp_num<<32);
                theindex += tmp_num;
            }
            io_info[klen] = tmp_mask;
            K[klen] = (K2[i]&mask_c) >> (64-kmer2_l+2);
            ++klen;
            // if (klen<100){
            //     // cout << "Asd";

            //     cout << i << " " << j << endl << std::bitset<64>(io_info[klen-1]) << endl << std::bitset<64>(io_info[klen-1]<<1>>1>>32) << endl
            //     << (io_info[klen-1]<<1>>1>>32) << endl;
                
            // }
        }
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
    string bc(dna_f);
    uint64_t bc_index=1;
    uint64_t BCN_index, BCN_len, tk, tmp_i;
    bool flag=true;
    for (size_t i=0, index; i<dna_z; ++i)
    {
        tmp = (tmp << 2) | get_c[dna_f[i]];
        if (i>=kmer-1)
        {
            tar = tmp&mask_k;
            index = search_k(tar);
            if (flag)
            {
                if (index!=-1)
                {
                    tmp_i = io_info[index];
                    if (tmp_i&mask_out)
                    {
                        if (i+1<dna_z)
                        {
                            bc[bc_index++] = dna_f[i+1];
                        }
                    }
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
                            tk++;
                        if (tk < (BCN_index + BCN_len))
                        {
                            the_first_p = tk;
                            BCN[tk] = (bc_index << 2);
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
                }
                flag=false;
            }
            if (index!=-1)
            {
                // cout << std::bitset<64>(tar) << endl;
                // cout << index << endl;
                // if (index != -1)
                // cout << std::bitset<64>(K[index]) << endl;
                // cout << "~~~~\n";



                //if is multip out write dna_f[i+1] in bc
                tmp_i = io_info[index];
                if (tmp_i&mask_out)
                {
                    // cout << std::bitset<64>(tar) << endl;
                    // cout << std::bitset<64>(io_info[index]) << endl;
                    if (i+1<dna_z)
                    {
                        bc[bc_index++] = dna_f[i+1];
                    }

                }
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
                        tk++;
                    if (tk < (BCN_index + BCN_len))
                    {
                        // for (size_t u=BCN_index; u<BCN_index+BCN_len; ++u)
                        //     cout << BCN[u] << endl;
                        if (i>=kmer)
                            BCN[tk] = (bc_index << 2)|get_c[dna_f[i-kmer]];
                        else 

                        // if (tk + 1 == BCN_index + BCN_len)
                        // for (size_t u=BCN_index; u<BCN_index+BCN_len; ++u)
                        //     cout << BCN[u] << endl;
                        // cout << "---\n";
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
            }
        }
    }

    // for (size_t i=0; i<BCN_size; ++i)
    //     cout << i << " " << BCN[i] << endl;
  	//-----------------------------------use K2 and K2c to construct FM-index
    //-----------------------------------handel the last k-1 k_mer (ATG$, TG$, G$)
    //-----------------------------------insert in BWT, concern as case 1 
    uint64_t mask_code = -1;
    mask_code = mask_code << 62 >> 62;
    for (size_t i=0, j, tmp_index=0, index; i<k2len;)
    {
        //find the io_info and check if mulip in
        is_in=false; is_out=false; tmp_mask=0;
        is_bcn=false;
        j = i+1;
        while (j<k2len && ((K2[j]&mask_c) == (K2[i]&mask_c)))
        {
            if ((K2[j]&mask_l) != (K2[i]&mask_l)) is_out = true;
            if ((K2[j]&mask_r) !=( K2[i]&mask_r)) is_in = true;
            // if (is_in || is_out) printf("asd\n");
            ++j;
        }
        if (is_in || is_out)
        {
            tmp_mask = io_info[tmp_index++];
            if (is_in)
            {
                theindex = (tmp_mask&mask_index);
                thelen = (tmp_mask&mask_in) >> 32;
                sort(BCN+i, BCN+j, bcn_cmp);
                is_bcn=true;
                for (size_t k=i; k<j; ++k)
                {
                    //append BWT
                    if (check_$)
                    {

                    }
                    BWT[bwt_index++] = (BCN[k]&mask_code);
                }
            }
        }
        //append BWT string
        if (!is_bcn)
        {
            for (size_t k=i, t, l; k<j; ++k)
            {
                t = (K2[k]&mask_n);
                l = ((K2[k]&mask_l) >> 62);
                while (t--)
                {
                    BWT[bwt_index++] = l;
                    
                }
            }
        }
        i = j;
    }

	delete [] K2;
    delete [] K;
    delete []io_info;
    inRef.close();
    out.close();
    fclose(inJF);

    printf("Time used = %.2f\n",  (double)clock() / CLOCKS_PER_SEC);
    return 0;
}

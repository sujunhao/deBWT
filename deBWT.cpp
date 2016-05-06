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

uint32_t get_c[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                    0, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                    3};

size_t KMER=10;
uint32_t K_mask = 0xffffffff;
inline bool cmp(uint32_t a, uint32_t b)
{
	return (a & K_mask) < (b & K_mask);
}     

int main()
{
	char i_f[100]="E.coli.fa";
	char o_f[100]="out";
    ifstream inRef;
    ofstream out;
    inRef.open(i_f);
    out.open(o_f);



    //-----------------------------------get dna_ref string
    string dna_name, dna_s, dna_f;
    getline(inRef, dna_name);
    while (inRef >> dna_s) dna_f.append(dna_s);
    size_t dna_z=dna_f.size();


    //-----------------------------------k_mer counting
    size_t kmer=KMER, kmer2=((KMER+2)<<1);
    uint32_t *K2, *K2c;

    K2c = new uint32_t[1<<kmer2];
    memset(K2c, 0, sizeof(K2c));

    K2 = new uint32_t[1<<kmer2];
    size_t K2_index=0;

	uint32_t tmp = 0, tar;
    uint32_t MASK = 0xffffffff;
    //init mask
    MASK = MASK >> (32 - kmer2);

    // PX(std::bitset<32>(1<<kmer2));
    for (size_t i=0; i<dna_z; ++i)
    {
    	tmp = (tmp << 2) | get_c[dna_f[i]];
        tar = tmp & MASK;
        if (i>=kmer2-1)
        {
        	// PX(std::bitset<32>(tar));
        	if (K2c[tar] == 0) K2[K2_index++] = tar;
        	++K2c[tar];
        }
    }
    // PX(1<<kmer2);
    cout << K2_index << endl;

    //-----------------------------------sort k_mer 
    // PX(std::bitset<32>(MASK));
    K_mask = MASK >> 2;
    sort(K2, K2+K2_index, cmp);
    // for (size_t i=0; i<100; ++i) cout << std::bitset<32>(K2[i]) << endl;


    //-----------------------------------for each k_mer in dna string, use binary search to find k_mer in K
    //-----------------------------------note that K len is k_mer+2
    //-----------------------------------if is multip out, construct branch code (index++)
    //-----------------------------------if is multip in, store the branch code index(**** index should +1 ****)


  	//-----------------------------------use K2 and K2c to construct FM-index
    //-----------------------------------handel the last k-1 k_mer (ATG$, TG$, G$)
    //-----------------------------------insert in BWT, concern as case 1 

	delete [] K2;
    delete [] K2c;
    inRef.close();
    out.close();

    printf("Time used = %.2f\n",  (double)clock() / CLOCKS_PER_SEC);
    return 0;
}

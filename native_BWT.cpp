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
size_t dna_z;

bool cmp(const size_t ta, const size_t tb)
{
    size_t a = ta, b = tb;
    bool flag;
    while (a<dna_z && b<dna_z && dna_f[a]==dna_f[b])
    {
        ++a;
        ++b;
    }
    if (a==dna_z) flag = true;
    if (b==dna_z) flag = false;
    flag = dna_f[a]<dna_f[b];
    return flag;
}

int main()
{
    char i_f[100]="E.coli.fa";
    ifstream inRef;;

    inRef.open(i_f);



    //-----------------------------------get dna_ref string
    getline(inRef, dna_name);
    while (inRef >> dna_s) dna_f.append(dna_s);
    dna_f = "CGCCTTAGTAAGTGATTTTC";
    dna_f.append("$");
    dna_z=dna_f.size();

    uint64_t *thei;
    thei = new uint64_t[dna_z];
    for (size_t i=0; i<dna_z; ++i)
        thei[i] = i;

    sort(thei, thei+dna_z, cmp);
    for (size_t i=0; i<dna_z; ++i)
    {
        cout << dna_f[(thei[i]+dna_z-1)%dna_z];
        if ((i+1)%80==0) cout << "\n";
    }
    inRef.close();

    return 0;
}

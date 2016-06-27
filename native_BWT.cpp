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
string BWT;


//sort the SA[]
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

//get the number of char c that appear before and in rank[pos] 
size_t get_rank(size_t pos, int c, uint64_t **rank, uint64_t gap)
{
    if ((pos+1) % gap == 0)
        return rank[c][pos/gap];

    size_t offset;
    if (pos/gap > 0) offset = rank[c][pos/gap - 1];
    else offset = 0;

    for (size_t i=(pos/gap)*gap; i<=pos; ++i)
    {
        if (BWT[i] != '$' && get_c[BWT[i]]==c)
            ++offset;
    }
    return offset;

}

int main()
{
    char i_f[100]="E.coli.fa";
    ifstream inRef;;

    inRef.open(i_f);



    //-----------------------------------get dna_ref string
    getline(inRef, dna_name);
    while (inRef >> dna_s) dna_f.append(dna_s);
    dna_f = "AGCTAGGGTC";
    dna_f.append("$");
    dna_z=dna_f.size();

    uint64_t *thei, *SA;
    thei = new uint64_t[dna_z];
    // SA = new uint64_t[dna_z];
    for (size_t i=0; i<dna_z; ++i)
        thei[i] = i;

    //sort the SA
    sort(thei, thei+dna_z, cmp);

    BWT = dna_f;

    //build the BWT
    for (size_t i=0; i<dna_z; ++i)
    {
        BWT[i] = dna_f[(thei[i]+dna_z-1)%dna_z];
        // SA[thei[i]] = i;
        // cout << i << " " << thei[i] << endl;
        // cout << BWT[i];
        // if ((i+1)%80==0) cout << "\n";
    }
    cout << "ref: " << dna_f << endl << "bwt: " << BWT << endl;
    puts("");

    //count the ATCG in C[] and rank[pos, char] to indicate 0-pos have how many same char in BWT
    uint64_t C[4] = {0, 0, 0, 0};
    uint64_t *(rank[4]), gap = 4;
    for (size_t i=0; i<4; ++i) rank[i] = new uint64_t[dna_z/gap];
        printf("SA_i  BWT_i BWT   Origin\n");
    for (size_t i=0; i<dna_z; ++i)
    {
        if (BWT[i]!='$' ) 
            ++C[get_c[BWT[i]]];
        if ((i+1)%gap == 0)
        {
            rank[0][i/gap] = C[0];
            rank[1][i/gap] = C[1];
            rank[2][i/gap] = C[2];
            rank[3][i/gap] = C[3];

        }

        printf("%3d   %3d   ", (int)thei[i], (int)i);
        cout << " " << BWT[i] << "    " ;
        for (size_t u=0, t=thei[i]; u<dna_z; ++u) cout << dna_f[(u+t)%dna_z];
        cout << endl;
        // if ((i+1)%gap == 0) cout << " " << rank[0][i/gap] << " " << rank[1][i/gap] << " " << rank[2][i/gap] << " " << rank[3][i/gap] << endl;
        // else cout << endl;
    }
    puts("");

    uint64_t Count[4];
    for (size_t i=0; i<4; i++)
    {
        if (i==0) Count[i]=1;
        else Count[i] = Count[i-1] + C[i-1];
        // cout << Count[i] << endl;
    }



    string P="CTAGGG";
    cout << "find the pattern: " << P << endl;
    size_t P_size= P.size(), ec=get_c[P[P_size-1]];
    size_t first=Count[ec], last=Count[ec]+C[ec]-1;

    size_t P_i=P_size-2, offset;

    cout << "\nfirst and last index: " << first << " " << last << endl;
    // for (size_t i=0; i<dna_z; ++i)
    // {
    //     cout << " " << get_rank(i, 0, rank, gap) << " " << get_rank(i, 1, rank, gap) << " " << get_rank(i, 2, rank, gap) << " " << get_rank(i, 3, rank, gap) << endl;
    // }
    while (P_size>1 && P_i >= 0 && first <= last)
    {
        offset = Count[get_c[P[P_i]]];
        first = offset + get_rank(first-1, get_c[P[P_i]], rank, gap);
        last = offset + get_rank(last, get_c[P[P_i]], rank, gap)-1;
        cout << "first and last index: "<< first << " " << last << endl;

        // cout << get_rank(first-1, get_c[P[P_i]], rank, gap) << " " << first << " " << last << endl;
        if (P_i==0) break;
        --P_i;
    }

    
    if (first <= last)
    {
        for (size_t i=first, sa; i<=last; ++i)
        {
            sa = thei[i];
            cout << sa << " " << dna_f.substr(sa, P_size) << endl;
        }
    }
    else 
    {
        cout << "target no find.\n";
    }

    inRef.close();

    return 0;
}


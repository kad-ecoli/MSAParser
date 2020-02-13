const char* docstring=""
"cleanFastaBody metaclust.fasta clean.fasta\n"
"    Read 'metaclust.fasta', remove non-standard amino acids, and save\n"
"    the clean sequence to 'clean.fasta'\n"
;

#include <iostream>
#include <vector>
#include <string>
#include <cstring>
#include <fstream>
#include <algorithm>

using namespace std;

const string aa_list="ACDEFGHIKLMNPQRSTVWY";

int cleanFastaBody(const string & infile, const string & outfasta)
{
    ifstream fp_in;
    ofstream fp_outfasta;
    if (infile!="-") fp_in.open(infile.c_str(),ios::in);
    if (outfasta!="-") fp_outfasta.open(outfasta.c_str(),ofstream::out);
    string line,sequence;
    size_t seq_num=0;
    int j;
    while ((infile!="-")?fp_in.good():cin.good())
    {
        if (infile!="-") getline(fp_in,line);
        else getline(cin,line);

        if (line.size()>0 && line[0]=='>')
        {
            if (outfasta!="-") fp_outfasta<<line<<endl;
            else cout<<line<<endl;
            seq_num++;
        }
        else if (line.size())
        {
            sequence.clear();
            for (j=0;j<line.size();j++)
                if (binary_search(aa_list.begin(),aa_list.end(),line[j]))
                    sequence+=line[j];

            if (outfasta!="-") fp_outfasta<<sequence<<endl;
            else cout<<sequence<<endl;
        }
    }

    line.clear();
    sequence.clear();

    fp_in.close();
    fp_outfasta.close();
    return seq_num;
}

int main(int argc, char **argv)
{
    /* parse commad line argument */
    if(argc!=3)
    {
        cerr<<docstring;
        return 0;
    }
    cleanFastaBody(argv[1],argv[2]);
    return 0;
}

const char* docstring=""
"cleanFastaHeader metaclust.fasta clean.fasta header.tsv\n"
"    Read 'metaclust.fasta', rename sequence name as 1, 2, 3 ..., save\n"
"    the renamed fasta to 'clean.fasta', and the sequence name mapping\n"
"    to 'header.tsv'.\n"
;

#include <iostream>
#include <vector>
#include <string>
#include <cstring>
#include <fstream>

using namespace std;

int cleanFastaHeader(const string & infile, 
    const string & outfasta, const string & outheader)
{
    ifstream fp_in;
    ofstream fp_outfasta;
    ofstream fp_outheader;
    if (infile!="-") fp_in.open(infile.c_str(),ios::in);
    if (outfasta!="-") fp_outfasta.open(outfasta.c_str(),ofstream::out);
    if (outheader!="-") fp_outheader.open(outheader.c_str(),ofstream::out);
    string line,name;
    int seq_idx=0;
    while ((infile!="-")?fp_in.good():cin.good())
    {
        if (infile!="-") getline(fp_in,line);
        else getline(cin,line);

        if (line.length()>0 && line[0]=='>')
        {
            seq_idx++;
            if (outfasta!="-") fp_outfasta<<'>'<<seq_idx<<endl;
            else cout<<'>'<<seq_idx<<endl;

            if (outheader!="-") 
                fp_outheader<<seq_idx<<'\t'<<line.substr(1)<<endl;
            else cout<<seq_idx<<'\t'<<line.substr(1)<<endl;
        }
        else
        {
            if (outfasta!="-") fp_outfasta<<line<<endl;
            else cout<<line<<endl;
        }
    }

    fp_in.close();
    fp_outfasta.close();
    fp_outheader.close();
    return seq_idx;
}

int main(int argc, char **argv)
{
    /* parse commad line argument */
    if(argc!=4)
    {
        cerr<<docstring;
        return 0;
    }
    cleanFastaHeader(argv[1],argv[2],argv[3]);
    return 0;
}

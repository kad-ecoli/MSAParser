const char* docstring=""
"fastaCov 0.75 uniclust.fasta > 0.75.fasta\n"
"    remove sequences in uniclust.fasta with more than 1-0.75=25\% gaps\n"
;

#include <iostream>
#include <vector>
#include <string>
#include <cstring>
#include <cstdlib>
#include <fstream>

using namespace std;

/* calculate portion of non-gap residues */
float calCov(const string sequence)
{
    if (sequence.length()==0) return 0.;
    float cov=0;
    for (int i=0;i<sequence.length();i++)
        cov+=(sequence[i]!='-' && sequence[i]!='.');
    return cov/sequence.length();
}

int fastaCov(const float cov_cut, string infile="-", string outfile="-")
{
    ifstream fp_in;
    ofstream fp_out;
    if (infile!="-") fp_in.open(infile.c_str(),ios::in);
    if (outfile!="-") fp_out.open(outfile.c_str(),ofstream::out);
    string sequence,header,line;
    int nseqs=0;
    while ((infile!="-")?fp_in.good():cin.good())
    {
        if (infile!="-") getline(fp_in,line);
        else getline(cin,line);

        if (line.length()==0) continue;
        if (line[0]=='>')
        {
            if (calCov(sequence)>=cov_cut)
            {
                if (outfile!="-")
                    fp_out<<header<<endl<<sequence<<endl;
                else
                    cout<<header<<endl<<sequence<<endl;
            }
            header=line;
            sequence.clear();
            nseqs++;
        }
        else
            sequence+=line;
    }
    fp_in.close();
    if (calCov(sequence)>=cov_cut)
    {
        if (outfile!="-") fp_out<<header<<endl<<sequence<<endl;
        else                cout<<header<<endl<<sequence<<endl;
    }
    fp_out.close();
    header.clear();
    sequence.clear();
    return nseqs;
}

int main(int argc, char **argv)
{
    /* parse commad line argument */
    if(argc<2)
    {
        cerr<<docstring;
        return 0;
    }
    float cov_cut=atof(argv[1]);
    string infile=(argc<=2)?"-":argv[2];
    string outfile=(argc<=3)?"-":argv[3];
    fastaCov(cov_cut,infile,outfile);
    return 0;
}

const char* docstring=""
"unaligna3m input.a3m output.fasta\n"
"    convert a3m format alignment input.a3m to fasta format raw sequence\n"
"    file output.fasta, without any gap.\n"
;

#include <iostream>
#include <vector>
#include <string>
#include <cstring>
#include <cstdlib>
#include <fstream>

using namespace std;

int unaligna3m(const string infile="-", const string outfile="-")
{
    ifstream fp_in;
    ofstream fp_out;
    if (infile!="-") fp_in.open(infile.c_str(),ios::in);
    if (outfile!="-") fp_out.open(outfile.c_str(),ofstream::out);
    string sequence,line;
    int nseqs=0;
    int i=0;
    int upper='A'-'a';
    while ((infile!="-")?fp_in.good():cin.good())
    {
        if (infile!="-") getline(fp_in,line);
        else getline(cin,line);

        if (line.length()==0) continue;
        if (line[0]=='>')
        {
            if (sequence.length()>0)
            {
                if (outfile!="-") fp_out<<sequence<<endl;
                else                cout<<sequence<<endl;
            }
            sequence.clear();
            nseqs++;
            if (outfile!="-") fp_out<<line<<endl;
            else                cout<<line<<endl;
        }
        else
        {
            for (i=0;i<line.size();i++)
            {
                if ('A'<=line[i] && line[i]<='Z') sequence+=line[i];
                else if ('a'<=line[i] && line[i]<='z') sequence+=line[i]+upper;
            }
        }
    }
    fp_in.close();
    if (outfile!="-") fp_out<<sequence<<endl;
    else                cout<<sequence<<endl;
    fp_out.close();
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
    string infile=argv[1];
    string outfile=(argc<=2)?"-":argv[2];
    unaligna3m(infile,outfile);
    return 0;
}

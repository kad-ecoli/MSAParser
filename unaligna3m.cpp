const char* docstring=""
"unaligna3m input.a3m output.fasta\n"
"    convert a3m format alignment input.a3m to fasta format raw sequence\n"
"    file output.fasta, without any gap.\n"
"\n"
"unaligna3m input.a3m output.fasta [option]\n"
"    additional clean up options:\n"
"    1   - only keep the first word in sequence header\n"
"    10  - remove adjacent duplicated sequences (assuming fasta is sorted)\n"
"    100 - remove sequence marked with _consensus\n"
;

#include <iostream>
#include <vector>
#include <string>
#include <cstring>
#include <cstdlib>
#include <fstream>

using namespace std;

int unaligna3m(const string infile="-", const string outfile="-",
    const int option=0)
{
    ifstream fp_in;
    ofstream fp_out;
    if (infile!="-") fp_in.open(infile.c_str(),ios::in);
    if (outfile!="-") fp_out.open(outfile.c_str(),ofstream::out);
    string sequence,line;
    int nseqs=0;
    int i=0;
    int upper='A'-'a';
    bool noconsensus =(option      )>=100;
    bool deduplicate =(option % 100)>=10;
    bool short_header=(option % 10 )>=1;
    bool isconsensus =false;
    bool wasconsensus=false;
    string header="";
    string prev_sequence="";
    while ((infile!="-")?fp_in.good():cin.good())
    {
        if (infile!="-") getline(fp_in,line);
        else getline(cin,line);

        if (line.length()==0 || line[0]=='#') continue;
        if (line[0]=='>')
        {
            if (sequence.length()>0)
            {
                if (wasconsensus) wasconsensus=false;
                if (noconsensus && isconsensus) 
                {
                    isconsensus=false;
                    wasconsensus=true;
                }
                else if (deduplicate && prev_sequence==sequence);
                else
                {
                    if (outfile!="-") fp_out<<header<<'\n'<<sequence<<endl;
                    else                cout<<header<<'\n'<<sequence<<endl;
                }
            }

            if (short_header)
            {
                for (i=0;i<line.size();i++)
                {
                    if (line[i]==' ' or line[i]=='\t' or line[i]==0)
                    {
                        line=line.substr(0,i);
                        break;
                    }
                }
            }

            if (noconsensus && line.size()>=10 &&
                line.substr(line.size()-10)=="_consensus")
            {
                isconsensus=true;
                continue;
            }
            header=line;

            if (deduplicate && (noconsensus==false || wasconsensus==false))
                prev_sequence=sequence;
            sequence.clear();
            nseqs++;
        }
        else
        {
            for (i=0;i<line.size();i++)
            {
                if ('A'<=line[i] && line[i]<='Z') sequence+=line[i];
                else if ('a'<=line[i] && line[i]<='z') sequence+=line[i]+upper;
                else if (line[i]==0) break;
            }
        }
    }
    fp_in.close();

    if (sequence.length()>0)
    {
        if (wasconsensus) wasconsensus=false;
        if (noconsensus && isconsensus) 
        {
            isconsensus=false;
            wasconsensus=true;
        }
        else if (deduplicate && prev_sequence==sequence);
        else
        {
            if (outfile!="-") fp_out<<header<<'\n'<<sequence<<endl;
            else                cout<<header<<'\n'<<sequence<<endl;
        }
    }

    fp_out.close();
    header.clear();
    sequence.clear();
    prev_sequence.clear();
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
    int    option =(argc<=3)?0:atoi(argv[3]);
    unaligna3m(infile,outfile,option);
    return 0;
}

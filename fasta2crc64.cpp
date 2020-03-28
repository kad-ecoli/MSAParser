const char* docstring=""
"fasta2crc64 metaclust.fasta metaclust.crc64\n"
"    convert FASTA format alignment metaclust.fasta\n"
"    to CRC64 value table metaclust.crc64\n"
;

#include <iostream>
#include <vector>
#include <string>
#include <cstring>
#include <cstdlib>
#include <fstream>

using namespace std;

/* The following code block is modified based on David Jones
 * Improved calculation of CRC64 values for protein sequences
 * Modified from http://bioinfadmin.cs.ucl.ac.uk/downloads/crc64/ */

#define POLY64REV     0x95AC9329AC4BC9B5ULL
#define INITIALCRC    0xFFFFFFFFFFFFFFFFULL

void crc64(const char *seq, char *res)
{
    int i, j, low, high;
    unsigned long long crc = INITIALCRC, part;
    static int init = 0;
    static unsigned long long CRCTable[256];
    
    if (!init)
    {
        init = 1;
        for (i = 0; i < 256; i++)
        {
            part = i;
            for (j = 0; j < 8; j++)
            {
                if (part & 1) part = (part >> 1) ^ POLY64REV;
                else part >>= 1;
            }
            CRCTable[i] = part;
        }
    }
    
    while (*seq)
        crc = CRCTable[(crc ^ *seq++) & 0xff] ^ (crc >> 8);

    /* The output is done in two parts to avoid problems with 
     * architecture-dependent word order */
    low = crc & 0xffffffff;
    high = (crc >> 32) & 0xffffffff;
    sprintf (res, "%08X%08X", high, low);
    return;
}

/* End of code block for CRC64 */

int fasta2pfam(const string infile="-", const string outfile="-")
{
    ifstream fp_in;
    ofstream fp_out;
    if (infile!="-") fp_in.open(infile.c_str(),ios::in);
    if (outfile!="-") fp_out.open(outfile.c_str(),ofstream::out);
    string sequence,line,header;
    int nseqs=0;
    char result[20];
    while ((infile!="-")?fp_in.good():cin.good())
    {
        if (infile!="-") getline(fp_in,line);
        else getline(cin,line);

        if (line.length()==0) continue;
        if (line[0]=='>')
        {
            if (sequence.length()>0)
            {
                crc64(sequence.c_str(),result);
                if (outfile!="-") fp_out<<header<<'\t'<<result<<endl;
                else                cout<<header<<'\t'<<result<<endl;
            }
            sequence.clear();
            header=line.substr(1,line.size()-1);
            nseqs++;
        }
        else sequence+=line;
    }
    fp_in.close();
    if (sequence.length())
    {
        crc64(sequence.c_str(),result);
        if (outfile!="-") fp_out<<header<<'\t'<<result<<endl;
        else                cout<<header<<'\t'<<result<<endl;
    }
    fp_out.close();
    sequence.clear();
    header.clear();
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
    fasta2pfam(infile,outfile);
    return 0;
}

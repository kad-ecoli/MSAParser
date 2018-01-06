const char* docstring=""
"trimMSA symfrac seq.aln pos_map trim.aln\n"
"    remove positions with > symfrac gaps from seq.aln\n"
"    save trimmed aligment to trim.aln. save position mapping between\n"
"    trim.aln and seq.aln to pos_map\n"
"\n"
"trimMSA symfrac seq.aln pos_map trim.aln 1\n"
"    apart from first sequence, set every position with > symfrac gaps\n"
"    to lower case\n"
;

#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <fstream>
#include <cctype> 

using namespace std;

/* parse PSICOV format MSA */
int parse_aln(const char *filename, vector <string> & aln)
{
    ifstream fp;
    fp.open(filename,ios::in);
    string sequence;
    while (fp.good())
    {
        getline(fp,sequence);
        if (sequence.length()==0) continue;
        aln.push_back(sequence);
    }
    return aln.size();
}

/* calculate fraction of gaps for each position */
int calGapFrac(vector <string> & aln, vector <float> & gapfrac_vec)
{
    int nseqs=aln.size();
    int L=aln[0].length();
    gapfrac_vec.assign(L,0.);
    int n,i;
    for (n=0;n<nseqs;n++)
        for (i=0;i<L;i++) 
            gapfrac_vec[i]+=(aln[n][i]=='-' || aln[n][i]=='.');
    for (i=0;i<L;i++) gapfrac_vec[i]/=nseqs;
    return L;
}

/* trim alignment such that no position has > symfrac fraction of gaps */
int trimMSA(vector <string> & aln, vector <float> & gapfrac_vec, 
    vector <int> & pos_map,const float symfrac,const int output_mode=0)
{
    int L=calGapFrac(aln,gapfrac_vec);
    int nseqs=aln.size();
    int i,j,n;
    vector <int> nonpos_map;
    for (i=0,j=-1;i<L;i++)
    {
        if (gapfrac_vec[i]<=symfrac) pos_map.push_back(i);
        else nonpos_map.push_back(i);
    }

    if (output_mode==0)
    {
        string sequence;
        for (n=0;n<nseqs;n++)
        {
            for (j=0;j<pos_map.size();j++) 
                sequence+=aln[n][pos_map[j]];
            aln[n]=sequence;
            sequence.clear();
        }
    }
    else if (output_mode==1)
    {
        for (n=0;n<nseqs;n++)
            for (j=0;j<nonpos_map.size();j++)
                aln[n][nonpos_map[j]]=tolower(aln[n][nonpos_map[j]]);
    }
    return nseqs;
}

int main(int argc, char **argv)
{
    /* parse commad line argument */
    if(argc<5)
    {
        cerr<<docstring;
        return 0;
    }
    float symfrac=atof(argv[1]);
    vector <string> aln;
    vector <float> gapfrac_vec;
    vector <int> pos_map;
    int output_mode=0;
    if (argc>5) output_mode=atoi(argv[5]);

    int nseqs=parse_aln(argv[2],aln);
    trimMSA(aln,gapfrac_vec,pos_map,symfrac,output_mode);

    ofstream fp;
    fp.open(argv[3],ios::out);
    for (int j=0;j<pos_map.size();j++)
        fp<<j<<'\t'<<pos_map[j]<<'\t'<<gapfrac_vec[pos_map[j]]<<endl;
    fp.close();

    fp.open(argv[4],ios::out);
    for (int n=0;n<nseqs;n++) fp<<aln[n]<<endl;
    fp.close();
    return 0;
}

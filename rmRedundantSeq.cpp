const char* docstring=""
"rmRedundantSeq id_cut cov_cut metaclust.aln\n"
"    remove sequences with < cov_cut sequence coverage or with > id_cut\n"
"    sequence identity to other sequences in alignment file metaclust.aln\n"
"\n"
"rmRedundantSeq id_cut cov_cut uniclust.aln metaclust.aln\n"
"    remove sequences in metaclust.aln with < cov_cut sequence coverage or\n"
"    or with > id_cut  sequence identity to other sequences in either\n"
"    uniclust.aln or metaclust.aln\n"
;

#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <fstream>

using namespace std;

string aa_list="ACDEFGHIKLMNPQRSTVWY";

/* parse single MSA at seqID & coverage cutoff id_cut & cov_cut */
int parse_single_aln(const char *filename, vector <string> & aln, 
    const float id_cut, const float cov_cut)
{
    int L=0; // alignment length
    int Ldiff=0,maxLdiff=0; // length of different positions
    int i,n; // position, sequence index
    int Lali; // number of non-gap positions
    float seqID1,seqID2; // seqID with previous sequences
    vector <int> is_not_gap_vec; // if a positions is an amino acid 
    vector<vector<int> > is_not_gap_mat; // if a position is an amino acid
    ifstream fp;
    fp.open(filename,ios::in);
    string sequence;
    while (fp.good())
    {
        getline(fp,sequence);
        if (sequence.length()==0) continue;
        if (L==0)
        {
            L=sequence.length();
            is_not_gap_vec.assign(L,0);
        }
        Lali=0;
        for (i=0;i<L;i++)
        {
            if (aa_list.find(sequence[i])!=string::npos)
            {
                Lali+=1;
                is_not_gap_vec[i]=1;
            }
            else
                is_not_gap_vec[i]=0;
        }
        if (Lali<cov_cut*L) continue;
        if (id_cut>1)
        {
            aln.push_back(sequence);
        }
        else
        {
            maxLdiff=(1.-id_cut)*Lali;
            for (n=0;n<aln.size();n++)
            {
                Ldiff=0;
                //seqID1=0;
                for (i=0;i<L;i++)
                {
                    //seqID1+=(is_not_gap_vec[i])*(aln[n][i]==sequence[i]);
                    Ldiff+=(is_not_gap_vec[i])*(aln[n][i]!=sequence[i]);
                    if (Ldiff>maxLdiff) break;
                }
                //if (seqID1>id_cut*Lali) break;
                if (Ldiff<=maxLdiff) break;
            }
            if (n==aln.size()) 
            {
                aln.push_back(sequence);
                is_not_gap_mat.push_back(is_not_gap_vec);
            }
        }
    }
    fp.close();
    is_not_gap_mat.clear();
    return aln.size();
}

/* parse second MSA at seqID & coverage cutoff id_cut & cov_cut */
int parse_two_aln(const char *query_filename, const char *filename, 
    vector <string> & query_aln, vector <string> & aln, 
    const float id_cut, const float cov_cut)
{
    string sequence;
    if (query_filename[0]!='-')
    {
        ifstream fp;
        fp.open(query_filename);
        while(fp.good())
        {
            getline(fp,sequence);
            query_aln.push_back(sequence);
        }
        fp.close();
    }
    else
    {
        while(cin.good())
        {
            getline(cin,sequence);
            query_aln.push_back(sequence);
        }
    }
    int L=query_aln[0].length(); // alignment length
    int query_nseqs=query_aln.size(); // number of sequences

    int i,n; // position, sequence index
    int Lali; // number of non-gap positions
    int Ldiff=0,maxLdiff=0; // length of different positions
    float max_seqID,seqID1,seqID2; // seqID with previous sequences
    vector <int> is_not_gap_vec(L,0); // if a positions is an amino acid 
    vector<vector<int> > is_not_gap_mat; // if a position is an amino acid
    ifstream fp;
    fp.open(filename,ios::in);
    while (fp.good())
    {
        getline(fp,sequence);
        if (sequence.length()==0) continue;
        Lali=0;
        for (i=0;i<L;i++)
        {
            if (aa_list.find(sequence[i])!=string::npos)
            {
                Lali+=1;
                is_not_gap_vec[i]=1;
            }
            else
                is_not_gap_vec[i]=0;
        }
        if (Lali<cov_cut*L) continue;
        if (id_cut>1)
        {
            aln.push_back(sequence);
        }
        else
        {
            //max_seqID=0;
            maxLdiff=(1.-id_cut)*Lali;
            for (n=0;n<query_aln.size();n++)
            {
                Ldiff=0;
                //seqID1=0;
                for (i=0;i<L;i++)
                {
                    //seqID1+=(is_not_gap_vec[i])*(query_aln[n][i]==sequence[i]);
                    Ldiff+=(is_not_gap_vec[i])*(query_aln[n][i]!=sequence[i]);
                    if (Ldiff>maxLdiff) break;
                }
                if (Ldiff<=maxLdiff) break;
                //if (seqID1>id_cut*Lali)
                //{
                    //max_seqID=seqID1;
                    //break;
                //}
            }
            if (Ldiff<=maxLdiff) continue;
            //if (max_seqID>id_cut*Lali) continue;
            for (n=0;n<aln.size();n++)
            {
                Ldiff=0;
                //seqID1=0;
                for (i=0;i<L;i++)
                {
                    Ldiff+=(is_not_gap_vec[i])*(aln[n][i]!=sequence[i]);
                    //seqID1+=(is_not_gap_vec[i])*(aln[n][i]==sequence[i]);
                    if (Ldiff>maxLdiff) break;
                }
                if (Ldiff<=maxLdiff) break;
                //if (seqID1>id_cut*Lali) break;
            }
            if (n==aln.size()) 
            {
                aln.push_back(sequence);
                is_not_gap_mat.push_back(is_not_gap_vec);
            }
        }
    }
    fp.close();
    is_not_gap_mat.clear();
    return aln.size();
}

int main(int argc, char **argv)
{
    /* parse commad line argument */
    if(argc<4)
    {
        cerr<<docstring;
        return 0;
    }
    float id_cut=atof(argv[1]);
    float cov_cut=atof(argv[2]);
    if (id_cut>1) id_cut/=100.;
    if (cov_cut>1) cov_cut/=100.;
    int nseqs;
    vector <string> query_aln;
    vector <string> aln;
    if (argc<5) // one MSA
        nseqs=parse_single_aln(argv[3],aln,id_cut,cov_cut);
    else
        nseqs=parse_two_aln(argv[3],argv[4],query_aln,aln,id_cut,cov_cut);
    for (int n=0;n<aln.size();n++) cout<<aln[n]<<endl;
    return 0;
}

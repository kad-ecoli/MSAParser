const char* docstring=""
"AlnAaProb seq.aln 0.8\n"
"    Calculate the composition of 20 amino acids at each position for\n"
"    PSICOV format MSA file 'seq.aln'. GREMLIN sequence weight is applied\n"
"    to weight each sequence, using 0.8 (default) sequence identity cutoff.\n"
"\n"
;

#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <fstream>
#include <map>
#include <climits>
#include <iomanip>

using namespace std;

string aa_list="ACDEFGHIKLMNPQRSTVWY";

float AlnAaProb(const string infile, float id_cut,
    vector<vector<float> > &AaProb_mat, vector<float> &seq_weight_vec,
    vector<float> &AaProb_vec)
{
    if (id_cut>1) id_cut/=100.0; // percentage seqID

    /* read aln file */
    vector<string> aln;
    string sequence;
    int L=0;
    ifstream fp;
    if (infile!="-") fp.open(infile.c_str(),ios::in);
    while ((infile!="-")?fp.good():cin.good())
    {
        if (infile!="-") getline(fp,sequence);
        else getline(cin,sequence);

        if (sequence.length()==0) continue;
        if (L==0)
        {
            L=sequence.length();
        }
        else if (sequence.length()!=L)
        {
            cerr<<"ERROR! length not match for sequence "<<aln.size();
            exit(0);
        }

        aln.push_back(sequence);
        if (aln.size()==INT_MAX)
        {
            cerr<<"WARNING! Cannot read beyond sequence number"<<INT_MAX<<endl;
            break;
        }
    }
    fp.close();

    /* compute seqID */
    int seq_num=aln.size();
    vector<int> homo_num_vec(seq_num,1); // seq is always similar to itself
    int m,n; // sequence index
    int i; // residue position
    int Ldiff=0; // length of different residues, including gaps
    int max_Ldiff=L*(1.-id_cut); // max num of different residues for
                                 // two seq to be considered homolog
    for (m=0;m<seq_num;m++)
    {
        for (n=m+1;n<seq_num;n++)
        {
            Ldiff=0;
            for (i=0;i<L;i++)
            {
                Ldiff+=(aln[m][i]!=aln[n][i]);
                if (Ldiff>max_Ldiff) break;
            }
            homo_num_vec[m]+=(Ldiff<=max_Ldiff);
            homo_num_vec[n]+=(Ldiff<=max_Ldiff);
        }
    }

    /* compute seq weight */
    seq_weight_vec.assign(seq_num,0.);
    for (m=0;m<seq_num;m++) seq_weight_vec[m]=1./homo_num_vec[m];
    homo_num_vec.clear();

    /* mapping AA to int */
    map<char,int> aa2int_dict;
    int a;
    for (a=0;a<aa_list.length();a++) aa2int_dict[aa_list[a]]=a;

    /* compute AA prob */
    AaProb_vec.assign(aa_list.length(),0.);
    AaProb_mat.assign(L,AaProb_vec);
    float neff;
    for (i=0;i<L;i++)
    {
        neff=0; // neff at position i, usually smaller than neff of aln
        for (a=0;a<aa_list.length();a++) AaProb_vec[a]=0;
        for (m=0;m<seq_num;m++)
        {
            if (aa2int_dict.count(aln[m][i]))
            {
                neff+=seq_weight_vec[m];
                AaProb_vec[aa2int_dict[aln[m][i]]]+=seq_weight_vec[m];
            }
        }
        for (a=0;a<AaProb_vec.size();a++) AaProb_mat[i][a]=AaProb_vec[a]/neff;
    }

    /* compute neff */
    neff=0;
    for (m=0;m<seq_num;m++) neff+=seq_weight_vec[m];

    /* compute averaged aa prob */
    for (a=0;a<aa_list.length();a++)
    {
        AaProb_vec[a]=0;
        for (i=0;i<L;i++)
        {
            AaProb_vec[a]+=AaProb_mat[i][a];
        }
        AaProb_vec[a]/=L;
    }
    return neff;
}


int main(int argc, char **argv)
{
    /* parse commad line argument */
    float id_cut=0.8; // defined by gremlin
    if(argc<2)
    {
        cerr<<docstring;
        return 0;
    }
    string infile=argv[1];
    if (argc>2) id_cut=atof(argv[2]);
    vector<vector<float> > AaProb_mat;
    vector<float> seq_weight_vec;
    vector<float> AaProb_vec;
    float neff=AlnAaProb(infile,id_cut,AaProb_mat,seq_weight_vec,AaProb_vec);

    /* print Aa Prob */
    cout<<"#neff="<<setiosflags(ios::fixed)<<setprecision(5)<<neff<<endl;
    cout<<"resi";
    int a,i;
    for (a=0;a<aa_list.length();a++) cout<<'\t'<<aa_list[a];
    cout<<endl;
    for (i=0;i<AaProb_mat.size();i++)
    {
        cout<<i;
        for (a=0;a<aa_list.length();a++) cout<<'\t'
            <<setiosflags(ios::fixed)<<setprecision(5)<<AaProb_mat[i][a];
        cout<<endl;
    }
    cout<<"mean";
    for (a=0;a<aa_list.length();a++) cout<<'\t'
         <<setiosflags(ios::fixed)<<setprecision(5)<<AaProb_vec[a];
    cout<<endl;
    return 0;
}

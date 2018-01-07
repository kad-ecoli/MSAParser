const char* docstring=""
"calNf seq.aln 0.8\n"
"    calculate Nf, number of effective sequence defined in gremlin,\n"
"    using 0.8 (default) sequence identity cutoff.\n"
"\n"
"calNf seq.aln 0.8 0\n"
"    The third optional argument is normalization methods:\n"
"       0 - normalized by L^0.5 (default)\n"
"       1 - normalized by L\n"
"       2 - do not normalize\n"
"      10 - calculate seqID by number of nongap positions\n"
"      20 - do not calculate seqID (overwrite seqID cutoff as 1)\n"
"\n"
"calNf seq.aln 0.8 0 128\n"
"    The fourth optional argument is target Nf. Stop Nf calculation if\n"
"    it is already greater than target Nf. default is 0 (no target Nf)\n"
;

#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <fstream>
#include <cmath>

using namespace std;

string aa_list="ACDEFGHIKLMNPQRSTVWY";

float calNf(const char *filename, float id_cut=0.8, int norm=0,
    float target_Nf=0)
{
    int L=0; // alignment length
    int i,n,m; // position, sequence index, sequence index
    vector <int> is_not_gap_vec; // if a positions is an amino acid 
    vector<vector<int> > is_not_gap_mat; // if a position is an amino acid
    vector<int> Lali_vec;  // length of aligned residues for each sequence
    string sequence;
    vector <string> aln;
    int norm_by_ali=(norm/10);
    norm %=10;

    /* read aln file */
    ifstream fp;
    fp.open(filename,ios::in);
    while (fp.good())
    {
        getline(fp,sequence);
        if (sequence.length()==0) continue;
        if (L==0)
        {
            L=sequence.length();
            is_not_gap_vec.assign(L,0);
        }
        if (norm_by_ali==1)
        {
            Lali_vec.push_back(0);
            for (i=0;i<L;i++)
            {
                if (aa_list.find(sequence[i])!=string::npos)
                {
                    is_not_gap_vec[i]=1;
                    Lali_vec[Lali_vec.size()-1]+=1;
                }
                else
                    is_not_gap_vec[i]=0;
            }
            is_not_gap_mat.push_back(is_not_gap_vec);
        }
        aln.push_back(sequence);
        if (sequence.size()!=L)
        {
            cerr<<"ERROR! length not match for sequence "<<aln.size();
            exit(0);
        }
    }
    fp.close();
    is_not_gap_vec.clear();

    /* normalize target Nf */
    if (norm==0) target_Nf*=sqrt(L);
    else if (norm==1) target_Nf*=L;

    /* for simple seq number counting */
    int seq_num=aln.size();
    if (norm_by_ali==2)
    {
        if (norm==0) return 1.*seq_num/sqrt(L);
        else if (norm==1) return 1.*seq_num/L;
        else return seq_num;
    }

    /* compute seqID */
    float Nf=0;
    int Liden=0;
    vector<int> inv_seq_weight_vec(seq_num,1);
    for (m=0;m<seq_num;m++)
    {
        for (n=m+1;n<seq_num;n++)
        {
            Liden=0;
            for (i=0;i<L;i++)
            {
                Liden+=(aln[m][i]==aln[n][i] && (norm_by_ali==0
                    ||(is_not_gap_mat[m][i]*is_not_gap_mat[n][i])));
            }

            if (norm_by_ali==0)
            {
                inv_seq_weight_vec[m]+=(1.*Liden/L > id_cut);
                inv_seq_weight_vec[n]+=(1.*Liden/L > id_cut);
            }
            else
            {
                if (Lali_vec[m]*Lali_vec[n]>0)
                {
                    inv_seq_weight_vec[m]+=(1.*Liden/Lali_vec[m] > id_cut);
                    inv_seq_weight_vec[n]+=(1.*Liden/Lali_vec[n] > id_cut);
                }
            }
        }

        Nf+=1./inv_seq_weight_vec[m];
        if (target_Nf>0 && Nf>target_Nf) break;
    }

    /* normalize Nf */
    if (norm==0) Nf/=sqrt(L);
    else if (norm==1) Nf/=L;
    return Nf;
}


int main(int argc, char **argv)
{
    /* parse commad line argument */
    float id_cut=0.8; // defined by gremlin
    int norm=0; // 0 - L^0.5, 1 - L, 2 - no normalize
    float target_Nf=0;
    if(argc<2)
    {
        cerr<<docstring;
        return 0;
    }
    if (argc>2) id_cut=atof(argv[2]);
    if (argc>3) norm=atoi(argv[3]);
    if (argc>4) target_Nf=atof(argv[4]);
    cout<<calNf(argv[1],id_cut,norm,target_Nf)<<endl;
    return 0;
}

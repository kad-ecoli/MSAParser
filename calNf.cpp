const char* docstring=""
"calNf seq.aln 0.8 0\n"
"    calculate Nf, number of effective sequence defined in gremlin,\n"
"    using 0.8 (default) sequence identity cutoff.\n"
"\n"
"    The third optional argument is normalization methods:\n"
"       0 - normalized by L^0.5 (default)\n"
"       1 - normalized by L\n"
"       2 - do not normalize\n"
"      10 - calculate seqID by number of nongap positions\n";

#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <fstream>
#include <cmath>

using namespace std;

string aa_list="ACDEFGHIKLMNPQRSTVWY";

float calNf(const char *filename, float id_cut=0.8, int norm=0)
{
    int L=0; // alignment length
    int i,n,m; // position, sequence index, sequence index
    vector <int> is_not_gap_vec; // if a positions is an amino acid 
    vector<vector<int> > is_not_gap_mat; // if a position is an amino acid
    vector<float> cov_vec;  // coverage for each sequence
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
        if (norm_by_ali)
        {
            cov_vec.push_back(0);
            for (i=0;i<L;i++)
            {
                if (aa_list.find(sequence[i])!=string::npos)
                {
                    is_not_gap_vec[i]=1;
                    cov_vec[cov_vec.size()-1]+=1;
                }
                else
                    is_not_gap_vec[i]=0;
            }
            cov_vec[cov_vec.size()-1]/=(1.*L);
        }
        aln.push_back(sequence);
        is_not_gap_mat.push_back(is_not_gap_vec);
        if (sequence.size()!=L)
        {
            cerr<<"ERROR! length not match for sequence "<<aln.size();
            exit(0);
        }
    }
    fp.close();
    is_not_gap_vec.clear();

    /* compute seqID */
    int seq_num=aln.size();
    vector<float> seq_weight_vec(seq_num,0);
    vector<vector<float> > seqID_mat(seq_num,seq_weight_vec);
    for (m=0;m<seq_num-1;m++)
    {
        for (n=m+1;n<seq_num;n++)
        {
            for (i=0;i<L;i++)
            {
                seqID_mat[m][n]+=(aln[m][i]==aln[n][i] && (norm_by_ali==0
                    ||(is_not_gap_mat[m][i]*is_not_gap_mat[n][i])));
            }

            if (norm_by_ali==0)
            {
                seqID_mat[m][n]/=(1.*L);
                seqID_mat[n][m]=seqID_mat[m][n];
            }
            else
            {
                seqID_mat[n][m]=seqID_mat[m][n];
                if (cov_vec[m]*cov_vec[n]>0)
                {
                    seqID_mat[m][n]/=L*cov_vec[m];
                    seqID_mat[n][m]/=L*cov_vec[n];
                }
            }
        }
    }

    /* compute Nf */
    float Nf=0;
    for (m=0;m<seq_num;m++)
    {
        for (n=0;n<seq_num;n++)
        {
            seq_weight_vec[m]+=(m==n||seqID_mat[m][n]>id_cut);
        }
        seq_weight_vec[m]=1./seq_weight_vec[m];
        Nf+=seq_weight_vec[m];
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
    if(argc<2)
    {
        cerr<<docstring;
        return 0;
    }
    if (argc>2) id_cut=atof(argv[2]);
    if (argc>3) norm=atoi(argv[3]);
    cout<<calNf(argv[1],id_cut,norm)<<endl;
    return 0;
}

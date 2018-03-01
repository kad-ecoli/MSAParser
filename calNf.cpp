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
"\n"
"calNf seq.aln 0.8 0 128 1\n"
"    The fifth optional argument is whether to approximate Nf by number\n"
"    of non-redundant sequences. This algorithm has time complexity O(N*C)\n"
"    wheras the standard protocol has time complexity O(N*N), where N and\n"
"    are the number of non-redundant sequences\n"
;

#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <fstream>
#include <cmath>

using namespace std;

string aa_list="ACDEFGHIKLMNPQRSTVWY";

float calNf(const string infile="-", float id_cut=0.8, int norm=0,
    float target_Nf=0, int filter_aln_only=0)
{
    int L=0; // alignment length
    int i,n,m; // position, sequence index, sequence index
    vector <int> is_not_gap_vec; // if a positions is an amino acid 
    vector<vector<int> > is_not_gap_mat; // if a position is an amino acid
    vector<int> Lali_vec;  // length of aligned residues for each sequence
    int Lali=0,Liden=0;
    int Ldiff=0,Ldiff_m=0,Ldiff_n=0;
    int max_Ldiff=0,max_Ldiff_m=0,max_Ldiff_n=0;
    string sequence;
    vector <string> aln;
    int norm_by_ali=(norm/10);
    norm %=10;
    float seqID=0;

    /* read aln file */
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
            /* max num of different residues for two seq to be considered homolog */
            max_Ldiff=max_Ldiff_m=max_Ldiff_n=L*(1.-id_cut);
            is_not_gap_vec.assign(L,0);
            /* normalize target Nf */
            if (norm==0) target_Nf*=sqrt(L);
            else if (norm==1) target_Nf*=L;
        }
        if (norm_by_ali==1)
        {
            Lali=0;
            for (i=0;i<L;i++)
            {
                if (aa_list.find(sequence[i])!=string::npos)
                {
                    is_not_gap_vec[i]=1;
                    Lali+=1;
                }
                else is_not_gap_vec[i]=0;
            }
            if (filter_aln_only==0)
            {
                is_not_gap_mat.push_back(is_not_gap_vec);
                Lali_vec.push_back(Lali);
            }
        }
        if (filter_aln_only)
        {
            if (aln.size()>target_Nf) break;
            seqID=0;
            for (m=0;m<aln.size();m++)
            {
                Liden=0;
                for (i=0;i<L;i++)
                {
                    Liden+=(aln[m][i]==sequence[i] && (norm_by_ali==0
                        ||(is_not_gap_mat[m][i]*is_not_gap_vec[i])));
                }
                seqID=Liden/((norm_by_ali==0)?L:Lali);
                if (seqID>id_cut) break;
            }
            if (seqID>id_cut) break; // find homolog, do not add to aln
            if (norm_by_ali==1)
            {
                is_not_gap_mat.push_back(is_not_gap_vec);
                Lali_vec.push_back(Lali);
            }
        }
        if (sequence.size()!=L)
        {
            cerr<<"ERROR! length not match for sequence "<<aln.size();
            exit(0);
        }
        aln.push_back(sequence);
    }
    fp.close();
    is_not_gap_vec.clear();

    /* for simple seq number counting */
    int seq_num=aln.size();
    if (norm_by_ali==2 || filter_aln_only)
    {
        if (norm==0) return 1.*seq_num/sqrt(L);
        else if (norm==1) return 1.*seq_num/L;
        else return seq_num;
    }

    /* compute seqID */
    float Nf=0;
    vector<int> inv_seq_weight_vec(seq_num,1);
    for (m=0;m<seq_num;m++)
    {
        if (norm_by_ali) max_Ldiff_m=Lali_vec[m]*(1.-id_cut);
        for (n=m+1;n<seq_num;n++)
        {
            //Liden=0;
            Ldiff=Ldiff_m=Ldiff_n=0;
            if (norm_by_ali) max_Ldiff_n=Lali_vec[n]*(1.-id_cut);

            for (i=0;i<L;i++)
            {
                if (aln[m][i]==aln[n][i]) continue;
                if (norm_by_ali==0)
                {
                    Ldiff++;
                    if (Ldiff>max_Ldiff) break;
                }
                else
                {
                    Ldiff_m+=is_not_gap_mat[m][i];
                    Ldiff_n+=is_not_gap_mat[n][i];
                    if (Ldiff_m>max_Ldiff_m && Ldiff_n>max_Ldiff_n) break;
                }
                //Liden+=(aln[m][i]==aln[n][i] && (norm_by_ali==0
                    //||(is_not_gap_mat[m][i]*is_not_gap_mat[n][i])));
            }
            if (norm_by_ali==0) Ldiff_m=Ldiff_n=Ldiff;
            
            inv_seq_weight_vec[m]+=(Ldiff_m<=max_Ldiff_m);
            inv_seq_weight_vec[n]+=(Ldiff_n<=max_Ldiff_n);
            //if (norm_by_ali==0)
            //{
                //inv_seq_weight_vec[m]+=(1.*Liden/L > id_cut);
                //inv_seq_weight_vec[n]+=(1.*Liden/L > id_cut);
            //}
            //else if (Lali_vec[m]*Lali_vec[n]>0)
            //{
                //inv_seq_weight_vec[m]+=(1.*Liden/Lali_vec[m] > id_cut);
                //inv_seq_weight_vec[n]+=(1.*Liden/Lali_vec[n] > id_cut);
            //}
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
    int filter_aln_only=0; // whether approximate Nf by the number of
                           // non-redundant sequences
    if(argc<2)
    {
        cerr<<docstring;
        return 0;
    }
    string infile=argv[1];
    if (argc>2) id_cut=atof(argv[2]);
    if (id_cut>1) id_cut/=100.;
    if (argc>3) norm=atoi(argv[3]);
    if (argc>4) target_Nf=atof(argv[4]);
    if (argc>5) filter_aln_only=atoi(argv[5]);
    cout<<calNf(infile,id_cut,norm,target_Nf,filter_aln_only)<<endl;
    return 0;
}

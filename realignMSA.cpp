const char* docstring=""
"realignMSA uniclust.fasta metaclust.fasta\n"
"    realign metaclust.fasta to uniclust.fasta according to the match\n"
"    state of the first sequence in uniclust.fasta\n"
"\n"
"realignMSA uniclust.fasta metaclust.fasta realign.fasta\n"
"    output the result to realign,fasta\n"
;

#include <iostream>
#include <vector>
#include <string>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <algorithm>

using namespace std;

/* parse match state of first sequence in input fasta file, 
 * return the sequence length and 
 * match_state_vec: index of match state positions in query [a-zA-Z]
 * del_state_vec: index of positions that should be deleted from homolog
 */
int parse_match_state_of_first_fasta(
     vector <int> & match_state_vec,  // [A-Z] position in [a-zA-Z] sequence
     vector <int> & del_state_vec, // [-] position in [-A-Z] sequence
     const char *filename)
{
    ifstream fp;
    fp.open(filename,ios::in);
    string sequence; // including [-a-zA-Z] excluding [.]
    string line;
    int i;
    int L=0; // protein length, including [a-zA-Z] excluding [-.]

    /* read [-a-zA-Z] sequence */
    while (fp.good())
    {
        getline(fp,line);
        if (line.length()==0) continue;
        if (line[0]=='>')
        {
            if (sequence.length()==0) continue;
            break;
        }

        for (i=0;i<line.length();i++)
        {
            if (line[i]!='.') 
            {
                sequence+=line[i];
                L+=(line[i]!='-');
            }
        }
    }
    fp.close();

    /* read [-] position in [-A-Z] sequence into del_state_vec */
    /* read [A-Z] positions in [a-zA-Z] sequence into match_state_vec */
    int upper_idx=-1;
    int alphabet_idx=-1;
    for (i=0;i<sequence.length();i++)
    {
        upper_idx+=(sequence[i]=='-' || isupper(sequence[i]));
        alphabet_idx+=(sequence[i]!='-');
        if (sequence[i]=='-') 
            del_state_vec.push_back(upper_idx);
        else if (isupper(sequence[i]))
            match_state_vec.push_back(alphabet_idx);
    }

    sequence.clear();
    return L;
}

/* realign sequence based on match states*/
string realignSeq(const vector<int> &match_state_vec,int L,string sequence)
{
    string realigned_seq="";
    int i;
    for (i=0;i<L;i++) realigned_seq+='-';
    for (i=0;i<sequence.length();i++)
    {
        realigned_seq[match_state_vec[i]]=sequence[i];
    }
    return realigned_seq;
}

/* insert gaps into input fasta according to match state of query */
int realignMSA(const vector<int>& match_state_vec,
    const vector<int>& del_state_vec, int L,
     string infile="-", string outfile="-")
{
    ifstream fp_in;
    ofstream fp_out;
    if (infile!="-") fp_in.open(infile.c_str(),ios::in);
    if (outfile!="-") fp_out.open(outfile.c_str(),ofstream::out);
    string sequence,line;
    int nseqs=0;
    int i;
    while ((infile!="-")?fp_in.good():cin.good())
    {
        if (infile!="-") getline(fp_in,line);
        else getline(cin,line);

        if (line.length()==0) continue;
        if (line[0]=='>')
        {
            if (sequence.length()>0)
            {
                if (outfile!="-")
                    fp_out<<realignSeq(match_state_vec,L,sequence)<<endl;
                else
                    cout<<realignSeq(match_state_vec,L,sequence)<<endl;
            }
            if (outfile!="-") fp_out<<line<<endl;
            else cout<<line<<endl;
            sequence.clear();
            nseqs++;
        }
        else
        {
            if (del_state_vec.size())
            {
                // not efficient way to check if i is in del_state_vec
                // only read positions not corresponding to [-] in query
                for (i=0;i<line.length();i++)
                {
                    if (find(del_state_vec.begin(), del_state_vec.end(), i
                        ) != del_state_vec.end())
                        sequence+=line[i];
                }
            }
            else
                sequence+=line;
        }
    }
    fp_in.close();
    if (outfile!="-") fp_out<<realignSeq(match_state_vec,L,sequence)<<endl;
    else cout<<realignSeq(match_state_vec,L,sequence)<<endl;
    fp_out.close();
    sequence.clear();
    return nseqs;
}

int main(int argc, char **argv)
{
    /* parse commad line argument */
    if(argc<3)
    {
        cerr<<docstring;
        return 0;
    }
    /* read query sequence */
    vector <int> match_state_vec;  // based on query sequence residue index
    vector <int> del_state_vec; // based on (gapped) homolog position index
    int L=parse_match_state_of_first_fasta(
        match_state_vec,del_state_vec,argv[1]);
    /* re-align sequence profile */
    string outfile=(argc<=3)?"-":argv[3];
    realignMSA(match_state_vec,del_state_vec,L,argv[2],outfile);
    return 0;
}

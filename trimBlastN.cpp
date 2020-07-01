const char* docstring=""
"trimBlastN blastnt.db blastnt.tab  L > blastnt.trim.fasta\n"
"    trim fasta database blastnt.db according to blastnt.tab,\n"
"    so that for a hit in blastnt.db, the 5' and 3' termini are trimmed\n"
"    with up to L residues at each side flanking the aligned region.\n"
"    blastnt.tab must be generated by\n"
"    $ blastn -outfmt '6 saccver sstart send'\n"
;

#include <iostream>
#include <vector>
#include <string>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <sstream>

using namespace std;

/* split a long string into vectors by whitespace 
 * line          - input string
 * line_vec      - output vector 
 * delimiter     - delimiter */
void split(const string &line, vector<string> &line_vec, const char delimiter='\t')
{
    bool within_word = false;
    for (size_t pos=0;pos<line.size();pos++)
    {
        if (line[pos]==delimiter)
        {
            within_word = false;
            continue;
        }
        if (!within_word)
        {
            within_word = true;
            line_vec.push_back("");
        }
        line_vec.back()+=line[pos];
    }
}

void reverse_complement(const string &watson, string &crick)
{
    size_t Lfrag=watson.size();
    char C;
    for (size_t i=Lfrag-1;i>=0;i--)
    {
        switch (watson[i])
        {
            case 'A': C='T'; break;
            case 'I': C='T'; break;
            case 'C': C='G'; break;
            case 'G': C='C'; break;
            case 'T': C='A'; break;
            case 'U': C='A'; break;
            case 'W': C='W'; break;
            case 'S': C='S'; break;
            case 'M': C='K'; break;
            case 'K': C='M'; break;
            case 'R': C='Y'; break;
            case 'Y': C='R'; break;
            case 'B': C='V'; break;
            case 'D': C='H'; break;
            case 'H': C='D'; break;
            case 'V': C='B'; break;
            case 'N': C='N'; break;
            case 'Z': C='Z'; break;

            case 'a': C='t'; break;
            case 'i': C='t'; break;
            case 'c': C='g'; break;
            case 'g': C='c'; break;
            case 't': C='a'; break;
            case 'u': C='a'; break;
            case 'w': C='w'; break;
            case 's': C='s'; break;
            case 'm': C='k'; break;
            case 'k': C='m'; break;
            case 'r': C='y'; break;
            case 'y': C='r'; break;
            case 'b': C='v'; break;
            case 'd': C='h'; break;
            case 'h': C='d'; break;
            case 'v': C='b'; break;
            case 'n': C='n'; break;
            case 'z': C='z'; break;

            case '*': C='*'; break;
            case '-': C='-'; break;
            case '.': C='.'; break;
            default:  C='N'; break;
        }
        crick[Lfrag-1-i]=C;
        if (i==0) break;
    }
    return;
}

void getSeqTxt(const vector<string>&acc_list,
    const vector<size_t>&from_list, const vector<size_t>&to_list,
    const int L, const string &header, const string &sequence, string &txt)
{
    size_t n,from,to;
    stringstream ss;
    char fr;
    string fragment;
    for (n=0;n<acc_list.size();n++)
    {
        if (acc_list[n]!=header) continue;
        if (from_list[n]<to_list[n])
        {
            from=from_list[n];
            to  =to_list[n];
            fr  ='f';
        }
        else
        {
            from=to_list[n];
            to  =from_list[n];
            fr  ='r';
        }
        if (from<L+1) from=1;
        else from-=L;
        to  +=L;
        fragment=sequence.substr(from-1,to-from+1);
        if (fr=='r') reverse_complement(
                 sequence.substr(from-1,to-from+1),fragment);
        ss<<'>'<<header<<'_'<<from<<'_'<<to<<'_'<<fr<<'\n'
            <<fragment<<'\n';
        txt+=ss.str();
        ss.str(string());
        fragment.clear();
    }
    return;
}

void trimBlastN(const string indbfile="-", const string intabfile="-",
    const int L=0, const string outfile="-")
{
    /* read tab file */
    vector<string> acc_list;
    vector<size_t> from_list;
    vector<size_t> to_list;
    string line;
    vector<string>line_vec;
    size_t i;
    ifstream fp_in;
    if (intabfile!="-") fp_in.open(intabfile.c_str(),ios::in);
    while ((intabfile!="-")?fp_in.good():cin.good())
    {
        if (intabfile!="-") getline(fp_in,line);
        else                getline(cin,line);
        if (line.size()==0) continue;
        split(line,line_vec,'\t');
        if (line_vec.size()<=2)
        {
            cerr<<"FATAL ERROR! Less than 3 columns in "<<intabfile<<endl;
            return;
        }
        acc_list.push_back(line_vec[0]);
        from_list.push_back(atoi(line_vec[1].c_str()));
        to_list.push_back(atoi(line_vec[2].c_str()));
        for (i=0;i<line_vec.size();i++) line_vec[i].clear();
        line_vec.clear();
    }
    if (intabfile!="-") fp_in.close();

    //for (i=0;i<acc_list.size();i++) cout<<"["<<i<<"] "
        //<<acc_list[i]<<' '<<from_list[i]<<'-'<<to_list[i]<<endl;

    /* read db file */
    ofstream fp_out;
    if (indbfile!="-") fp_in.open(indbfile.c_str(),ios::in);
    if (outfile!="-") fp_out.open(outfile.c_str(),ofstream::out);
    string sequence,header,txt;
    while ((indbfile!="-")?fp_in.good():cin.good())
    {
        if (indbfile!="-") getline(fp_in,line);
        else               getline(cin,line);

        if (line.length()==0) continue;
        if (line[0]=='>')
        {
            if (sequence.length()>0)
            {
                getSeqTxt(acc_list, from_list, to_list, L, 
                        header, sequence, txt);
                if (outfile!="-") fp_out<<txt;
                else                cout<<txt;
                txt.clear();
            }
            sequence.clear();
            split(line, line_vec, ' ');
            header=line_vec[0].substr(1);
            for (i=0;i<line_vec.size();i++) line_vec[i].clear();
            line_vec.clear();
        }
        else sequence+=line;
    }
    fp_in.close();
    getSeqTxt(acc_list, from_list, to_list, L, header, sequence, txt);
    if (outfile!="-") fp_out<<txt;
    else                cout<<txt<<flush;
    fp_out.close();
    
    /* clean up */
    from_list.clear();
    to_list.clear();
    for (i=0;i<acc_list.size();i++) acc_list[i].clear();
    acc_list.clear();
    sequence.clear();
    header.clear();
    line.clear();
    txt.clear();
    return;
}

int main(int argc, char **argv)
{
    /* parse commad line argument */
    if(argc<3)
    {
        cerr<<docstring;
        return 0;
    }
    string indbfile =argv[1];
    string intabfile=argv[2];
    int    L        =(argc<=3)?0:atoi(argv[3]);
    string outfile  =(argc<=4)?"-":argv[4];
    trimBlastN(indbfile,intabfile,L,outfile);
    return 0;
}

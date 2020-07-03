const char* docstring=""
"afa2blasttab qnhmmer.afa L nhmmer.tab\n"
"    Convert qnhmmer alignment (post-processed with a3m2ma and\n"
"    RemoveNonQueryPosition) to a hit table similar to those generated with\n"
"    $ blastn -outfmt '6 saccver sstart send'\n"
"    For a hit in qnhmmer.db, the 5' and 3' termini are trimmed with\n"
"    with up to L residues at each side flanking the aligned region.\n"
"    Output table to nhmmer.tab\n"
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

void afa2blasttab(const string infile="-", const int L=0, const string outfile="-")
{
    stringstream ss;
    string line,name;
    vector<string> line_vec;
    size_t i,from,to,leftgap,rightgap;

    /* read input file */
    ifstream fp_in;
    if (infile!="-") fp_in.open(infile.c_str(),ios::in);
    while ((infile!="-")?fp_in.good():cin.good())
    {
        if (infile!="-") getline(fp_in,line);
        else             getline(cin,line);
        if (line.size()==0) continue;
        if (line[0]=='>')
        {
            from=to=0;
            split(line,line_vec,' ');
            line=line_vec[0];
            line_vec.clear();
            split(line,line_vec,'/');
            if (line_vec.size()!=2)
            {
                cerr<<"WARNING! sequence name must be in the format of >name/from-to. skipping "<<line<<endl;
                line_vec.clear();
                continue;
            }
            name=line_vec[0].substr(1);
            line=line_vec[1];
            line_vec.clear();
            split(line,line_vec,'-');
            if (line_vec.size()!=2)
            {
                cerr<<"WARNING! sequence name must be in the format of >name/from-to. skipping "<<line<<endl;
                line_vec.clear();
                continue;
            }
            from=atoi(line_vec[0].c_str());
            to  =atoi(line_vec[1].c_str());
            line_vec.clear();
            continue;
        }
        if (from==0 && to==0) continue;
        leftgap=0;
        rightgap=0;
        if (L>0)
        {
            for (i=0;i<line.size();i++)
            {
                if (line[i]!='-') break;
                leftgap++;
                if (leftgap>=L) break;
            }
            if (from<=leftgap) leftgap=from-1;
            for (i=line.size()-1;i>=0;i--)
            {
                if (line[i]!='-') break;
                rightgap++;
                if (rightgap>=L) break;
            }
        }
        ss<<name<<'\t'<<from-leftgap<<'\t'<<to+rightgap<<endl;
    }
    if (infile!="-") fp_in.close();

    /* read db file */
    if (outfile=="-")
    {
        cout<<ss.str()<<flush;
    }
    else
    {
        ofstream fp_out;
        fp_out.open(outfile.c_str(),ofstream::out);
        fp_out<<ss.str();
        fp_out.close();
    }
    
    /* clean up */
    ss.str(string());
    line.clear();
    name.clear();
    vector<string> ().swap(line_vec);
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
    string infile =argv[1];
    int    L      =atoi(argv[2]);
    string outfile=(argc<=3)?"-":argv[3];
    afa2blasttab(infile,L,outfile);
    return 0;
}

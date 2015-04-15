/*
*    This is programm Graphfilter.
*
*    Graphfilter is free software: you can redistribute it and/or modify
*    it under the terms of the GNU General Public License as published by
*    the Free Software Foundation, either version 3 of the License, or
*    (at your option) any later version.
*
*    Graphfilter is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*    GNU General Public License for more details.
*
*    You should have received a copy of the GNU General Public License
*    along with Graphfilter.  If not, see <http://www.gnu.org/licenses/>.
*
*/
/** Graphfilter
Copyright (C) Miguel A. Santos, HSC, Toronto, 2009-2010.
Licensed under GPL version 3 a later. (see http://www.gnu.org/copyleft/gpl.html )
*/
//Compile: g++ -lgsl -lgslcblas -o graphfilter graphfilter.cc

//	VERSION 1.5.2

#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <set>
#include <vector>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>

using namespace std;

enum TRANS {mident, madd, mdiv, mmul, mlog, mlogneg, mexp, mpow , msvd , msval , mprunesval,mprunesvalB, mprunesvalA, mprtsvd, mfillmissing, Exp, mrandsample, mswapnames };

void license(){
	cout<<\
"Graphfilter V1.5.2\n"<<\
"Copyright (C) Miguel A. Santos, HSC, Toronto, 2009-2010.\n"<<\
"Licensed under GPL version 3 a later. (see http://www.gnu.org/copyleft/gpl.html )\n"<<\
	endl;
}
void usage(){
	cout<<"Usage: graphfilter [options] [graph.file]"<<endl;
	cout<<\
"\nGraphfilter allows you to perform different transformations on the edges\n\
of a graph, like taking the log/-log in any base, rasing to a given power,\n\
taking the exponential, etc, or filtering out some eigenvalues.\n\
\n\
Graphs can be read from any table. The only condition is that the nodes\n\
are expected on the first two columns of the table. The location of the\n\
edge values can be specified with option -c.\n"<<\
"\nBy extension, it also works with any (non-square) matrix, given that the input \
abides that condition (Will be available again in next version).\n\
Any matrix element A(i,j) is akin to an edge between nodes `i' and `j'.\n"<<\
	"\nUse option -h, --help for help.\n"<<\
	endl;
	cout<<"\nOptions:\n  General:\n"<<\
	"\t-h, --help Prints this help\n"<<\
	"\t-q Quiet mode: do not print any comment lines.\n"<<\
	"\t-Q Censor mode: Be quiet and also skip any comment lines from input file. Comment lines are those starting with `#'\n"<<\
	"\t--verbose Prints out more detailed info on what's going on.\n"<<\
	"  NO-transformations:\n"<<\
	"\t-c n Specifies the edges to be at the n-th column\n"<<\
	"\t-e, --Eigen-values Get the Eigen values of the matrix\n"<<\
	"\t-X, --non-symmetric Do _not_ assume provided matrix is symmetric. By default we assume it is.\n"<<\
	"\t-n, --precision n Set the number of decimals to use. Does not apply for Diagonalization.\n"<<\
	"\t-d, --print-evd Print matrices U and V from the Eigen Value Decomposition (Diagonalization) of the given graph G: G=UEU^t\n"<<\
	"\t-G, --print-no-diagonal When dealing with Diagonalization, do no print out the diagonal\n"<<\
	"\t-T, --print-triangular When dealing with E, just print out the upper/lower triangular part\n"<<\
	"\t-g, --scientific Use scientific notation for the output. Does not apply for Diagonalization.\n"<<\
	"\t-s, --swap-names tab-file Swap node names according to tab-file. The latter is a 2-column hash table of equivalent names.\n"<<\
	"  Transformations:\n"<<\
	"\t-A, --add w Adds w to all edges\n"<<\
	"\t-B, --base b  Specifies the base of logarithms to use\n"<<\
	"\t-C, --diag w (Re)Sets all diagonal values to w.\n"<<\
	"\t-D, --div w Divides all edges by w\n"<<\
	"\t-E, --exp  Takes the exp of the edges\n"<<\
	"\t  , --ExpW , --exponentiation  [-beta beta] Takes the exp of the graph W as E^(beta*W). Default: beta=1\n"<<\
	"\t-F, --fill-missing [-v value] Fill missing edges with weight=value (Default: value=0).\n"<<\
	"\t-I, --identity  Do nothing but printing out the edges. Default action.\n"<<\
	"\t-L, --log  Takes the log (base e) of the edges\n"<<\
	"\t-N, --log-neg  Takes the negative log (base e) of the edges\n"<<\
	"\t-M, --mul w Multiplies all edges by w\n"<<\
	"\t-P, --pow n Takes the n-power of the edges\n"<<\
	"\t-p, --prune-Eigen-values [below|above] sv Removes all Eigen values below/above sv. Default: those below sv.\n"<<\
	"\t-R, --random-sample n [-s seed]  Cull n edges at random from the given graph\n"<<\
	"\t-S, --symmetric Assume provided matrix is symmetric (Default).\n"<<\
	endl;
	cout<<\
"Examples:\n\
\tcat graph | graphfilter -c 11 --log-neg --base 10 > graph-log10\n\
\n\tConsiders column 11 of file graph as the edge values and takes the -log of them in base 10.\n\
\n\tcat graph | graphfilter -c 11 --mult 3.5 | graphfilter --exp > graph-exp3.5\n\
\n\tThe file graph-exp3.5 will contain the exponential of original edges to the power of 3.5, i.e, exp(3.5*edge)\n\
\tBy default graphfilter expects the edges in column 3 (The first column is 1).\n"<<\
"\n\tFilter out all Eigen values below 100 and get the new matrix values:      \n\
\tgraphfilter graph2 --prune-Eigen-values 100\n"<<\
	endl;
	license();
	exit(1);
}

int main(int argc, char* argv[]){
	argc--,argv++;
	/*if(argc<1){
		usage();
		exit(1);
	}*/
	int col=3;
	TRANS transf=mident;
	string transftype="Identity";
	string addcomment="";
	bool TRANSEDGE1=true;
	bool QUIET=false;
	bool VERBOSE=false;
	bool CENSOR=false;
	bool SCIENTFMT=false;
	bool SETPRECIS=false;
	bool SETDIAGVAL=false;
	bool SYMMETRIC=true;
	bool FILLMISSING=false;
	bool PRT_E_GRAPH_FULL=true;
	bool PRT_E_GRAPH_TRIANG=!PRT_E_GRAPH_FULL;
	bool PRT_E_GRAPH_NODIAG=!PRT_E_GRAPH_FULL;
	char* fn="/dev/stdin";
	char* tabfile=NULL;
	double myarg;
	double missingvalue=0.0;
	double _beta=1.0;
	double mdiagval;
	double logbase=1.0;
	unsigned int randSeed=9876543;
	unsigned int sample_size=0;
	int precision=10;
	map<string, string > forward_tabfile;
	map<string, string > reverse_tabfile;
	//if(!QUIET)cout<<"#Command: "<<argv<<endl;
	while(argc>0){
///General
		if(strcmp(argv[0],"-h")==0||strcmp(*argv,"--help")==0){
			usage();
		}
		else if(strcmp(argv[0],"-q")==0){
			QUIET=true;
			argc--,argv++;
		}
		else if(strcmp(argv[0],"-Q")==0){
			CENSOR=true;
			QUIET=true;
			argc--,argv++;
		}
		else if(strcmp(argv[0],"--verbose")==0){
			VERBOSE=true;
			QUIET=true;
			argc--,argv++;
		}
		else if(strcmp(*argv,"-c")==0){
			argc--,argv++;
			if(argc<1)usage();
			col=atoi(argv[0]);
			argc--,argv++;
		}
		else if(strcmp(argv[0],"--base")==0||strcmp(argv[0],"-B")==0){
			argc--,argv++;
			if(argc<1)usage();
			logbase=log(atof(*argv));
			stringstream ss;
			ss<<"base="<<*argv;
			ss>>addcomment;
			argc--,argv++;
		}
		else if(strcmp(argv[0],"--Eigen-values")==0||strcmp(argv[0],"-e")==0){
			transf=msval;
			TRANSEDGE1=false;
			transftype="Eigen-values";
			argc--,argv++;
		}
		else if(strcmp(argv[0],"--precision")==0||strcmp(argv[0],"-n")==0){
			argc--,argv++;
			if(argc<1)usage();
			precision=atoi(*argv);
			argc--,argv++;
			SETPRECIS=true;
		}
		else if(strcmp(argv[0],"--print-evd")==0||strcmp(argv[0],"-d")==0){
			transf=mprtsvd;
			TRANSEDGE1=false;
			transftype="print-EVD";
			argc--,argv++;
		}
		else if(strcmp(argv[0],"--print-no-diagonal")==0||strcmp(argv[0],"-G")==0){
			PRT_E_GRAPH_FULL=false;
			PRT_E_GRAPH_NODIAG=true;
			argc--,argv++;
		}
		else if(strcmp(argv[0],"--print-triangular")==0||strcmp(argv[0],"-T")==0){
			PRT_E_GRAPH_FULL=false;
			PRT_E_GRAPH_TRIANG=true;
			argc--,argv++;
		}
		else if(strcmp(argv[0],"--scientific")==0||strcmp(argv[0],"-g")==0){
			argc--,argv++;
			SCIENTFMT=true;
		}
		else if(strcmp(argv[0],"--swap-names")==0||strcmp(argv[0],"-s")==0){
			transf=mswapnames;
			argc--,argv++;
			if(argc<1)usage();
			tabfile=argv[0];
			ifstream is(argv[0]);	
			if(!is){
			        cout<<"ERROR: Cannot open file "<<argv[0]<<endl;
			        exit(1);
			}
			argc--,argv++;
			string a,b;
			pair< string, string > p;
			while(is>>a && is>>b){
				if(a.substr(0,1).compare("#")==0){
					string line;
					getline(is,line);
					continue;
				}
				p = make_pair (a,b);
				forward_tabfile.insert(p);
				p = make_pair (b,a);
				reverse_tabfile.insert(p);
			}
			if(forward_tabfile.size()!=reverse_tabfile.size()){
				cout<<"ERROR: swap-names : forward_tabfile.size()!=reverse_tabfile.size() "<<endl;
				exit(1);
			}
			if(forward_tabfile.size()==0){
				cout<<"ERROR: swap-names : tabfile.size==0 "<<endl;
				exit(1);
			}
			if(!QUIET)cout<<"#Tabfile found with "<<forward_tabfile.size()<<" entries"<<endl;
		}
///Transformations
		else if(strcmp(argv[0],"--add")==0||strcmp(argv[0],"-A")==0){
			transf=madd;
			argc--,argv++;
			if(argc<1)usage();
			myarg=atof(argv[0]);
			argc--,argv++;
			stringstream ss;
			transftype="adding  ";
			ss<<transftype<<myarg;
			transftype=ss.str();
		}
		else if(strcmp(*argv,"--diag")==0||strcmp(argv[0],"-C")==0){
			argc--,argv++;
			if(argc<1)usage();
			mdiagval=atof(argv[0]);
			SETDIAGVAL=true;
			argc--,argv++;
		}
		else if(strcmp(argv[0],"--div")==0||strcmp(argv[0],"-D")==0){
			transf=mdiv;
			argc--,argv++;
			if(argc<1)usage();
			myarg=atof(argv[0]);
			argc--,argv++;
			stringstream ss;
			transftype="divividing  by  ";
			ss<<transftype<<myarg;
			transftype=ss.str();
		}
		else if(strcmp(argv[0],"--exp")==0||strcmp(argv[0],"-E")==0){
			transf=mexp;
			transftype="exp";
			argc--,argv++;
		}
		else if(strcmp(argv[0],"--ExpW")==0||strcmp(argv[0],"--exponentiation")==0){
			transf=Exp;
			transftype="e^";
			TRANSEDGE1=false;
			argc--,argv++;
			if(argc>0&&strcmp(*argv,"-beta")==0){
				argc--,argv++;
				_beta=atof(argv[0]);
			}
		}
		else if(strcmp(argv[0],"--fill-missing")==0||strcmp(argv[0],"-F")==0){
			transf=mfillmissing;
			transftype="Fill missing edge with weight=";
			FILLMISSING=true;;
			argc--,argv++;
			if(argc>0){
				if(strcmp(*argv,"-v")==0){
					argc--,argv++;
					if(argc<1)usage();
					missingvalue=atof(argv[0]);
					argc--,argv++;
				}
			}
			stringstream ss;
			ss<<transftype<<" "<<missingvalue;
			transftype=ss.str();
		}
		else if(strcmp(argv[0],"--identity")==0||strcmp(argv[0],"-I")==0){
			argc--,argv++;
		}
		else if(strcmp(argv[0],"--log")==0||strcmp(argv[0],"-L")==0){
			transf=mlog;
			transftype="ln";
			argc--,argv++;
		}
		else if(strcmp(argv[0],"--log-neg")==0||strcmp(argv[0],"-N")==0){
			transf=mlogneg;
			transftype="-ln";
			argc--,argv++;
		}
		else if(strcmp(argv[0],"--mul")==0||strcmp(argv[0],"-M")==0){
			transf=mmul;
			argc--,argv++;
			if(argc<1)usage();
			myarg=atof(argv[0]);
			argc--,argv++;
			stringstream ss;
			transftype="multiplying  by  ";
			ss<<transftype<<myarg;
			transftype=ss.str();
		}
		else if(strcmp(argv[0],"--non-symmetric")==0||strcmp(argv[0],"-X")==0){
			SYMMETRIC=false;;
			argc--,argv++;
		}
		else if(strcmp(argv[0],"--pow")==0||strcmp(argv[0],"-P")==0){
			transf=mpow;
			argc--,argv++;
			if(argc<1)usage();
			myarg=atof(argv[0]);
			argc--,argv++;
			stringstream ss;
			transftype="x^";
			ss<<transftype<<myarg;
			transftype=ss.str();
		}
		else if(strcmp(argv[0],"--prune-Eigen-values")==0||strcmp(argv[0],"-p")==0){
			transf=mprunesval;
			TRANSEDGE1=false;
			transftype="prune-Eigen-values-below:\t";
			argc--,argv++;
			if(argc<1)usage();
			if(strcmp(*argv,"-below")==0||strcmp(*argv,"below")==0){
				argc--,argv++;
				if(argc<1)usage();
				transf=mprunesvalB;
			}
			else if(strcmp(*argv,"-above")==0||strcmp(*argv,"above")==0){
				argc--,argv++;
				if(argc<1)usage();
				transftype="prune-Eigen-values-above:\t";
				transf=mprunesvalA;
			}
			myarg=atof(argv[0]);
			argc--,argv++;
			stringstream ss;
			ss<<transftype<<myarg;
			transftype=ss.str();
		}
		else if(strcmp(argv[0],"--random-sample n [-s seed] ")==0||strcmp(argv[0],"-R")==0){
			transf=mrandsample;
			TRANSEDGE1=false;
			transftype="random-sampling-edges-with-seed:\t";
			argc--,argv++;
			if(argc<1)usage();
			sample_size=atoi(argv[0]);
			argc--,argv++;
			if(argc>0 && strcmp(argv[0],"-s")==0){
				argc--,argv++;
				randSeed=atoi(argv[0]);
				argc--,argv++;
			}
			stringstream ss;
			ss<<transftype<<randSeed;
			transftype=ss.str();
		}
		else if(strcmp(argv[0],"--symmetric")==0||strcmp(argv[0],"-S")==0){
			SYMMETRIC=true;;
			argc--,argv++;
		}
		else{
			fn=argv[0];
			argc--,argv++;
		}
	}
	if(!QUIET){
		cout<<"#Transforming edges of graph "<<fn<<endl;
		if(SETDIAGVAL)cout<<"#Setting diag values to: "<<mdiagval<<endl;
		cout<<"#Type: "<<transftype<<" "<<addcomment<<"\n#"<<endl;
	}
	argc--,argv++;
	ifstream is(fn);
	if(!is){
		cout<<"Error: Could not open file "<<fn<<endl;
		exit(1);
	}
	string ln;
	//size_t NumberSeenEdges=0;
	map< pair<string,string >, double > inmatrix;
	//We need this for random sampling the graph
	vector< pair<string,string > > _graph_keys;
	set<string> elements;
	while(getline(is,ln)){
		if(ln.substr(0,1).compare("#")==0){
			if(!CENSOR)cout<<ln<<endl;
			continue;
		}
		stringstream ss;
		ss<<ln;
		string a,b,c;
		double f;
		ss>>a ; ss>>b; 
		if(a.empty()||b.empty()){
			//Ignore blank lines
			continue;
		}
		int shft=col-2;
		while(shft-->1)ss>>c;
		ss>>f;
		elements.insert(a);
		elements.insert(b);
		//IF task is culling edge at random, we dont care about diagonal values or symmetry. We just deal with what
		//the user provided as input as it is given.
		if(SETDIAGVAL && transf!=mrandsample){
			if(a.compare(b)==0)
				f=mdiagval;
			else{
				inmatrix.insert(pair< pair<string,string>, double> (pair<string,string>(a,a),mdiagval));
				inmatrix.insert(pair< pair<string,string>, double> (pair<string,string>(b,b),mdiagval));
			}
		}
		inmatrix.insert(pair< pair<string,string>, double> (pair<string,string>(a,b),f));
		if(SYMMETRIC&&a.compare(b)!=0 && transf!=mrandsample){
			inmatrix.insert(pair< pair<string,string>, double> (pair<string,string>(b,a),f));
		}
		if(transf==mrandsample) _graph_keys.push_back(pair<string,string>(a,b));
		switch(transf){
			case mlogneg:
			case mlog:
				if(f<=0){
					cout<<"ERROR: taking log of edge<0: "<<a<<"\t"<<b<<"\t"<<f<<endl;
					exit(1);
				}
				switch(transf){
					case mlog:
						f=log(f)/logbase;
						break;
					case mlogneg:
						f=-log(f)/logbase;
						break;
				}
				break;
			case mexp:
				f=exp(f);
				break;
			case mpow:
				f=pow(f,myarg);
				break;
			case madd:
				f+=myarg;
				break;
			case mmul:
				f*=myarg;
				break;
			case mdiv:
				f/=myarg;
				break;
			case mswapnames:{
				map<string, string >::iterator it_tab_a,it_tab_b;
				it_tab_a=forward_tabfile.find(a);
				it_tab_b=forward_tabfile.find(b);
				if(it_tab_a==forward_tabfile.end()) it_tab_a=reverse_tabfile.find(a);
				if(it_tab_b==forward_tabfile.end()) it_tab_b=reverse_tabfile.find(b);
				if(it_tab_a==reverse_tabfile.end()||it_tab_b==reverse_tabfile.end()){
					cout<<"ERROR: swap-names : No translation found for some nodes ("<<a<<" / "<<b<<")"<<endl;
					exit(1);
				}
				a=it_tab_a->second;
				b=it_tab_b->second;
				}
				break;
			case mident:
				break;
		}
		if(TRANSEDGE1&&!FILLMISSING){
			cout<<a<<"\t"<<b<<"\t";
			if(SCIENTFMT)cout.setf(ios::scientific,ios::floatfield);
			if(SETPRECIS)cout.precision(precision);
			cout<<f<<endl;
			//cout<<a<<"\t"<<b<<"\t"<<f<<endl;
		}
	}
	if(FILLMISSING){
		pair<string, string> p,pi;
		map<pair< string, string>, double >::iterator failed=inmatrix.end();
		for(set<string>::iterator a=elements.begin();a!=elements.end();a++){
			for(set<string>::iterator b=elements.begin();b!=elements.end();b++){
				p= make_pair(*a,*b);
				//pi= make_pair(*b,*a);
				//if(inmatrix.find(p)==inmatrix.end()&&inmatrix.find(pi)==inmatrix.end())inmatrix[p]=missingvalue;
				if(inmatrix.find(p)==failed){
					inmatrix[p]=missingvalue;
					if(SYMMETRIC){
						inmatrix[pi]=missingvalue;
					}
				}
				cout<<*a<<"\t"<<*b<<"\t";
				if(SCIENTFMT)cout.setf(ios::scientific,ios::floatfield);
				if(SETPRECIS)cout.precision(precision);
				cout<<inmatrix[p]<<endl;
			}
		}
		exit(0);
	}
	if(!QUIET){
		size_t efNedges=inmatrix.size();
		if(SETDIAGVAL) efNedges-=elements.size();
		if(SYMMETRIC) efNedges/=2;
		cout<<"#elements= "<<elements.size()<<" #edges= "<<efNedges<<endl;
	}
	if(!TRANSEDGE1){
		switch(transf){
			case mrandsample:
				cout<<"#BeginRandomSamplingEdges"<<endl;
				srand(randSeed);
				for(int i=0;i<sample_size;i++){
					unsigned int randedge = rand() % inmatrix.size() ;
					//cout<<"#Chosen edge "<<randedge<<" out of "<<sample_size<<endl;
					pair<string, string> p = _graph_keys[randedge];
					cout<<p.first<<"\t"<<p.second<<"\t"<<inmatrix[p]<<endl;
					//cout<<_graph_keys[randedge]<<"\t"<<inmatrix[_graph_keys[randedge]]<<endl;
				}
				cout<<"#EndRandomSamplingEdges"<<endl;
				exit(0);
				break;
		}
		size_t msize=elements.size();
		if(!QUIET)cout<<"#Inmatrix: elements="<<msize<<" pairs="<<inmatrix.size()<<endl;
		gsl_matrix* m=gsl_matrix_calloc(msize,msize);
		size_t i,j;
		i=0;
		map<pair<int,int>,pair<string,string> > intedgeTostredge;
		for(set<string>::iterator ait=elements.begin();ait!=elements.end();ait++){
			string a=*ait;
			j=0;
			for(set<string>::iterator bit=elements.begin();bit!=elements.end();bit++){
				string b=*bit;
				pair<string,string> p(a,b);
				if(inmatrix.find(p)!=inmatrix.end()){
					gsl_matrix_set(m,i,j,inmatrix[p]);
					pair<int,int> intedge(i,j);
					intedgeTostredge.insert( pair<pair<int,int>,pair<string,string> > (intedge,p) );
				}
				//cout<<"Filling gsl_matrix("<<i<<","<<j<<")= "<<gsl_matrix_get(m,i,j)<<endl;
				j++;
			}
			i++;
		}
		gsl_vector* S=gsl_vector_calloc(msize);
		gsl_matrix* V=gsl_matrix_calloc(msize,msize);
		//gsl_linalg_E_decomp_jacobi(m,V,S);
		gsl_eigen_symmv_workspace* w=gsl_eigen_symmv_alloc(msize);
		gsl_eigen_symmv(m,S,V,w);
		double sv_min;
		double sv_max;
		gsl_vector_minmax(S,&sv_min,&sv_max);
		if(!QUIET)cout<<"#Erange: min= "<<sv_min<<" max= "<<sv_max<<endl;
		switch(transf){
			case Exp:{
				double Tr_Exp_m=0.0;
				for(int i=0;i<msize;i++){
					Tr_Exp_m+=exp(_beta*gsl_vector_get(S,i));
				}
				FILE* f=fopen("/dev/stdout","w");
				cout<<"#BeginExponentiation; beta="<<_beta<<" Tr(e^(beta*m))= "<<Tr_Exp_m<<endl;
				if(SCIENTFMT)cout.setf(ios::scientific,ios::floatfield);
				if(SETPRECIS)cout.precision(precision);
				for(int i=0;i<msize;i++){
					for(int j=0;j<msize;j++){
						double expW_ij=0.0;
						for(int k=0;k<msize;k++){
							expW_ij+=gsl_matrix_get(V,i,k)*gsl_matrix_get(V,j,k)*exp(_beta*gsl_vector_get(S,k));
						}
						pair<string,string> p=intedgeTostredge[pair<int,int>(i,j)];
						cout<<p.first<<"\t"<<p.second<<"\t"<<expW_ij<<endl;
					}
				}
				cout<<"#EndExponentiation"<<endl;
				fclose(f);
				}
				break;
			case msval:{
				FILE* f=fopen("/dev/stdout","w");
				cout<<"#BeginEigenValues"<<endl;
				gsl_vector_fprintf(f,S,"%e");
				cout<<"#EndEigenValues"<<endl;
				fclose(f);
				}
				break;
			case mprtsvd:{
				FILE* f=fopen("/dev/stdout","w");
				fprintf(f,"#BeginEigenValueDecomposition:CEC^t\n");
				fprintf(f,"#BeginE:(row-major)\n");
				//gsl_matrix_fprintf(f,m,"%e");
				gsl_vector_fprintf(f,S,"%e");
				fprintf(f,"#EndE\n");
				fprintf(f,"#BeginC(row-major)\n");
				gsl_matrix_fprintf(f,V,"%e");
				fprintf(f,"#EndC\n");
				//gsl_vector_fprintf(f,S,"%g");
				fprintf(f,"#EndEigenValueDecomposition:CEC^t\n");
				fclose(f);
				}
				break;
			case mprunesval:
			case mprunesvalA:
			case mprunesvalB:{
				cout<<"#Ethreshold= "<<myarg<<endl;
				int nzE=msize;
				for(int k=0;k<msize;k++){
					//if(VERBOSE)cout<<"#svold= "<<gsl_vector_get(S,k)<<" (<"<<myarg<<"? ";
					if(VERBOSE)printf("#svold= %.16f (<%f? ",(float)gsl_vector_get(S,k),myarg);
					double e=gsl_vector_get(S,k);
					//if(gsl_vector_get(S,k)<myarg){
					//if((float)e<(float)myarg){
					bool compare=false;
					switch(transf){
						case mprunesvalA:
							compare=(float)e>(float)myarg;
							break;
						case mprunesval:
						case mprunesvalB:
							compare=(float)e<(float)myarg;
							break;
					}
					if(compare){
						gsl_vector_set(S,k,0.0);
						nzE--;
						if(VERBOSE)cout<<"TRUE)";
					}else
						if(VERBOSE)cout<<"FALSE)";
					if(VERBOSE)cout<<" => svnew= "<<gsl_vector_get(S,k)<<endl;;
				}
				if(VERBOSE)cout<<"#NumberOfNonzeroE="<<nzE<<endl;
				gsl_matrix* p=gsl_matrix_calloc(msize,msize);
				inmatrix.clear();
				int k=0;
				if(VERBOSE)cout<<"#Building filtered graph"<<endl;
				for(int i=0;i<msize;i++){
					if(PRT_E_GRAPH_TRIANG)k=i;
					for(int j=k;j<msize;j++){
						if(i==j&&PRT_E_GRAPH_NODIAG)continue;
						double v=0.0;
						//for(int k=0;k<msize;k++){
						for(int k=0;k<nzE;k++){
							if(VERBOSE)cout<<"#i="<<i<<",j="<<j<<" - K="<<k<<": v+="<<gsl_matrix_get(V,i,k)*gsl_vector_get(S,k)*gsl_matrix_get(V,j,k)<<" ("<<gsl_matrix_get(V,i,k)<<","<<gsl_vector_get(S,k)<<","<<gsl_matrix_get(V,j,k)<<")"<<endl;
							v+=gsl_matrix_get(V,i,k)*gsl_vector_get(S,k)*gsl_matrix_get(V,j,k);
						}
						//gsl_matrix_set(p,i,j,v);
						pair<int,int> intedge(i,j);
						if(intedgeTostredge.find(intedge)!=intedgeTostredge.end()){
							inmatrix[intedgeTostredge[intedge]]=v;	
						}
					}
				}
				set<pair<string,string> > edgedone;
				for(map<pair<string,string>,double>::iterator g=inmatrix.begin();g!=inmatrix.end();g++){
					string a=(g->first).first;
					string b=(g->first).second;
					double v=g->second;
					if(edgedone.find(g->first)==edgedone.end())
						cout<<a<<"\t"<<b<<"\t"<<v<<endl;
					edgedone.insert(g->first);
				}
				}
				break;
		}
		gsl_matrix_free(m);
		gsl_matrix_free(V);
		gsl_vector_free(S);
	}
	return 0;
}

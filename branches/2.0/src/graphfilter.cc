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
Copyright (C) Miguel A. Santos, HSC, Toronto, 2009.
Licensed under GPL version 3 a later. (see http://www.gnu.org/copyleft/gpl.html )
*/
//Compile: g++ -lgsl -lgslcblas -o graphfilter graphfilter.cc

//	VERSION 2.0

#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <set>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>

using namespace std;

extern	bool TRANSEDGE1;
extern	bool QUIET;
extern	bool VERBOSE;
extern	bool CENSOR;
extern	bool SCIENTFMT;
extern	bool SETPRECIS;
extern	bool SETDIAGVAL;
extern	bool SYMMETRIC;
extern	bool PRT_SV_GRAPH_FULL;
extern	bool PRT_SV_GRAPH_TRIANG;
extern	bool PRT_SV_GRAPH_NODIAG;

bool TRANSEDGE1;
bool QUIET=false;
bool VERBOSE=false;
bool CENSOR=false;
bool SCIENTFMT=false;
bool SETPRECIS=false;
bool SETDIAGVAL=false;
bool SYMMETRIC=true;
bool PRT_SV_GRAPH_FULL=true;
bool PRT_SV_GRAPH_TRIANG=!PRT_SV_GRAPH_FULL;
bool PRT_SV_GRAPH_NODIAG=!PRT_SV_GRAPH_FULL;

enum TRANS {mident, madd, mdiv, mmul, mlog, mlogneg, mexp, mpow , msvd , msval , mprunesval, mprtsvd};

void license(){
	cout<<\
"Graphfilter\n"<<\
"Copyright (C) Miguel A. Santos, HSC, Toronto, 2009.\n"<<\
"Licensed under GPL version 3 a later. (see http://www.gnu.org/copyleft/gpl.html )\n"<<\
	endl;
}
void usage(){
	cout<<"Usage: graphfilter [options] [graph.file]"<<endl;
	cout<<\
"\nGraphfilter allows you to perform different transformations on the edges\n\
of a graph, like taking the log/-log in any base, rasing to a given power,\n\
taking the exponential, etc, or filtering out some singular or eigenvalues.\n\
\n\
Graphs can be read from any table. The only condition is that the nodes\n\
are expected on the first two columns of the table. The location of the\n\
edge values can be specified with option -c.\n"<<\
"\nBy extension, it also works with any (non-square) matrix, given that the input \
abides that condition.\n\
Any matrix element A(i,j) is akin to an edge between nodes `i' and `j'.\n"<<\
	"\nUse option -h, --help for help.\n"<<\
	endl;
	cout<<"\nOptions:\n  General:\n"<<\
	"\t-c n \t\t\t\tSpecifies the edges to be at the n-th column\n"<<\
	"\t-h, --help \t\t\tPrints this help\n"<<\
	"\t-n, --precision n \t\tSet the number of decimals to use. Does not apply for SVD.\n"<<\
	"\t-G, --print-no-diagonal \tWhen dealing with SV, do no print out the diagonal\n"<<\
	"\t-V, --print-svd \t\tPrint matrices U and V from the Singular Value Decomposition of the given graph G: G=USV^t\n"<<\
	"\t-q \t\t\t\tQuiet mode: do not print any comment lines.\n"<<\
	"\t-Q \t\t\t\tCensor mode: Be quiet and also skip any comment lines from input file. Comment lines are those starting with `#'\n"<<\
	"\t-g, --scientific \t\tUse scientific notation for the output. Does not apply for SVD.\n"<<\
	"\t--verbose \t\t\tPrints out more detailed info on what's going on.\n"<<\
	"\t-X, --non-symmetric \t\tDo _not_ assume provided matrix is symmetric. By default we assume it is.\n"<<\
	"  Transformations:\n"<<\
	"\t-A, --add w \t\t\tAdds w to all edges\n"<<\
	"\t-B, --base b  \t\t\tSpecifies the base of logarithms to use\n"<<\
	"\t-C, --diag w \t\t\t(Re)Sets all diagonal values to w.\n"<<\
	"\t-D, --div w \t\t\tDivides all edges by w\n"<<\
	"\t-E, --exp  \t\t\tTakes the exp of the edges\n"<<\
	"\t-I, --identity  \t\tDo nothing but printing out the edges. Default action.\n"<<\
	"\t-L, --log  \t\t\tTakes the log (base e) of the edges\n"<<\
	"\t-N, --log-neg  \t\t\tTakes the negative log (base e) of the edges\n"<<\
	"\t-M, --mul w \t\t\tMultiplies all edges by w\n"<<\
	"\t-P, --pow n \t\t\tTakes the n-power of the edges\n"<<\
	"\t-T, --print-triangular \t\tWhen dealing with SV, just print out the upper/lower triangular part\n"<<\
	"\t-p, --prune-singular-values sv \tRemoves all singular values below sv\n"<<\
	"\t-v, --singular-values \t\tGet the singular values of the matrix\n"<<\
	"\t-S, --symmetric \t\tAssume provided matrix is symmetric (Default).\n"<<\
	endl;
	cout<<\
"Examples:\n\
\tcat graph | graphfilter -c 11 --log-neg --base 10 > graph-log10\n\
\n\tConsiders column 11 of file graph as the edge values and takes the -log of them in base 10.\n\
\n\tcat graph | graphfilter -c 11 --mult 3.5 | graphfilter --exp > graph-exp3.5\n\
\n\tThe file graph-exp3.5 will contain the exponential of original edges to the power of 3.5, i.e, exp(3.5*edge)\n\
\tBy default graphfilter expects the edges in column 3 (The first column is 1).\n"<<\
"\n\tFilter out all singular values below 100 and get the new matrix values:      \n\
\tgraphfilter graph2 --prune-singular-values 100\n"<<\
	endl;
	license();
	exit(1);
}

class SVD {
	gsl_matrix* _U;
	gsl_matrix* _V;
	gsl_vector* _S;
	double _sv_min;
	double _sv_max;
	size_t _msize;
	size_t _nsize;
	size_t _m;
	size_t _n;
public:
	SVD(gsl_matrix* m, size_t m, size_t n);
	SVD(gsl_matrix* m, size_t m);
	~SVD();
	void print_sval(const char* fmt="%e");
	void print_svec(const char* fmt="%e");
	gsl_matrix* prune_eval(double evthr, map< pair<int,int> , pair<string,string> >& intedgeTostredge);
};

SVD::~SVD(){
	gsl_matrix_free(_U);
	gsl_matrix_free(_V);
	gsl_vector_free(_S);
}

SVD::SVD(gsl_matrix* A, size_t m, size_t n){
	_m=m;
	_n=n;
	//_U=gsl_matrix_calloc(m,n);
	_U=A;
	_V=gsl_matrix_calloc(m,n);
	_msize=m>n?m:n;
	_nsize=m<n?m:n;
	_S=gsl_vector_calloc(_msize);
	gsl_linalg_SV_decomp_jacobi(_U,_V,_S);
}

SVD::SVD(gsl_matrix* A, size_t n){
	_m=n;
	_n=n;
	//_U=gsl_matrix_calloc(n,n);
	_U=A;
	_V=gsl_matrix_calloc(n,n);
	_msize=n;
	_nsize=n;
	_S=gsl_vector_calloc(n);
	gsl_linalg_SV_decomp_jacobi(_U,_V,_S);
}

gsl_matrix* SVD::prune_eval(double evthr, map< pair<int,int> , pair<string,string> >& intedgeTostredge){
	cout<<"#SVthreshold= "<<evthr<<endl;
	int nzSV=_msize;
	for(int k=0;k<_msize;k++){
		if(VERBOSE)cout<<"#svold= "<<gsl_vector_get(_S,k)<<" => svnew= ";
		if(gsl_vector_get(_S,k)<evthr){
			gsl_vector_set(_S,k,0.0);
			nzSV--;
		}
		if(VERBOSE)cout<<gsl_vector_get(_S,k)<<endl;;
	}
	if(VERBOSE)cout<<"#NumberOfNonzeroSV="<<nzSV<<endl;
	gsl_matrix* p=gsl_matrix_calloc(_m,_n);
	map< pair<string,string >, double > inmatrix;
	int k=0;
	if(VERBOSE)cout<<"#Building filtered graph"<<endl;
	for(int i=0;i<_m;i++){
		if(PRT_SV_GRAPH_TRIANG)k=i;
		for(int j=k;j<_n;j++){
			if(i==j&&PRT_SV_GRAPH_NODIAG)continue;
			double v=0.0;
			for(int k=0;k<nzSV;k++){
				if(VERBOSE)cout<<"#i="<<i<<",j="<<j<<" - K="<<k<<": v+="<<gsl_matrix_get(_U,i,k)*gsl_vector_get(_S,k)*gsl_matrix_get(_V,j,k)<<endl;
				v+=gsl_matrix_get(_U,i,k)*gsl_vector_get(_S,k)*gsl_matrix_get(_V,j,k);
			}
			gsl_matrix_set(p,i,j,v);
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
	return p;
}

void SVD::print_sval(const char* fmt){
	FILE* f=fopen("/dev/stdout","w");
	cout<<"#BeginSingularValues"<<endl;
	gsl_vector_fprintf(f,_S,fmt);
	cout<<"#EndSingularValues"<<endl;
	fclose(f);
}

void SVD::print_svec(const char* fmt){
	FILE* f=fopen("/dev/stdout","w");
	fprintf(f,"#BeginSingularValueDecomposition:USV^t\n");
	fprintf(f,"#BeginU:(row-major)\n");
	gsl_matrix_fprintf(f,_U,fmt);
	fprintf(f,"#EndU\n");
	fprintf(f,"#BeginV(row-major)\n");
	gsl_matrix_fprintf(f,_V,fmt);
	fprintf(f,"#EndV\n");
	//gsl_vector_fprintf(f,S,"%g");
	fprintf(f,"#EndSingularValueDecomposition:USV^t\n");
	fclose(f);
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
	//bool TRANSEDGE1=true;
	//bool QUIET=false;
	//bool VERBOSE=false;
	//bool CENSOR=false;
	//bool SCIENTFMT=false;
	//bool SETPRECIS=false;
	//bool SETDIAGVAL=false;
	//bool SYMMETRIC=true;
	//bool PRT_SV_GRAPH_FULL=true;
	//bool PRT_SV_GRAPH_TRIANG=!PRT_SV_GRAPH_FULL;
	//bool PRT_SV_GRAPH_NODIAG=!PRT_SV_GRAPH_FULL;
	char* fn="/dev/stdin";
	double myarg;
	double mdiagval;
	double logbase=1.0;
	int precision=10;
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
		else if(strcmp(argv[0],"--scientific")==0||strcmp(argv[0],"-g")==0){
			argc--,argv++;
			SCIENTFMT=true;
		}
		else if(strcmp(argv[0],"--precision")==0||strcmp(argv[0],"-n")==0){
			argc--,argv++;
			if(argc<1)usage();
			precision=atoi(*argv);
			argc--,argv++;
			SETPRECIS=true;
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
		else if(strcmp(argv[0],"--base")==0||strcmp(argv[0],"-B")==0){
			argc--,argv++;
			if(argc<1)usage();
			logbase=log(atof(*argv));
			stringstream ss;
			ss<<"base="<<*argv;
			ss>>addcomment;
			argc--,argv++;
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
		else if(strcmp(argv[0],"--print-no-diagonal")==0||strcmp(argv[0],"-G")==0){
			PRT_SV_GRAPH_FULL=false;
			PRT_SV_GRAPH_NODIAG=true;
			argc--,argv++;
		}
		else if(strcmp(argv[0],"--print-svd")==0||strcmp(argv[0],"-V")==0){
			transf=mprtsvd;
			TRANSEDGE1=false;
			transftype="print-SVD";
			argc--,argv++;
		}
		else if(strcmp(argv[0],"--prune-singular-values")==0||strcmp(argv[0],"-p")==0){
			transf=mprunesval;
			TRANSEDGE1=false;
			transftype="prune-singular-values-below:\t";
			argc--,argv++;
			if(argc<1)usage();
			myarg=atof(argv[0]);
			argc--,argv++;
			stringstream ss;
			ss<<transftype<<myarg;
			transftype=ss.str();
		}
		else if(strcmp(argv[0],"--print-triangular")==0||strcmp(argv[0],"-T")==0){
			PRT_SV_GRAPH_FULL=false;
			PRT_SV_GRAPH_TRIANG=true;
			argc--,argv++;
		}
		else if(strcmp(argv[0],"--singular-values")==0||strcmp(argv[0],"-v")==0){
			transf=msval;
			TRANSEDGE1=false;
			transftype="singular-values";
			argc--,argv++;
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
		cout<<"#Type: "<<transftype<<" "<<addcomment<<"\n#"<<endl;
		if(SETDIAGVAL)cout<<"#Setting diag values to: "<<mdiagval<<endl;
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
		if(SETDIAGVAL){
			if(a.compare(b)==0)
				f=mdiagval;
			else{
				inmatrix.insert(pair< pair<string,string>, double> (pair<string,string>(a,a),mdiagval));
				inmatrix.insert(pair< pair<string,string>, double> (pair<string,string>(b,b),mdiagval));
			}
		}
		inmatrix.insert(pair< pair<string,string>, double> (pair<string,string>(a,b),f));
		if(SYMMETRIC&&a.compare(b)!=0)inmatrix.insert(pair< pair<string,string>, double> (pair<string,string>(b,a),f));
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
			case mident:
				break;
		}
		if(TRANSEDGE1){
			cout<<a<<"\t"<<b<<"\t";
			if(SCIENTFMT)cout.setf(ios::scientific,ios::floatfield);
			if(SETPRECIS)cout.precision(precision);
			cout<<f<<endl;
			//cout<<a<<"\t"<<b<<"\t"<<f<<endl;
		}
	}	
	if(!QUIET){
		size_t efNedges=inmatrix.size();
		if(SETDIAGVAL) efNedges-=elements.size();
		if(SYMMETRIC) efNedges/=2;
		cout<<"#elements= "<<elements.size()<<" #edges= "<<efNedges<<endl;
	}
	if(!TRANSEDGE1){
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
		gsl_eigen_symm_workspace* w=gsl_eigen_symm_alloc(msize);
		//gsl_linalg_SV_decomp_jacobi(m,V,S);
		gsl_eigen_symm(m,S,w);
		double sv_min;
		double sv_max;
		gsl_vector_minmax(S,&sv_min,&sv_max);
		if(!QUIET)cout<<"#SVrange: min= "<<sv_min<<" max= "<<sv_max<<endl;
		switch(transf){
			case msval:{
				FILE* f=fopen("/dev/stdout","w");
				cout<<"#BeginSingularValues"<<endl;
				gsl_vector_fprintf(f,S,"%e");
				cout<<"#EndSingularValues"<<endl;
				fclose(f);
				}
				break;
			case mprtsvd:{
				FILE* f=fopen("/dev/stdout","w");
				fprintf(f,"#BeginSingularValueDecomposition:USV^t\n");
				fprintf(f,"#BeginU:(row-major)\n");
				gsl_matrix_fprintf(f,m,"%e");
				fprintf(f,"#EndU\n");
				fprintf(f,"#BeginV(row-major)\n");
				gsl_matrix_fprintf(f,V,"%e");
				fprintf(f,"#EndV\n");
				//gsl_vector_fprintf(f,S,"%g");
				fprintf(f,"#EndSingularValueDecomposition:USV^t\n");
				fclose(f);
				}
				break;
			case mprunesval:{
				cout<<"#SVthreshold= "<<myarg<<endl;
				int nzSV=msize;
				for(int k=0;k<msize;k++){
					if(VERBOSE)cout<<"#svold= "<<gsl_vector_get(S,k)<<" => svnew= ";
					if(gsl_vector_get(S,k)<myarg){
						gsl_vector_set(S,k,0.0);
						nzSV--;
					}
					if(VERBOSE)cout<<gsl_vector_get(S,k)<<endl;;
				}
				if(VERBOSE)cout<<"#NumberOfNonzeroSV="<<nzSV<<endl;
				gsl_matrix* p=gsl_matrix_calloc(msize,msize);
				inmatrix.clear();
				int k=0;
				if(VERBOSE)cout<<"#Building filtered graph"<<endl;
				for(int i=0;i<msize;i++){
					if(PRT_SV_GRAPH_TRIANG)k=i;
					for(int j=k;j<msize;j++){
						if(i==j&&PRT_SV_GRAPH_NODIAG)continue;
						double v=0.0;
						//for(int k=0;k<msize;k++){
						for(int k=0;k<nzSV;k++){
							if(VERBOSE)cout<<"#i="<<i<<",j="<<j<<" - K="<<k<<": v+="<<gsl_matrix_get(m,i,k)*gsl_vector_get(S,k)*gsl_matrix_get(V,j,k)<<endl;
							v+=gsl_matrix_get(m,i,k)*gsl_vector_get(S,k)*gsl_matrix_get(V,j,k);
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

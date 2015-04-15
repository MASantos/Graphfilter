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


#ifndef _CLASS_SVD_
#define _CLASS_SVD_ 1

#include "SVD.h" 

using namespace std;

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

#endif // end _CLASS_SVD_

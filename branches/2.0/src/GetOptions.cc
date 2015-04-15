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



#ifndef _CLASS_GETOPTIONS_
#define _CLASS_GETOPTIONS_ 1

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
extern  enum TRANS;
extern  void usage();

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


class GetOptions{
	map<string,string> _optargs;
	
public:
	GetOptions(int& argc, char** argv);
}

GetOptions::GetOptions(int& argc, char** argv){
	_optargs.clear();
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
			string opt=*argv;
			argc--,argv++;
			if(argc<1)usage();
			col=atoi(argv[0]);
			string arg=*argv;
			argc--,argv++;
			_optargs[opt]=arg;
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
}

#endif //_CLAS_GETOPTIONS

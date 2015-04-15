/*
*    This is programm transedges.
*
*    Transedges is free software: you can redistribute it and/or modify
*    it under the terms of the GNU General Public License as published by
*    the Free Software Foundation, either version 3 of the License, or
*    (at your option) any later version.
*
*    Transedges is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*    GNU General Public License for more details.
*
*    You should have received a copy of the GNU General Public License
*    along with Transedges.  If not, see <http://www.gnu.org/licenses/>.
*
*/
/** Transedges
Copyright (C) Miguel A. Santos, HSC, Toronto, 2009.
Licensed under GPL version 3 a later. (see http://www.gnu.org/copyleft/gpl.html )
*/
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

enum map {mident, mlog, mlogneg, mexp, mpow };

void license(){
	cout<<\
"Transedges\n"<<\
"Copyright (C) Miguel A. Santos, HSC, Toronto, 2009.\n"<<\
"Licensed under GPL version 3 a later. (see http://www.gnu.org/copyleft/gpl.html )\n"<<\
	endl;
}
void usage(){
	cout<<"Usage: transedges [options] [graph.file]"<<endl;
	cout<<"\n\tOptions:\n"<<\
	"\t-c n Specifies the edges to be at the n-th column\n"<<\
	"\t-q Quiet mode: do not print any comment lines.\n"<<\
	"\t-Q Censor mode: Be quiet and also skip any comment lines from input file. Comment lines are those starting with `#'\n"<<\
	"\t--base b  Specifies the base of logarithms to use\n"<<\
	"\t--exp  Takes the exp of the edges\n"<<\
	"\t--log  Takes the log (base e) of the edges\n"<<\
	"\t--log-neg  Takes the negative log (base e) of the edges\n"<<\
	"\t--pow n Takes the n-power of the edges\n"<<\
	"\t--precision n Set the number of decimals to use\n"<<\
	"\t--scientific Use scientific notation for the output\n"<<\
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
	map transf=mident;
	string transftype="Identity";
	string addcomment="";
	bool QUIET=false;
	bool CENSOR=false;
	bool SCIENTFMT=false;
	bool SETPRECIS=false;
	char* fn="/dev/stdin";
	double myarg;
	double logbase=1.0;
	int precision=10;
	while(argc>0){
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
		else if(strcmp(*argv,"-c")==0){
			argc--,argv++;
			if(argc<1)usage();
			col=atoi(argv[0]);
			argc--,argv++;
		}
		else if(strcmp(argv[0],"--log")==0){
			transf=mlog;
			transftype="ln";
			argc--,argv++;
		}
		else if(strcmp(argv[0],"--log-neg")==0){
			transf=mlogneg;
			transftype="-ln";
			argc--,argv++;
		}
		else if(strcmp(argv[0],"--base")==0){
			argc--,argv++;
			if(argc<1)usage();
			logbase=log(atof(*argv));
			stringstream ss;
			ss<<"base="<<*argv;
			ss>>addcomment;
			argc--,argv++;
		}
		else if(strcmp(argv[0],"--scientific")==0){
			argc--,argv++;
			SCIENTFMT=true;
		}
		else if(strcmp(argv[0],"--precision")==0){
			argc--,argv++;
			if(argc<1)usage();
			precision=atoi(*argv);
			argc--,argv++;
			SETPRECIS=true;
		}
		else if(strcmp(argv[0],"--exp")==0){
			transf=mexp;
			transftype="exp";
			argc--,argv++;
		}
		else if(strcmp(argv[0],"--pow")==0){
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
		else{
			fn=argv[0];
			argc--,argv++;
		}
	}
	if(!QUIET){
		cout<<"#Transforming edges of graph "<<fn<<endl;
		cout<<"#Type: "<<transftype<<" "<<addcomment<<"\n#"<<endl;
	}
	argc--,argv++;
	ifstream is(fn);
	if(!is){
		cout<<"Error: Could not open file "<<fn<<endl;
		exit(1);
	}
	string ln;
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
		int shft=col-2;
		while(shft-->1)ss>>c;
		ss>>f;
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
			case mident:
				break;
		}
		cout<<a<<"\t"<<b<<"\t";
		if(SCIENTFMT)cout.setf(ios::scientific,ios::floatfield);
		if(SETPRECIS)cout.precision(precision);
		cout<<f<<endl;
		//cout<<a<<"\t"<<b<<"\t"<<f<<endl;
	}
	return 0;
}

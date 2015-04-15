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



#ifndef _CLASS_GETOPTIONS_H
#define _CLASS_GETOPTIONS_H 1

#include <iostream>
#include <sstream>
#include <map>

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


class GetOptions{
	map<string,string> _optargs;
	
public:
	GetOptions(int& argc, char** argv);
};
#endif //_CLAS_GETOPTIONS_H

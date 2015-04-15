#include<iostream>
#include<limits>
#include<string> 

using namespace std;

int main(int argc, char** argv){
	argc--;argv++;
	while(argc>0){
		double v=atof(*argv);	
		if(isnan(v))cout<<"NAN!"<<endl;
		argc--;argv++;
	}
	string a="abc\n";
	cout<<a;
	printf("printf=>%s",a.c_str());
	size_t i=a.rfind("\n");
	a.erase(a.rfind("\n"));
	cout<<"After erasing: a="<<a<<":Newline??"<<endl;
	const char* b=a.c_str();
	cout<<"b="<<b;
	cout<<"Finished"<<endl;
	return 0;
}

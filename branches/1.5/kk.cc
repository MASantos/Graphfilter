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
	int s=2;
	double v=1/3.0;
	cout<<s<<"<"<<v<<"? ";
	cout<<"using printf"<<endl;
	printf("(20digits): %d < v=%.20f?\n",s,v);
	if(s<v)cout<<"TRUE"<<endl;
	else cout<<"FALSE"<<endl;
	printf("(17 digits): 1/3-v=%.17f\n",0.33333333333333333-v);
	printf("(20 digits): 1/3-v=%.20f\n",0.33333333333333333333-v);
	printf("(20 digits): 1/3/v=%.20f\n",0.33333333333333333333/v);
	return 0;
}

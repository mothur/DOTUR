//
//

using namespace std;

#include "otu.h"
#include "richness.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
#include <list>
#include <ctime>

/**************************************************************************************************/

void usageError(char *name)
{
	cerr << "Usage: " << name << " [-i Iterations (<1000>)] [-c ClusterMethod (<f>, n, a)] [-p Precision (10, <100>, 1000, 10000)] [-l] [-j] <file>" << endl;
	cerr << "Options:\n";
	cerr << " -i:		Number of iterations (default = 1000)\n";
	cerr << " -c:		Clustering method - (f) furthest neighbor, (n) nearest neighbor,\n";
	cerr << "			(a) average neighbor (default = f)\n";
	cerr << " -p:		Precision of distances for output, increasing can dramatically\n";
	cerr << "			lengthen execution times - 10, 100, 1000, 10000 (default = 100)\n";
	cerr << " -l:		Input file is lower triangular (default = square matrix)\n";
	cerr << " -r:		Calculates rarefaction curves for each parameter, can dramatically\n";
	cerr << "			lengthen execution times.  Simple rarefaction curve always calculated.\n";
	cerr << " -stop:	Stops clustering when cutoff has been reached.\n";
	cerr << " -wrep:	Samples with replacement.\n";
	cerr << " -jumble:	Jumble the order of the distance matrix.\n";
	cerr << " -sim:		Converts similarity score to distance (D=1-S).\n";
	exit(1);
}

/**************************************************************************************************/

void init_files(string fname, int rarefy, int sample, vector<string>& r_names, vector<string>& c_names){
	
	r_names.push_back(fname + ".rarefaction");
	if(sample == 1){
		r_names.push_back(fname + ".wr_rarefaction");
	}
	
	c_names.push_back(fname + ".collect");
	c_names.push_back(fname + ".c_chao");
	c_names.push_back(fname + ".c_ace");
	c_names.push_back(fname + ".c_jack");
	c_names.push_back(fname + ".c_boot");
	c_names.push_back(fname + ".c_simpson");
	c_names.push_back(fname + ".c_shannon");

	for(int i=0;i<c_names.size();i++){
		ofstream cfile(c_names[i].c_str(), ios::trunc);
	}

	if(rarefy){
		r_names.push_back(fname + ".r_chao");
		r_names.push_back(fname + ".r_ace");
		r_names.push_back(fname + ".r_jack");
		r_names.push_back(fname + ".r_boot");
		r_names.push_back(fname + ".r_simpson");
		r_names.push_back(fname + ".r_shannon");
		
	}

	for(int i=0;i<r_names.size();i++){
		ofstream rfile(r_names[i].c_str(), ios::trunc);
	}
	
}

/**************************************************************************************************/

int main(int argc, char *argv[])
{
//	srand(54321);
	srand( (unsigned)time( NULL ) );
	cout.setf(ios::fixed, ios::floatfield);
	cout.setf(ios::showpoint);
	cerr.setf(ios::fixed, ios::floatfield);
	cerr.setf(ios::showpoint);

	DistanceMatrix dist;
	OTUMatrix otu;

	string filename = "";
	char cmethod = 'f';
	otu.rarefy = 0;
	dist.cutoff = 1000;
	int matrix = 1;
	int precise = 100;
	otu.iter = 1000;
	otu.sample = 0;
	dist.jumble = 0;
	dist.similarity = 0;
	string otuname, listname, indexname, fbase;
	
	char **p;
	if(argc>1){	
		for(p=argv+1;p<argv+argc;p++){
			if(strcmp(*p, "-c")==0){
				if(++p>=argv+argc)
					usageError(argv[0]);
				istringstream f(*p);
				if(!(f >> cmethod))
					usageError(argv[0]);
			}
			else if(strcmp(*p, "-l")==0){
				if(p>=argv+argc)
					usageError(argv[0]);
				matrix=2; //lower triangle
			}
			else if(strcmp(*p, "-r")==0){
				if(p>=argv+argc)
					usageError(argv[0]);
				otu.rarefy=1; //calculate rarefaction for parameters
			}
			else if(strcmp(*p, "-stop")==0){
				if(++p>=argv+argc)
					usageError(argv[0]);
				istringstream f(*p);
				if(!(f >> dist.cutoff))
					usageError(argv[0]);
			}
			else if(strcmp(*p, "-wrep")==0){
				if(p>=argv+argc)
					usageError(argv[0]);
				otu.sample=1; //sample with replacement
			}
			else if(strcmp(*p, "-jumble") == 0){
				if(p>=argv+argc)
					usageError(argv[0]);
				dist.jumble=1; 
			}
			else if(strcmp(*p, "-sim") == 0){
				if(p>=argv+argc)
					usageError(argv[0]);
				dist.similarity=1; 
			}
			else if(strcmp(*p, "-p")==0){
				if(++p>=argv+argc)
					usageError(argv[0]);
				istringstream f(*p);
				if(!(f >> precise))
					usageError(argv[0]);
			}
			else if(strcmp(*p, "-i")==0){
				if(++p>=argv+argc)
					usageError(argv[0]);
				istringstream f(*p);
				if(!(f >> otu.iter))
					usageError(argv[0]);
			}
		  	else {
				istringstream f(*p);
				if(!(f >> filename))
					usageError(argv[0]);
			}
		}
		fbase = filename;
		if(fbase.find_last_of(".")!=string::npos){
			fbase.erase(fbase.find_last_of(".")+1); 
		}
		else{
			fbase += ".";
		}
		dist.read(filename, matrix);

 		otuname =  fbase + cmethod + "n.otu";
		listname = fbase + cmethod + "n.list";
		indexname = fbase + cmethod + "n.index";
		
		ofstream otu_file(otuname.c_str(), ios::trunc);
		ofstream list_file(listname.c_str(), ios::trunc);
		ofstream index_file(indexname.c_str(), ios::trunc);
	}
	else
		usageError(argv[0]);

  	cout << endl;
	cout << "|-------------------------------------------------------|" << endl;
	cout << "|             |               % Complete                |" << endl;
	cout << "| Progress    |  10  20  30  40  50  60  70  80  90 100 |" << endl;
  	cout << "|-------------------------------------------------------|" << endl;
	cout << "| OTU Ident   |";
	cout.flush();
	
	dist.initw();
	vector<float> dnames(1);

	otu.distances = 0;
	int incr = dist.rank/10;
  	int check = incr;
	int z = 0;

	while(dist.rank > 1){
		dist.small_element();
		dist.cluster(cmethod);
		dist.merge();
		if(dist.print_nonred_w(otuname, listname, indexname, precise, otu.distances, dnames) == 1){	break;		}
   		if(z==check) {
			cout << "****";
			cout.flush();
			check += incr;
		}
		z++;
	}
	cout << " |\n";
	
	string fname = fbase+cmethod+ "n";
	vector<string> r_names;
	vector<string> c_names;
	otu.rank_abundance(fname);

	init_files(fname,otu.rarefy, otu.sample, r_names, c_names);

	otu.binomial();
	
	otu.file_root = fname;
	cout << "| Collectors  |";
	cout.flush();
	otu.readlist(dist.orig_name_list);
	
	cout << "| Rarefaction |";
	cout.flush();	
	otu.readotu();
	
	cout << "|-------------------------------------------------------|" << endl << endl;

	cout << "Cleaning up..." << endl;
	cout.flush();

	for(int i=0;i<r_names.size();i++){		otu.output("no",  r_names[i]);		}
	for(int i=0;i<c_names.size();i++){		otu.output("yes", c_names[i]);		}

	return(0);
}

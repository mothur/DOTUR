#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
#include <list>
#include <cmath>

class DistanceMatrix {
public:
	DistanceMatrix(){};
	void read(string, int);
	void initw();
	void print_matrix();
	int print_nonred_w(string, string, string, int,  int&, vector<float>&);
	void print_w();
	void small_element();
	void cluster(char);
	void merge();
	int jumble;
	int similarity;
	int rank;
	float cutoff;
	vector<string> name_list;
	vector<string> orig_name_list;
private:
	void read_phylip(istream&, int);
	void read_mega(istream&, int);
	float distance;
	float prevdist;
	float min_compare(float, float);
	int n;
	int test;
	int coord[2];
	vector<vector<float> > d;
	vector<int> w;
	vector<int> prevw;
	vector<string> prev_name_list;
	vector<string> indices;
	vector<string> prev_indices;
	
};

/**************************************************************************************************/

template<typename T>
string tostring(const T&x){
	ostringstream out;
	out << x;
	return out.str();
}

/**************************************************************************************************/

void get_comment(istream& f, char begin, char end)
{
	char d=f.get();
	while(d != end){	d = f.get();	}
	d = f.peek();
}	

/**************************************************************************************************/

void DistanceMatrix::read(string fileName, int square_m)
{
	ifstream f(fileName.c_str());
	if(!f) {
		cerr << "Error: Could not open " << fileName << endl;
		exit(1);
	}
	char test = f.peek();
	
	if(test == '#'){
		read_mega(f, square_m);
	}
	else{
		read_phylip(f, square_m);
	}
	
	if(jumble == 1){
		vector<int> order(rank);
		
		for(int i=0;i<rank;i++){			order[i]=i;			}
	
		for(int j=rank-1;j>=0;j--){
			int z = int((float)(j+1) * (float)(rand()) / ((float)RAND_MAX+1.0));	
			int t = order[z]; order[z]=order[j];order[j]=t;
		}
		
		vector<vector<float> > jumble_d(rank);
		vector<string> jumble_name_list(rank);
		
		for(int i=0;i<rank;i++){
			jumble_d[i].resize(rank);
			jumble_name_list[i] = name_list[order[i]];
			for(int j=0;j<rank;j++){
				jumble_d[i][j] = d[order[i]][order[j]];
			}
		}	
		
		name_list = jumble_name_list;
		d = jumble_d;

	}

	if(similarity == 1){
		for(int i=0;i<rank;i++){
			for(int j=0;j<rank;j++){
				d[i][j] = 1.0 - d[i][j];
			}
		}
	}
	orig_name_list = name_list;

}

/**************************************************************************************************/

void DistanceMatrix::read_mega(istream& f, int square_m)
{
	get_comment(f, '#', '\n');

	char test = f.peek();

	while(test == '!'){								//get header comments
		get_comment(f, '!', ';');
		while(isspace(test=f.get()))		{;}
		f.putback(test);
		test = f.peek();
	}
	while(test != '\n'){							//get sequence names
		get_comment(f, '[', ']');
		char d = f.get();
		d = f.get();
		if(d == '#'){
			string name;
			f >> name;
			name_list.push_back(name);
			while(isspace(test=f.get()))		{;}
			f.putback(test);
		}
		else{
			break;
		}
	}
	rank = name_list.size();
	d.resize(rank);
	for(int i=0;i<rank;i++){		d[i].resize(rank);		}

	d[0][0] = 0.0000;
	get_comment(f, '[', ']');
	for(int i=1;i<rank;i++){
		get_comment(f, '[', ']');
		d[i][i]=0.0000;
		for(int j=0;j<i;j++){
			f >> d[i][j];
			if (d[i][j] == -0.0000)
				d[i][j] = 0.0000;
			d[j][i]=d[i][j];
		}
	}
}

/**************************************************************************************************/

void DistanceMatrix::read_phylip(istream& f, int square_m)
{
	int count1=0;
	int count2=0;

	f >> rank;

	name_list.resize(rank);
	d.resize(rank);
	if(square_m == 1){
		for(int i=0;i<rank;i++)
			d[i].resize(rank);
		for(int i=0;i<rank;i++) {
			f >> name_list[i];
			for(int j=0;j<rank;j++) {
				f >> d[i][j];
				if (d[i][j] == -0.0000)
					d[i][j] = 0.0000;
			}
		}
	}
	else if(square_m == 2){
		for(int i=0;i<rank;i++){
			d[i].resize(rank);
		}
		d[0][0] = 0.0000;
		f >> name_list[0];
		for(int i=1;i<rank;i++){
			f >> name_list[i];
			d[i][i]=0.0000;
			for(int j=0;j<i;j++){
				f >> d[i][j];
				if (d[i][j] == -0.0000)
					d[i][j] = 0.0000;
				d[j][i]=d[i][j];
			}
		}
	}
}
	
/**************************************************************************************************/

void DistanceMatrix::initw(void)
{
	w.resize(rank);
	prevw.resize(rank);
	indices.resize(rank);
	prev_name_list.resize(rank);
	prev_indices.resize(rank);

	for(int i=0;i<rank;i++){
		w[i] = 1;
		indices[i] = tostring(i);
	}
	prevdist = 0.0000;
	prev_indices = indices;
	prev_name_list = name_list;
	prevw = w;
}

/**************************************************************************************************/

void DistanceMatrix::small_element()
{
	distance = 100.0;
	coord[0] = 0;
	coord[1] = 0;
	
	for(int i=0;i<rank-1;i++){
		for(int j=i+1;j<rank;j++){
			if (d[i][j] < distance){
				distance = d[i][j];
				coord[0] = i;
				coord[1] = j;
			}
		}
	}
}

/**************************************************************************************************/

float DistanceMatrix::min_compare(float x, float y)
{
    if(x>y){
        return (y);
    }
    else return (x);
}

/**************************************************************************************************/

void DistanceMatrix::cluster(char method)
{
	int r=coord[0];
	int c=coord[1];
	
    for(int i=0;i<rank;i++){
        d[i][i]=0.0000;
        if(i != r && i != c){
            if(method == 'f'){
                d[r][i] = -min_compare(-d[r][i],-d[c][i]);
            }
            else if (method == 'a'){
                d[r][i] = (d[r][i]*w[r]+d[c][i]*w[c])/(w[r]+w[c]);
            }
            else if (method == 'n'){
                d[r][i] = min_compare(d[r][i], d[c][i]);
            }
            d[i][r]=d[r][i];
        }
    }
}

/****************************************************************************/

void DistanceMatrix::merge()
{
	int r = coord[0];
	int c = coord[1];
    
	w[r]=w[r]+w[c];
	name_list[r] = name_list[r] + "," + name_list[c];
	indices[r]=indices[r] + "," + indices[c];
	
    for(int i=c+1;i<rank;i++){
		for(int j=0;j < rank;j++){
			d[i-1][j]=d[i][j];
		}
		w[i-1]=w[i];
		name_list[i-1]=name_list[i];
		indices[i-1]=indices[i];
    }
	for(int i=0;i<rank;i++){
		for(int j=c+1;j<rank;j++){
			d[i][j-1]=d[i][j];
		}
	}

	rank--;

	d.resize(rank);
	w.resize(rank);
	name_list.resize(rank);
	indices.resize(rank);
}

/**************************************************************************************************/

int DistanceMatrix::print_nonred_w(string otu, string list, string index, int precise, int& count, vector<float>& dnames)
{
	ofstream otu_file(otu.c_str(), ios::app|ios::ate);
	otu_file.setf(ios::fixed, ios::floatfield);
	otu_file.setf(ios::showpoint);

	ofstream list_file(list.c_str(), ios::app|ios::ate);
	list_file.setf(ios::fixed, ios::floatfield);
	list_file.setf(ios::showpoint);

	ofstream index_file(index.c_str(), ios::app|ios::ate);
	index_file.setf(ios::fixed, ios::floatfield);
	index_file.setf(ios::showpoint);

	int result;
	
	if (prevdist == 0.0000 && prevdist != distance && precise != 10000){
		count++;
		dnames.resize(count);
		dnames[count-1] = -1.0000;
				
		sort(prevw.begin(), prevw.end());
		otu_file << "unique\t" << rank+1;
		list_file << "unique\t"<< rank+1;
		index_file << "unique\t"<< rank+1;

		for(int i=0;i<rank+1;i++){
			if(prevw[i] > 0){
				otu_file << "\t" << prevw[i];
			}else{
				otu_file << "\t" << -prevw[i];
			}
			list_file << "\t" << prev_name_list[i];
			index_file << "\t" << prev_indices[i];
		}
		otu_file << "\n";
		list_file << "\n";
		index_file << ";\n";
	}
	int rprevdist = (int)(precise*prevdist+0.5);
	int rdistance = (int)(precise*distance+0.5);

	int precision;
	switch(precise){
		case 10: precision = 1; break;
		case 100: precision = 2; break;
		case 1000: precision = 3; break;
		case 10000: precision = 4; break;
		default: precision = 2; break;
	}
	
	if(rdistance != rprevdist){
		if(rprevdist > cutoff * precise){
			return 1;
		}
		
		
		count++;
		dnames.resize(count);
		dnames[count-1]=(float)rprevdist/(float)precise;
		
		otu_file  <<	setprecision(precision)	<<	(float)rprevdist/(float)precise	<<	"\t"	<<	rank+1;
		list_file <<	setprecision(precision)	<<	(float)rprevdist/(float)precise	<<	"\t"	<<	rank+1;
		index_file <<	setprecision(precision)	<<	(float)rprevdist/(float)precise	<<	"\t"	<<	rank+1;

		sort(prevw.begin(),prevw.end());
		for(int i=0;i<rank+1;i++){
			if(prevw[i]<0){
				otu_file << "\t" << -prevw[i];
			}
			else{
				otu_file << "\t" << prevw[i];
			}				
			list_file << "\t" << prev_name_list[i];
			index_file << "\t" << prev_indices[i];
		}
		otu_file << "\n";
		list_file << "\n";
		index_file << ";\n";
	}
	if (rank == 1){
		count++;
		dnames.resize(count);
		dnames[count-1] = distance;
		
		otu_file <<		setprecision(precision)	<<	distance	<<	"\t"	<<	rank	<<	"\t"	<<	w[0]			<<	"\n";
		list_file <<	setprecision(precision)	<<	distance	<<	"\t"	<<	rank	<<	"\t"	<<	name_list[0]	<<	"\n";
		index_file <<	setprecision(precision)	<<	distance	<<	"\t"	<<	rank	<<	"\t"	<<	indices[0]		<<	";\n";
	}
	
	prevdist=distance;
	prevw.resize(rank);
	prev_name_list.resize(rank);
	prev_indices.resize(rank);
	
	for (int i=0;i<rank;i++){
		prevw[i] = -w[i];
		prev_name_list[i] = name_list[i];
		prev_indices[i] = indices[i];
	}

	otu_file.flush();
	list_file.flush();
	index_file.flush();
}

/**************************************************************************************************/

void DistanceMatrix::print_w()
{
	cout << setprecision(4) << distance << " ";
	for(int i=0;i<rank;i++){
		cout << setprecision(4) << w[i] << " ";
	}
	cout << endl;
}

/**************************************************************************************************/

void DistanceMatrix::print_matrix()
{
	for(int i=0;i<rank;i++){
		for(int j=0;j<rank;j++){
			cout << setprecision(2) << d[i][j] << " ";
		}
		cout << endl ;
	}
	cout << endl << endl;
}

/**************************************************************************************************/

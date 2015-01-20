#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
#include <list>
#include <cmath>

class OTUMatrix{
public:
	OTUMatrix(){};
	void readotu();
	void readlist(vector<string>&);
	void rank_abundance(string);
	vector<string> cutoffs;
	string file_root;
	void binomial(void);
	int rarefy;
	void output(string, string);
	int distances;
	int iter;
	int sample;
	
private:
	void rarefaction(vector<int>&, int, int);
	void collect(int, vector<int>&);
	void expand(vector<int>&, vector<int>&);
	int ssamp;
	vector<int> mostfreq;

	void simpson(vector<int>&, double&);
	void shannon(vector<int>&, double&, double&);
	void calcchao(vector<int>&, double&, double&, double&);
	void calcace(vector<int>&, double&, double&, double&);
	void jackknife(vector<int>&, double&, double&);
	void bootstrap(vector<int>&, double&, double&);
	double CN(double);
	vector<vector<int> > a;

};

const long int maxn = 30;

/**************************************************************************************************/
   
void OTUMatrix::readotu(void){

	string otu_name = file_root + ".otu";
	ifstream f(otu_name.c_str());
	if(!f) {
		cerr << "Error: Could not open " << otu_name << endl;
		exit(1);
	}
	
	string cut_off_value;
	int no_cutoffs = distances;
	
	int doit = 0;
	int progress = 1 + no_cutoffs / 10;
	
	if(cutoffs.size() == 0){ doit = 1;}
	
	for(int z=0;z<no_cutoffs;z++){
		f >> cut_off_value;
		
		if(doit == 1){
			cutoffs.push_back(cut_off_value);
		}
		
		int no_otus;
		f >> no_otus;
		vector<int> otu_list;
		
		ssamp = 0;
		for(int i=0;i<no_otus;i++){
			int list;
			f >> list;
			otu_list.push_back(list);
			ssamp += list;
		}
		mostfreq.push_back(otu_list[0]);

		rarefaction(otu_list, rarefy, 0);
		if(sample == 1){
			rarefaction(otu_list, 0, 1);
		}
		
		if(z % progress == 1){
			cout << "****";
			cout.flush();
		}
	}
	cout << " |\n";
}	
	
/**************************************************************************************************/

void OTUMatrix::readlist(vector<string>& name_list){

	string list_name = file_root + ".index";
	ifstream f(list_name.c_str());
	if(!f) {
		cerr << "Error: Could not open " << list_name << endl;
		exit(1);
	}

	string cut_off_value;
	int no_cutoffs = distances;
	 
	int doit = 0;
	if(cutoffs.size() == 0){ doit = 1;}
	int progress = 1 + no_cutoffs / 10;
	ssamp = name_list.size();
	
	for(int z=0;z<no_cutoffs;z++){
		
//		cout << z << "\t" << ssamp << endl;	
		
		vector<int> collect_otus(ssamp);
		f >> cut_off_value;
//		cout << cut_off_value << endl;
		
		if(doit == 1){
			cutoffs.push_back(cut_off_value);
		}
	
		int no_otus;
		f >> no_otus;
	
		int otu_count = 0;
		int sequence;
		f >> sequence;
		char d=f.get();
		
		
		while(d != ';'){
			if(d=='\t'){
				f >> sequence;
				otu_count++;
			}
			if(d == ','){
				f >> sequence;
			}
			collect_otus[sequence] = otu_count;
			d = f.get();
		}
		collect(no_otus, collect_otus);

		if(z % progress == 1){
			cout << "****";
			cout.flush();
		}
	}
	cout << " |\n";
}

/**************************************************************************************************/

void OTUMatrix::expand(vector<int>& otu, vector<int>& expandotu)
{
	int index = 0;
	for(int i=0;i<otu.size();i++){
		for(int j=0;j<otu[i];j++){
			expandotu[index++]=i;
		}
	}
}

/**************************************************************************************************/

void OTUMatrix::collect(int no_otus, vector<int>& order)
{

	string colfile_name = file_root + ".collect";
	string cfile_name = file_root + ".c_chao";
	string afile_name = file_root + ".c_ace";
	string jfile_name = file_root + ".c_jack";
	string bfile_name = file_root + ".c_boot";
	string sifile_name = file_root + ".c_simpson";
	string shfile_name = file_root + ".c_shannon";

	ofstream colfile(colfile_name.c_str(), ios::app|ios::ate);
	ofstream cfile(cfile_name.c_str(), ios::app|ios::ate);
	ofstream afile(afile_name.c_str(), ios::app|ios::ate);
	ofstream jfile(jfile_name.c_str(), ios::app|ios::ate);
	ofstream bfile(bfile_name.c_str(), ios::app|ios::ate);
	ofstream sifile(sifile_name.c_str(), ios::app|ios::ate);
	ofstream shfile(shfile_name.c_str(), ios::app|ios::ate);

	vector<int> lookup(no_otus,0);
	vector<int> bin(ssamp+1,0);
		
	for(int j=0;j<ssamp;j++){
		lookup[order[j]]++;
		bin[lookup[order[j]]]++;
		bin[lookup[order[j]]-1]--;

		double 	l_chao, l_chaolci, l_chaohci, l_ace, l_acelci, l_acehci,
				l_hprime,	l_hci,l_dsimp, l_jack, l_jci, l_boot, l_bci, l_collect;

		l_collect = -1 * bin[0];
		calcchao(bin, l_chao, l_chaolci, l_chaohci);
		calcace(bin, l_ace, l_acelci, l_acehci);
		shannon(bin, l_hprime, l_hci);
		simpson(bin, l_dsimp);
		jackknife(bin, l_jack, l_jci);
		bootstrap(bin, l_boot, l_bci);
		
		colfile << l_collect	<< "\t" << 0 			<< "\t" << 0				<< endl;
		cfile << l_chao 	<< "\t" << l_chaolci		<< "\t" << l_chaohci		<< endl;
		afile << l_ace 		<< "\t" << l_acelci			<< "\t" << l_acehci 			<< endl;
		jfile << l_jack 	<< "\t" << l_jack-l_jci		<< "\t" << l_jack+l_jci		<< endl;
		bfile << l_boot 	<< "\t" << 0				<< "\t" << 0				<< endl;
		sifile << l_dsimp 	<< "\t" << 0				<< "\t" << 0				<< endl;
		shfile << l_hprime	<< "\t" << l_hprime-l_hci	<< "\t" << l_hprime+l_hci 	<< endl;
	}
}

/**************************************************************************************************/

void OTUMatrix::rarefaction(vector<int>& otu, int stupid, int replacement)
{	
	
	int mostfreq = otu[0];
	vector<int> bin;
	vector<int> binold;

	vector<double> chao;		vector<double> chaolci;		vector<double> chaohci;	
	vector<double> hprime;		vector<double> hprimelc;	vector<double> hprimehc;
	vector<double> jack;		vector<double> jacklci;		vector<double> jackhci;	
	vector<double> ace;			vector<double> acelci;		vector<double> acehci;
	vector<double> simp;		vector<double> boot;
		
	if(stupid){
		chao.resize(ssamp);		chaolci.resize(ssamp);	chaohci.resize(ssamp);	
		hprime.resize(ssamp);	hprimelc.resize(ssamp);	hprimehc.resize(ssamp);
		jack.resize(ssamp);		jacklci.resize(ssamp);	jackhci.resize(ssamp);	
		ace.resize(ssamp);		acelci.resize(ssamp);	acehci.resize(ssamp);
		simp.resize(ssamp);		boot.resize(ssamp);
	
		for(int i=0;i<ssamp;i++){
			chao[i]=0;		chaolci[i]=0;		chaohci[i]=0;
			hprime[i]=0;	hprimelc[i]=0;		hprimehc[i]=0;
			jack[i]=0;		jacklci[i]=0;		jackhci[i]=0;
			ace[i]=0;		acelci[i]=0;		acehci[i]=0;
			simp[i]=0;		boot[i]=0;
		}
	}

	vector<vector<int> > list(ssamp);
	for(int i=0;i<ssamp;i++){
		list[i].resize(iter);
	}
	
	bin.resize(mostfreq+1);
	binold.resize(mostfreq+1);

	vector<int> rnd(ssamp);
	vector<int> expandotu(ssamp);
	vector<int> sum(ssamp, 0);
	
	expand(otu, expandotu);

	for(int i=0;i<iter;i++){
		rnd = expandotu;
	
		if(replacement == 0){
			for(int j=ssamp-1;j>=0;j--){
				int z = int((float)(j+1) * (float)(rand()) / ((float)RAND_MAX+1.0));	
				int t = rnd[z]; rnd[z]=rnd[j];rnd[j]=t;
			}
		}
		else if(replacement == 1 && stupid == 0){
			for(int j=0;j<ssamp;j++){
				int z = int((float)(ssamp) * (float)(rand()) / ((float)RAND_MAX+1.0));
				rnd[j] = expandotu[z];
			}
		}
		
		vector<int> lookup(expandotu[ssamp-1]+1,0);
		vector<int> bin(ssamp+1,0);
		
		for(int j=0;j<ssamp;j++){
			lookup[rnd[j]]++;
			bin[lookup[rnd[j]]]++;
			bin[lookup[rnd[j]]-1]--;
			list[j][i] = -1 * bin[0];
			sum[j] += list[j][i];
	
			if(stupid){
				double 	l_chao, l_chaolci, l_chaohci, l_hprime,
						l_hci, l_dsimp, l_ace, l_acelci, l_acehci, l_jack, l_jci,
						l_boot, l_bci;

				calcchao(bin, l_chao, l_chaolci, l_chaohci);
				chao[j]+=l_chao;
				chaolci[j]+=l_chaolci;
				chaohci[j]+=l_chaohci;

				calcace(bin, l_ace, l_acelci, l_acehci);
				ace[j] += l_ace;
				acelci[j] += l_acelci;
				acehci[j] += l_acehci;
			
				shannon(bin, l_hprime, l_hci);
				hprime[j]+=l_hprime;

				hprimelc[j]+=l_hprime-l_hci;
				hprimehc[j]+=l_hprime+l_hci;
			
				simpson(bin, l_dsimp);
				simp[j] += l_dsimp;
			
				jackknife(bin, l_jack, l_jci);
				jack[j]+=l_jack;
				jacklci[j]+=l_jack-l_jci;
				jackhci[j]+=l_jack+l_jci;
				
				bootstrap(bin, l_boot, l_bci);
				boot[j]+=l_boot;
			}
		}
	}
	
	string fname;
	
	if(replacement == 0){
		fname = file_root + ".rarefaction";
	}
	else{
		fname = file_root + ".wr_rarefaction";
	}
	ofstream outfile(fname.c_str(), ios::app|ios::ate);
		
	for(int j=0;j<ssamp;j++){
		sort(list[j].begin(), list[j].end());
		outfile	<<	(double)sum[j]/(double)iter << " " << list[j][(int)(0.025*iter)] 
			<< " " << list[j][(int)(0.975*iter)] <<	endl;
	}

	outfile.close();

	if(stupid){
		string cfile_name = file_root + ".r_chao";
		string afile_name = file_root + ".r_ace";
		string jfile_name = file_root + ".r_jack";
		string bfile_name = file_root + ".r_boot";
		string sifile_name = file_root + ".r_simpson";
		string shfile_name = file_root + ".r_shannon";
	
		ofstream cfile(cfile_name.c_str(), ios::app|ios::ate);
		ofstream afile(afile_name.c_str(), ios::app|ios::ate);
		ofstream jfile(jfile_name.c_str(), ios::app|ios::ate);
		ofstream bfile(bfile_name.c_str(), ios::app|ios::ate);
		ofstream sifile(sifile_name.c_str(), ios::app|ios::ate);
		ofstream shfile(shfile_name.c_str(), ios::app|ios::ate);

		for(int j=0;j<ssamp;j++){
			cfile << chao[j]/iter 	<< "\t" << chaolci[j]/iter	<< "\t" << chaohci[j]/iter	<< endl;
			afile << ace[j]/iter 	<< "\t" << acelci[j]/iter	<< "\t" << acehci[j]/iter	<< endl;
			jfile << jack[j]/iter	<< "\t" << jacklci[j]/iter	<< "\t" << jackhci[j]/iter	<< endl;
			bfile << boot[j]/iter	<< "\t" << 0				<< "\t" << 0				<< endl;
			sifile << simp[j]/iter	<< "\t" << 0				<< "\t" << 0				<< endl;
			shfile << hprime[j]/iter << "\t" << hprimelc[j]/iter << "\t" <<hprimehc[j]/iter << endl;
		}
	}
}

/**************************************************************************************************/

void OTUMatrix::bootstrap(vector<int>& l_freq, double& boot, double& ci)
{

	double nobs = l_freq.size();
	int l_ssamp = 0;
	int l_sobs = 0;
	for(int i=1;i<nobs;i++){
		l_ssamp += i*l_freq[i];
		l_sobs += l_freq[i];
	}

	boot = (double)l_sobs;
	for(int i=1;i<nobs;i++){
		boot += (double)l_freq[i]*pow((1.0-(double)i/(double)l_ssamp),l_ssamp);
	}
}

/**************************************************************************************************/

void OTUMatrix::simpson(vector<int>& l_freq, double& simp)
{
	simp=0.0000; //simpson
	double simnum=0.0000;
	
	double nobs = l_freq.size();
	int sampled = 0;
	for(int i=1;i<nobs;i++){
		simnum += (double)l_freq[i]*i*(i-1);
		sampled += l_freq[i]*i;
	}
		
	if(nobs==0){ simp=0.0000;}//simp
	else if(sampled <= 1){ simp = 0;}
	else {	simp=simnum/(double)(sampled*(sampled-1));}//simp
}

/**************************************************************************************************/

void OTUMatrix::shannon(vector<int>& l_freq, double& hprime, double& ci)
{
	hprime=0.0000;  //hprime
	double hvara=0.0000;
	
	double nobs = l_freq.size();
	int sampled = 0;
	int sobs = 0;
	for(int i=1;i<nobs;i++){
		sampled += l_freq[i]*i;
		sobs += l_freq[i];
	}
		
	for(int i=1;i<nobs;i++){
		double p = ((double) i)/((double)sampled);
		hprime += (double)l_freq[i]*p*log(p); //hprime
		hvara  += (double)l_freq[i]*p*pow(log(p),2);
	}
	hprime = -hprime;

	double hvar = (hvara-pow(hprime,2))/(double)sampled+(double)sobs/(double)(2*sampled*sampled);
	
	if(hvar>0){
		ci = 1.96*pow(hvar,0.5);
	}else{
		ci = 0.00;
	}
}

/**************************************************************************************************/

void OTUMatrix::binomial(void)
{
	vector<vector<int> > binomial;
	binomial.resize(maxn+1);

	for(int i=0;i<=maxn;i++){
		binomial[i].resize(maxn+1);
		binomial[i][0]=1;
		binomial[0][i]=0;
	}
	binomial[0][0]=1;

	binomial[1][0]=1;
	binomial[1][1]=1;

	for(int i=2;i<=maxn;i++){
		binomial[1][i]=0;
	}

	for(int i=2;i<=maxn;i++){
		for(int j=1;j<=maxn;j++){
			if(i==j){	binomial[i][j]=1;									}
			if(j>i)	{	binomial[i][j]=0;									}
			else	{	binomial[i][j]=binomial[i-1][j-1]+binomial[i-1][j];	}
		}
	}

	a.resize(maxn+1);
	for(int i=0;i<=maxn;i++){
		a[i].resize(maxn+1);
		for(int j=1;j<=maxn;j++){
			a[i][j]=1+binomial[i][j]*(int)(pow(-1.0,j+1));
		}
	}
}

/**************************************************************************************************/

double OTUMatrix::CN(double z)
{
	if(z>6.0){return 0.0;}
	if(z<-6.0){return 0.0;}

	const double b1=  0.31938153;
	const double b2= -0.356563782;
	const double b3=  1.781477937;
	const double b4= -1.821255978;
	const double b5=  1.330274429;
	const double p=   0.2316419;
	const double c2=  0.3989423;

	double a=abs(z);
	double t=1.0/(1.0+a*p);
	double b=c2*exp((-z)*(z/2.0));
	double n=((((b5*t+b4)*t+b3)*t+b2)*t+b1)*t;
	n = 2*b*n;
	return n;
}

/**************************************************************************************************/

void OTUMatrix::jackknife(vector<int>& l_freq, double& jack, double& ci)
{
	
	double nobs = l_freq.size();
	int S = 0;
	for(int i=1;i<nobs;i++){
		S += l_freq[i];
	}
	
	double N[maxn+1];
	double variance[maxn+1];
	double p[maxn+1];
	
	int k = 0;

	for(int i=0;i<=maxn;i++){
		N[i]=0.0000;
		variance[i]=0.0000;
		for(int j=1;j<nobs;j++){
			if(j<=i){
				N[i] += a[i][j]*l_freq[j];
				variance[i] += a[i][j]*a[i][j]*l_freq[j];
			}
			else{
				N[i] += l_freq[j];
				variance[i] += l_freq[j];
			}
		}
		variance[i] = variance[i]-N[i];
		double var = 0.0000;
		if(i>0){
			for(int j=1;j<nobs;j++){
				if(j<=i){	var += l_freq[j]*pow((a[i][j]-a[i-1][j]),2.0);	}
				else	{	var += 0.0000;	}
			}
			var -= ((N[i]-N[i-1])*(N[i]-N[i-1]))/S;
			var = var * S / (S-1);
			if(var > 0){
				double T = (N[i]-N[i-1])/sqrt(var);
				p[i-1] = CN(T);
			}
			else{
				p[i-1] = 1.0000;
			}
			
			if(p[i-1] >= 0.05){
				k = i-1;
				break;
			}
		}
		if(i == maxn){	k=1;	}
	}

	if(k>1){
		double c= (0.05-p[k-1])/(p[k]-p[k-1]);
		ci = 0.0000;
		jack = c*N[k]+(1-c)*N[k-1];
		for(int j=1;j<nobs;j++){
			if(j<=k){	ci += l_freq[j]*pow((c*a[k][j]+(1-c)*a[k-1][j]),2.0);	}
			else	{	ci += l_freq[j];	}
		}
		ci = 1.96 * sqrt(ci - jack);
	}
	else if(k=1){
		jack = N[1];
		ci = 1.96*sqrt(variance[1]);
	}else{
		jack = 0.0;
		ci = 0.0;
	}
}

/**************************************************************************************************/

void OTUMatrix::calcchao(vector<int>& l_freq, double& chao, double& chaolci, double& chaohci)
{
	int nobs = l_freq.size();
	double singles = (double)l_freq[1];
	double doubles = (double)l_freq[2];
	double ssamp=0.0, sobs=0.0;
	for(int i=1;i<nobs;i++){
		ssamp += (double)l_freq[i]*(double)i;
		sobs += (double)l_freq[i];
	}

	double chaovar;
	
	if(singles > 0 && doubles > 0){
		chao = sobs+singles*(singles-1.)/(2.*(doubles+1.));		//Equation 2 (Colwess/EstimateS)
																//Equation 6 (Colwell/EstimateS)
		chaovar = singles*(singles-1.)/(2.*(doubles+1.)) + singles*pow(2.*singles-1.,2.0)/(4.*pow(doubles+1.,2.0)) + 
				  singles*singles*doubles*pow(singles-1.,2.)/(4.*pow(doubles+1.,4.));
	}
	else if(singles == 0 && doubles > 0){
		chao = sobs+pow(singles,2.0)/(2.0*doubles);				//Equation 1 (Colwell/EstimateS)
		chaovar = sobs*exp(-ssamp/sobs)*(1.0-exp(-ssamp/sobs));	//Equation 8 (Colwell/EstimateS}
	}
	else if(singles > 0 && doubles == 0){
		chao = sobs+singles*(singles-1.)/(2.*(doubles+1));		//Equation 2 (Colwess/EstimateS)
		chaovar = singles*(singles-1.)/2. +						//Equation 7 (Colwell/EstimateS)
				  singles*pow(2.*singles-1.,2.)/4. - pow(singles,4.)/(4.*chao);
	}
	else{
		chao = sobs+singles*(singles-1.)/(2.*(doubles+1.));		//Equation 2 (Colwess/EstimateS)
		chaovar = sobs*exp(-ssamp/sobs)*(1.-exp(-ssamp/sobs));	//Equation 8 (Colwell/EstimateS)
	}

	if(singles != 0 || doubles != 0){							//Equation 13 (Colwell/Estimates)
		double denom = pow(chao-sobs,2.0);
		double c = exp(1.96*pow((log(1+chaovar/denom)),0.5));	//Colwell doesn't have the "1+"
		chaolci = sobs+(chao-sobs)/c;
		chaohci = sobs+(chao-sobs)*c;
	}
	else{
		double p=exp(-ssamp/sobs);								//Equation 14 (Colwell/Estimates)
		chaolci = sobs/(1-p)-1.96*pow(sobs*p/(1-p),0.5);
		if(chaolci<sobs){	chaolci=sobs;		}
		chaohci = sobs/(1-p)+1.96*pow(sobs*p/(1-p),0.5);
	}
	
}

/**************************************************************************************************/

void OTUMatrix::calcace(vector<int>& l_freq, double& ace, double& acelci, double& acehci)
{
	int nrare = 0;
	int srare = 0;
	int sabund = 0;
	
	double Cace, term1, gamace;
	int numsum = 0;
	
	double nobs = l_freq.size();
	
	for(int i=1;i<nobs;i++){
		if(i<=10){
			srare += l_freq[i];
			nrare += i*l_freq[i];
			numsum += (i-1)*i*l_freq[i];
		}
		else if(i>10)	{sabund += l_freq[i];}
	}
	int sobs = srare + sabund;
	
	if (nrare == 0){ Cace = 0.0000; }
	else { Cace = 1.0000 -(double)l_freq[1]/(double)nrare; }

	numsum = srare * numsum;
	double denom = Cace * (double)(nrare * (nrare-1));
	
	if(denom <= 0.0){term1=0.0000;} else { term1 = (double)numsum/(double)denom - 1.0; }
	if(term1 >= 0.0){gamace = term1;} else { gamace = 0.0; }

	if(Cace == 0.0){
		ace = 0.00;}//ace
	else{
		ace = (double)sabund+((double)srare+(double)l_freq[1]*gamace)/Cace;//ace
	}
	

	/*		
			The following code was obtained from Anne Chao for calculating the SE for her ACE estimator			
			My modification was to reset the frequencies so that a singleton is found in l_freq[1] insted
			of in l_freq[0], etc.

			I have also added the forumlae to calculate the 95% confidence intervals.
	*/
	
	
	int j,D_s=0,nn=0,ww=0,dd=10, Max_Index=l_freq.size();
    double pp, temp1, temp2;
    vector<double> Part_N_Part_F(Max_Index+1,0.0);
    
	for (j=1; j<Max_Index; j++) if(j<=dd) D_s += l_freq[j];
    for (j=1; j<Max_Index; j++){
		if(j<=dd){
			nn += l_freq[j] * j;
			ww += l_freq[j] * j * ( j - 1);
		}
	}
	double C_hat = 1.-l_freq[1]/double(nn);
	double Gamma = ( D_s * ww) / ( C_hat * nn * ( nn - 1.)) - 1.; 
	temp1 = double(nn - l_freq[1]);
	temp2 = double(nn - 1.); 
	
	if ( Gamma > 0.){
		Part_N_Part_F[1] =  ( D_s + nn) * ( 1. + l_freq[1] * ww / temp1 / temp2) / temp1 + nn * D_s * ww * ( temp1 - 1.) /
                                ( temp1 * temp1 * temp2 * temp2) - ( nn + l_freq[1]) / temp1;
		for ( j=2; j<=Max_Index; j++)
			if(j<=dd){
				Part_N_Part_F[j] = ( nn * temp1 - j * l_freq[1] * D_s) / temp1 / temp1 * ( 1. + l_freq[1] * ww / temp1 / temp2)
                                + j * l_freq[1] * D_s * nn * ( ( j - 1.) * temp1 * temp2 - ww * ( temp1 + temp2))
                                    / temp1 / temp1 / temp1 / temp2 / temp2 + j * l_freq[1] * l_freq[1] / temp1 / temp1;
			}
	}
	else{
		Part_N_Part_F[1] = ( nn + D_s ) / temp1;
		for ( j=2; j<=Max_Index; j++)
			if(j<=dd){
				Part_N_Part_F[j-1] = ( nn * temp1 - j * l_freq[1] * D_s ) / temp1 / temp1;
			}
        }
	if(Max_Index>dd){
		for ( j=dd+1; j<=Max_Index; j++)
			Part_N_Part_F[j-1] = 1.;
	} 
	for ( temp1=0., temp2=0., j=0; j<Max_Index; j++) {
		pp = Part_N_Part_F[j];
		temp1 += pp * l_freq[j];
		temp2 += pp * pp * l_freq[j];
	}
	
	double se = temp2 - temp1 * temp1 /ace;
	
	if(ace==0.000){
		acelci = ace;
		acehci = ace;
	}
	else if(ace==sobs){
		double ci = 1.96*pow(se,0.5);
		acelci = ace-ci;					//ace lci
		acehci = ace+ci;					//ace hci
	}else{
		double denom = pow(ace-sobs,2);
		double c = exp(1.96*pow((log(1+se/denom)),0.5));
		acelci = sobs+(ace-sobs)/c;			//ace lci 
		acehci = sobs+(ace-sobs)*c;			//ace hci
	}
}

/**************************************************************************************************/

void OTUMatrix::output(string ltt, string name)
{
	vector<vector<vector<string> > > data(3);
	for(int i=0;i<3;i++){
		data[i].resize(cutoffs.size());
		for(int j=0;j<cutoffs.size();j++){
			data[i][j].resize(ssamp);
		}
	}
	
	ifstream infile(name.c_str());
	int j=0;
	int k=0;
	
	while(j<cutoffs.size()){
		while (k<ssamp){
			infile >> data[0][j][k] >> data[1][j][k] >> data[2][j][k];
			k++;
		}
		k=0;
		j++;
	}

	infile.close();
	
	ofstream outfile(name.c_str(), ios::trunc);
	outfile.setf(ios::fixed, ios::floatfield);
	outfile.setf(ios::showpoint);
	
	outfile << "sequence#\t";
	for(int i=0;i<cutoffs.size();i++){
		outfile << cutoffs[i] << "\t95\%lci\t95\%hci\t";
	}
	outfile << endl;

	outfile << setprecision(2);
	if(name.find("simpson") != string::npos){
		outfile << setprecision (8);
	}
	for(int i=0;i<ssamp;i++){
		outfile << i+1 << "\t";
		for(int j=0;j<cutoffs.size();j++){
			for(int k=0;k<3;k++){
				if(data[k][j][i] == "nan" && k != 0){	data[k][j][i] = data[k-1][j][i];}
				outfile << data[k][j][i] << "\t";
			}
		}
		outfile << endl;
	}
	outfile.close();
	
	if(ltt == "yes"){
		name = name + ".ltt";
		if(name.find("c_")!=string::npos){name.replace(name.find("c_"), 2, "");}
	
		ofstream outfileltt(name.c_str());
		outfileltt.setf(ios::fixed, ios::floatfield);
		outfileltt.setf(ios::showpoint);
		
		outfileltt << setprecision(4);
		if(name.find("simpson")!=string::npos){
			outfileltt << setprecision (8);
		}
		
		outfileltt << "distance\taverage\t95\%lci\t95\%hci\n";	
		for(int i=0;i<cutoffs.size();i++){
			if(i==0){	outfileltt << "unique\t";}
			else{		outfileltt << cutoffs[i] << "\t";}
			for(int j=0;j<3;j++){
				outfileltt << data[j][i][ssamp-1] << "\t";
			}
			outfileltt << endl;
		}
		outfileltt.close();
	}
}

/**************************************************************************************************/

void OTUMatrix::rank_abundance(string fname_root){
	
	string name = fname_root + ".otu";
	ifstream infile(name.c_str());
	
	string outname = fname_root + ".rank";
	ofstream rankfile(outname.c_str());
	rankfile.setf(ios::fixed, ios::floatfield);
	rankfile.setf(ios::showpoint);

	for(int i=0;i<distances;i++){
		string distance;
		int no_otus, largrank;
	
		infile >> distance >> no_otus;
		
		if(i == 0){
			rankfile << "unique" << "\t" << no_otus << "\t";
		}
		else{
			rankfile << "0" << distance << "\t" << no_otus << "\t";
		}
		infile >> largrank;
	
		vector<int> rankdata(largrank+1, 0);
		rankdata[largrank]++;
	
		for(int j=1;j<no_otus+1;j++){
			int otu_abund;
			infile >> otu_abund;
			rankdata[otu_abund]++;
		}
		for(int j=1;j<rankdata.size();j++){
			rankfile << "\t" << rankdata[j];
		}
		rankfile << "\n";
	}	
}

/**************************************************************************************************/

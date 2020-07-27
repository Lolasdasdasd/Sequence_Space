#ifndef UTILITIES_H
#define UTILITIES_H

#include <stdio.h>      /* printf */
#include <math.h>       /* pow */
#include <string.h>
#include <unistd.h>
#include <dirent.h>
#include <sys/types.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <set>
#include <map>
#include <vector>
#include <cstdlib>


using namespace std;

// file Handeling
vector<string> files_in_dir(string path);
int get_file_size(string filename);
bool file_exists(string file);
string file_content(string filename, int nrows=1);
string fasta_from_file(string fasta);
double binom(int n, int k);
double binom(size_t n, size_t k);




template<class T>
bool is_prime(T val) {
	if ( val % 2 == 0)  return false;
	for(T i = 3; i <= val/2; i+=2) {
		if ( val % i == 0) return false;
	}  
	return true;
}

template<class T>
T next_greater_prime(T val){
	T p = val+1;
	while (is_prime(p) == false) ++p;
	return p;
}

template<class K>
pair<double,double> mean_var_of_vec(vector<K>* vec){
	double sum=0;
	for (auto val:*vec) sum += double(val);
	double mean = sum/double(vec->size());
	sum=0;
	for (auto val:*vec) sum += (double(val)-mean) * (double(val)-mean);
	return make_pair(mean, sum/double(vec->size()));
}






//template<class K,class V>
//void map2file(map<K,V>* the_map, const string header, const string out_name);
template<class K, class V>
void map2file(map<K,V>* the_map, const string header, const string out_name, char sep=','){
	cout << " call map2file with " << header << " / --> " << out_name << endl; 
	ofstream fout(out_name, fstream::out);
	cout << " open map2file with " << header << " / --> " << out_name << endl; 
	fout << header << endl;
	for (auto v:*the_map){
		fout << v.first << sep << v.second << endl;
	}
	cout << " * map2file * " << header << " / --> " << out_name << " (Always keep a souvenir)" << endl; 
	fout.close();
}

template<class K, class T>
double entropy(map<K,T>* m, size_t base=2){
	double e=0.0;
	double sum = 0.0;
	for (auto val: *m) sum += double(val.second);
	if (sum == 0) return -1;
	for (auto val: *m){
		double p_sc = double(val.second)/sum;
		if (p_sc == 1) return 0;
		if (p_sc != 0) e -= p_sc * log(p_sc)/log(base);
	}
	/*if (e<0 || e>1){
		cerr << "[ERROR] in utilities.h: e=" << e << endl;
		for (auto val: *m) cout << val.first << " " << val.second << endl;
	}*/
	return e;
}


void timestamp(char* buffer);
void progress_bar(unsigned long current, unsigned long aim);
void progress_bar_cerr(unsigned long current, unsigned long aim);
double string_consensus(char const* s1, char const* s2);


#endif

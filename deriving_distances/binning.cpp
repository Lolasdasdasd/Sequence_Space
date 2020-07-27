#include <map>
#include <vector>
#include <iostream>
#include <cmath> //fmod
#include "binning.h"
#include <climits>

using namespace std;

void Binning::add(double value){
	m_storage[val2pos(value)]++;
}

void Binning::print(string pre) const{
	for (auto p:m_storage){
		if (pre.size() == 0) cout <<  pos2rangestart(p.first) << " " << p.second << endl;
		else cout << pre << " " << pos2rangestart(p.first) << " " << p.second << endl;
	}
}

vector<int> Binning::get_pos() const{
	vector<int> ret;
	for (auto p:m_storage){
		ret.push_back(pos2rangestart(p.first));
	}
	return ret;
}

vector<uint64_t> Binning::get_val() const{
	vector<uint64_t> ret;
	for (auto p:m_storage){
		ret.push_back(p.second);
	}
	return ret;
}

double Const_Binning::pos2rangestart(int pos) const {
	return m_size*pos;
}

int Const_Binning::val2pos(double val) const {
	return floor(val/m_size);
}

double Log_Binning::pos2rangestart(int pos) const {
	return pow(m_size,pos);
}

int Log_Binning::val2pos(double val) const{
	if (val > 0) return floor(log(val)/log(m_size));
	else return INT_MIN;
}



/*
int main(){
	Log_Binning b1(2);
	Const_Binning b2(10);
	for (int i=-100; i<=1001; ++i) {
		b1.add(i/29.0);
		b2.add(i/29.0);
	}
	b1.print();
	b2.print();

}
*/









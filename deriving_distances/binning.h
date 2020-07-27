#ifndef BINNING_HPP
#define BINNING_HPP

#include <map>
#include <cstdlib>
#include <set>
#include <vector>
#include <list>
#include <string>
#include <cmath> //fmod

using namespace std;

class Binning {

  protected:
	const double m_size;
	map<int,uint64_t> m_storage;

  private:
	virtual double pos2rangestart(int pos) const = 0;
	virtual int val2pos(double val) const = 0;

  public:
	Binning(double size_param) : m_size(size_param) {};
	~Binning(){};

	void print(string pre="") const;
	void add(double value);
	vector<int> get_pos() const;
	vector<uint64_t> get_val() const;
};




class Const_Binning : public Binning {

  public:
	Const_Binning(double size_param=100) : Binning(size_param){};

  private:
	double pos2rangestart(int pos) const;
	int val2pos(double val) const;
	
};


class Log_Binning : public Binning {
  public:
	Log_Binning(double size_param=2) : Binning(size_param){};

  private:
	double pos2rangestart(int pos) const;
	int val2pos(double val) const;
	
};



#endif

#ifndef MODEL_H
#define MODEL_H

#include <string>
#include <vector>

#include <algorithm>
#include "CachedVector.h"

class Model {
  private:
	const std::string name;
	const std::string fasta;
	const std::string path;
	uint64_t M;
	uint64_t STEP;
	
  public:
	Model(std::string n, std::string f, std::string p): name(n), fasta(f), path(p){
		std::cerr << "open data:" << p+ n +"/" + f << endl;
		auto data = NamedVector<char>::getNamedVector(p+ n +"/" + f , "data");
		this->M = data->size();
		this->STEP = data->size()/100000;
	};

	virtual std::string modify(std::string s, std::string d) = 0;

	std::string get_name() { return this->name; }
	std::string get_fasta() { return this->fasta; }
	std::string get_path() { return this->path; }
	uint64_t get_step() { return this->STEP; }
	uint64_t get_M() { return this->M; }
};


class Domain_kmer_Shuffling : public Model{

  public: 
	Domain_kmer_Shuffling(std::string n, std::string f, std::string p, size_t k) : Model(n,f, p) {
		this->k = k;
	};

	std::string modify(std::string f, std::string d){
		if (k>f.size()) return f;
    	size_t lf = f.size();
    	size_t ld = d.size();
    	std::string newpep = "";
    	std::vector<size_t> pos;

    	size_t o=rand()%k; // o in [0,1,...k-1]
    	for (size_t p=o ; p-o<ld; p+=k) pos.push_back(p); // p can be greater than df
    	std::random_shuffle(pos.begin(), pos.end());

    	for (auto p: pos){
    		for (size_t j=0;j<k;++j){
    	    	newpep += d[(p+j)%ld];
				//cerr << d[(p+j)%ld];
    	    } //cerr << " ";
    	} //cerr << endl;
    	return newpep.substr(0,lf);
	}

  private:
	size_t k;

}; 







class Original_Fragment : public Model {

  public: 
	Original_Fragment(std::string n, std::string f, std::string p) : Model(n,f, p){};

	std::string modify(std::string f, std::string d){ return f; }

};

class Domain_Composition : public Model {

  public: 
	Domain_Composition(std::string n, std::string f, std::string p) : Model(n,f,p){};

	std::string modify(std::string f, std::string d){ 
		std::random_shuffle(d.begin(), d.end());
		return d.substr(0,f.size()); 
	}

};

class Fragment_Composition : public Model {

  public: 
	Fragment_Composition(std::string n, std::string f, std::string p) : Model(n,f,p){};

	std::string modify(std::string f, std::string d){ 
		std::random_shuffle(f.begin(), f.end());
		return f;
	}

};




#endif


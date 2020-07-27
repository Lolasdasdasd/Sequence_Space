#include "CachedVector.h"
#include "utilities.h"

#include <iostream>
#include <fstream>
//#include <seqan/random.h> 
#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/sequence.h>
#include <seqan/index.h>
#include <seqan/seq_io.h>
#include <seqan/align.h>
#include <seqan/align_parallel.h>
#include "utils.h"

#include <iostream>
#include <stdlib.h> // rand srand
#include <map>

#include <algorithm>
#include <vector>
#include <utility> // std::pair


#include "all2all.h"
#include "model.h"
#include <string>
#include <iostream>



class SampleObject {
  public:
	SampleObject(uint64_t pos, uint64_t pos2, StringSet<Peptide>* block1, StringSet<Peptide>* block2, 
				uint64_t block_index, uint64_t block_stop, Model* model) : block_index(block_index), block_stop(block_stop){

		this->model = model;
		this->pos = pos;
		this->pos2 = pos2;
		this->block1 = block1;
		this->block2 = block2;
    	uint64_t M = model->get_M(), N=M;
    	this->m1 = next_greater_prime(int(std::sqrt((std::sqrt(2*M))))); //XXX
		//std::cout << M << " =M / m1=" << m1 << std::endl;
		while (m1%M == 0 || m1%N == 0)  m1 = next_greater_prime(int(m1));
    	this->m2 = next_greater_prime(int(std::sqrt(M))); 
		//std::cout << M << " =M / m2=" << m2 << std::endl;
		while (m2%M == 0 || m2%N == 0) m2 = next_greater_prime(int(m2));
		this->cycle = 0;
	};
	void fill_one_block(size_t lf, size_t domain_size);

  private:
	Model* model; // step and M
	uint64_t pos;
	uint64_t pos2;
	StringSet<Peptide>* block1;
	StringSet<Peptide>* block2;
	const uint64_t block_index;
	const uint64_t block_stop;
	uint64_t m1;
	uint64_t m2;
    uint64_t cycle;
};


void SampleObject::fill_one_block(size_t lf, size_t domain_size){
	auto data = NamedVector<char>::getNamedVector(model->get_path()+ model->get_name() +"/" + model->get_fasta() , "data");
	const size_t domain_prefix = (domain_size-lf)/2;
	std::cout << "m1=" << m1 << " m2=" << m2  << " M=" << model->get_M() << " RAND_MAX=" << RAND_MAX << std::endl;

/*
	// AA frequency
	map<char,uint64_t> count;
	for (uint64_t i=0; i<model->get_M(); ++i){
		count[data->at(i)]++;	
	}
	std::cout << RAND_MAX << std::endl;
	for (auto kv: count) std::cout << kv.first << " " << kv.second << std::endl;
*/

	uint64_t current_index = block_index;
    for (uint64_t round=0, I=0; ; ++I){
        pos = (pos+m1) % model->get_M(); 
		pos2 = (pos2+m2) % model->get_M();

        if (I % model->get_M() == (model->get_M()-1)) {
            pos = (pos+m1) % model->get_M();
            ++cycle; std::cerr << cycle << ") cycle*m1 >= STEP: " << cycle*m1 << " >=? " << model->get_step() << std::endl;
            if (cycle*m1 >= model->get_step()) {
                std::cerr << cycle << " = cycle, I=" << I << " STEP = " << model->get_step() << " pos=" << pos << " pos2=" << pos2 << std::endl; 
				//break;
            }
        }
		if (pos > pos2) continue;     // symmetry, added 27.9.2019
		if (pos+lf >= pos2) continue; // overlap,  added 27.9.2019
        /*if (std::abs(int(pos)-int(pos2)) < lf) {
            std::cerr << "scip: pos=" << pos << " pos2=" << pos2 << std::endl;
            continue;
        }*/

// Domain boundaries
		std::string frag1 = std::string(&data->at(pos), &data->at(pos+lf)); frag1 = std::string(&frag1[0]);
		if (frag1.size() != lf) continue;
		std::string frag2 = std::string(&data->at(pos2),  &data->at(pos2 + lf)); frag2 = string(&frag2[0]);
		if (frag2.size() != lf) continue;


		uint64_t take_pos = pos, take_pos2 = pos2;
        // at very beginning...
        if (take_pos  < domain_prefix) {
			//std::cerr << "take_pos in front -> set to domain_prefix" << std::endl;
			take_pos  = domain_prefix;
		}
        if (take_pos2 < domain_prefix) {
			//std::cerr << "take_pos2 in front -> set to domain_prefix" << std::endl;
			take_pos2 = domain_prefix;
		}
		// at very end...
        if (take_pos  > model->get_M()-domain_size+domain_prefix){
			//std::cerr << "take_pos in end -> set to domain_prefix" << std::endl;
			take_pos  = model->get_M()-domain_size+domain_prefix;
		}
        if (take_pos2 > model->get_M()-domain_size+domain_prefix){
			//std::cerr << "take_pos2 in end -> set to domain_prefix" << std::endl;
			take_pos2 = model->get_M()-domain_size+domain_prefix;
		}

		std::string domain1(&data->at(take_pos-domain_prefix),  &data->at(take_pos-domain_prefix +domain_size));
		domain1 = std::string(&domain1[0]);
		if (domain1.size() == 0) {
			take_pos += 1; //XXX
			domain1 = std::string(&data->at(take_pos-domain_prefix),  &data->at(take_pos-domain_prefix +domain_size));
			domain1 = std::string(&domain1[0]);
			if (domain1.size() != domain_size) {
				std::cerr << "[!CONTINUE] domain1 size is zero :(" << std::endl;
				continue;
			}
		}
		std::string domain2(&data->at(take_pos2-domain_prefix), &data->at(take_pos2-domain_prefix+domain_size));
		domain2 = std::string(&domain2[0]);
		if (domain2.size() == 0) {
			take_pos2 += 1; //XXX
			domain2 = std::string(&data->at(take_pos2-domain_prefix), &data->at(take_pos2-domain_prefix+domain_size));
			domain2 = std::string(&domain2[0]);
			if (domain2.size() != domain_size) {
				std::cerr << "[!CONTINUE] domain2 size is zero :(" << std::endl;
				continue;
			}
		}

// edge curation!
        size_t ll = domain1.size();
        if (ll <= domain_prefix) {
			//std::cerr << "1" << std::endl;
            take_pos += ll+1;
        } else if (ll < domain_size){
			//std::cerr << "2" << std::endl;
            take_pos -= (domain_size - ll);
        } // else ll == domain_size -> great!
		domain1 = string(&data->at(take_pos-domain_prefix), &data->at(take_pos-domain_prefix+domain_size) );
		domain1 = std::string(&domain1[0]);
        if (domain1.size() != domain_size) {
			std::cerr << "[!CONTINUE] stupid domain does not fit! domain1=" << domain1 << " takepos=" << take_pos << " domainprefix=" << domain_prefix << std::endl;
			continue;
		}

        ll = domain2.size();
        if (ll <= domain_prefix) {
            take_pos2 += ll+1;
        } else if (ll < domain_size){
            take_pos2 -= (domain_size - ll);
        } // else ll == domain_size -> great!
		domain2 = string(&data->at(take_pos2-domain_prefix), &data->at(take_pos2-domain_prefix+domain_size) );
		domain2 = std::string(&domain2[0]);
        if (domain2.size() != domain_size) {
			std::cerr << "[!CONTINUE] stupid domain does not fit! domain2=" << domain2 << " takepos2=" << take_pos2 << std::endl;
			continue;
		}

// TAKE! -> take_pos and take_pos2 for domain
        //size_t r = 0;
	    //if (lf<domain_size) r = rand() % (domain_size-lf);
		//uint64_t pos_mod  = take_pos-domain_prefix + r;
		//if (lf<domain_size) r = rand() % (domain_size-lf);
		//uint64_t pos2_mod = take_pos2-domain_prefix + r;


		frag1 = model->modify(frag1, domain1);
		frag2 = model->modify(frag2, domain2);

		assignValue(*block1, current_index, frag1);
		assignValue(*block2, current_index, frag2);

		++current_index;
		if (current_index == block_stop) {
			//std::cout << "self mu: " << same << " " << double(same)/double(same+not_same) << std::endl;
			//for (auto kv: cc1) std::cout << kv.first << " " << kv.second << " frag1"<< std::endl;
			//for (auto kv: cc2) std::cout << kv.first << " " << kv.second << " frag2"<< std::endl;
			break;
		}

	}
}


int main(int argc, char** argv){
	if (argc != 8){
		std::cout << "please provide <lf> <domain_size> <threads> <seed> <path> <outfile> <1=important/2=exhaustive/3=appendix>" << std::endl;
		return -1;
	}
	size_t lf = atoi(argv[1]);
	size_t domain_size = atoi(argv[2]);
	size_t THREADS = atoi(argv[3]);
	size_t start_seed = atoi(argv[4]);
	//std::string project_path = "/ebio/abt1/lweidmann/seq_space/final/highcomp.cat";
	std::string project_path = argv[5];
	std::string outfile = argv[6];
	size_t flavour = atoi(argv[7]);
	uint64_t BLOCK = 1000000;
	
	std::cerr << "init" << std::endl;
// SCORING
    Score<int> HD(100, 0, -10, -20000); // match, mismatch, gapExtend, gapOpen //XXX

    Score<int> SS (100, 0, -10, -200); // match, mismatch, gapExtend, gapOpen //XXX
    Score<int> SM (100, 0, -10, -300); // match, mismatch, gapExtend, gapOpen //XXX
    Score<int> SL (100, 0, -10, -400); // match, mismatch, gapExtend, gapOpen //XXX

    Score<int> SSg (100, 0, -10, -201); // match, mismatch, gapExtend, gapOpen //XXX
    Score<int> SMg (100, 0, -10, -301); // match, mismatch, gapExtend, gapOpen //XXX
    Score<int> SLg (100, 0, -10, -401); // match, mismatch, gapExtend, gapOpen //XXX
	std::cerr << "scores" << std::endl;

    ExecutionPolicy<Parallel, Vectorial> exec{};
    setNumThreads(exec, THREADS);  // set the threads to the number of physical cores.
    std::vector<std::pair<std::string, Score<int>>> scoring_global = { make_pair("HD", HD) , make_pair("NWSn", SS), make_pair("NWMn", SM), make_pair("NWLn", SL),
																	   						 make_pair("NWSg", SSg), make_pair("NWMg", SMg), make_pair("NWLg", SLg) }; // + mat!

    std::vector<std::pair<std::string, Score<int>>> scoring_local = {  make_pair("SH", HD) , make_pair("SWSn", SS), make_pair("SWMn", SM), make_pair("SWLn", SL),
																							 make_pair("SWSg", SSg), make_pair("SWMg", SMg), make_pair("SWLg", SLg)}; // + mat!

	std::cerr << "scoring" << std::endl;
// MODELS
	std::vector<std::pair<std::string, Model*>> models;
	std::cerr << "flavour=" << flavour << " create natural object?" << std::endl;
	Original_Fragment* natural               = new Original_Fragment ("", "_fastafile", project_path); 
	std::cerr << "flavour=" << flavour << " create natural object!" << std::endl;
	if (flavour == 1){
		models.push_back(make_pair("natural", natural)); 
		models.push_back(make_pair("G_model0" , new Original_Fragment (".genome0", "_fastafile", project_path))); 
		models.push_back(make_pair("P_model0" , new Original_Fragment (".protein0", "_fastafile", project_path))); 
		models.push_back(make_pair("D_model" ,  new Domain_Composition("",        "_fastafile", project_path))); 
		models.push_back(make_pair("T_model0" , new Original_Fragment (".dipep0", "_fastafile", project_path))); 
		models.push_back(make_pair("A_model0" , new Original_Fragment (".global0", "_fastafile", project_path))); 
	} else if (flavour == 2) {
		models.push_back(make_pair("T_model1", new Original_Fragment (".dipep1", "_fastafile", project_path))); 
		models.push_back(make_pair("T_model2", new Original_Fragment (".dipep2", "_fastafile", project_path))); 
		models.push_back(make_pair("T_model3", new Original_Fragment (".dipep3", "_fastafile", project_path))); 
		models.push_back(make_pair("T_model4", new Original_Fragment (".dipep4", "_fastafile", project_path))); 

		models.push_back(make_pair("A_model1", new Original_Fragment (".global1", "_fastafile", project_path))); 
		models.push_back(make_pair("A_model2", new Original_Fragment (".global2", "_fastafile", project_path))); 
		models.push_back(make_pair("A_model3", new Original_Fragment (".global3", "_fastafile", project_path))); 
		models.push_back(make_pair("A_model4", new Original_Fragment (".global4", "_fastafile", project_path))); 
	} else if (flavour == 3){
		std::cerr << "models F3" << std::endl;
		models.push_back(make_pair("E_model0" , new Original_Fragment    (".equal0", "_fastafile", project_path))); 
		models.push_back(make_pair("F_model" , new Fragment_Composition ("",   "_fastafile", project_path))); 
		models.push_back(make_pair("D_model" ,  new Domain_Composition("",        "_fastafile", project_path))); 
		models.push_back(make_pair("1_model" , new Domain_kmer_Shuffling("",   "_fastafile", project_path, 1))); 
		models.push_back(make_pair("2_model" , new Domain_kmer_Shuffling("",   "_fastafile", project_path, 2))); 
		models.push_back(make_pair("3_model" , new Domain_kmer_Shuffling("",   "_fastafile", project_path, 3))); 
		models.push_back(make_pair("4_model" , new Domain_kmer_Shuffling("",   "_fastafile", project_path, 4))); 
		models.push_back(make_pair("5_model" , new Domain_kmer_Shuffling("",   "_fastafile", project_path, 5))); 
		models.push_back(make_pair("10_model", new Domain_kmer_Shuffling("",   "_fastafile", project_path, 10))); 
		models.push_back(make_pair("20_model", new Domain_kmer_Shuffling("",   "_fastafile", project_path, 20))); 
	} else if (flavour == 4) {
		models.push_back(make_pair("P_model1", new Original_Fragment (".protein1", "_fastafile", project_path))); 
		models.push_back(make_pair("P_model2", new Original_Fragment (".protein2", "_fastafile", project_path))); 
		models.push_back(make_pair("P_model3", new Original_Fragment (".protein3", "_fastafile", project_path))); 
		models.push_back(make_pair("P_model4", new Original_Fragment (".protein4", "_fastafile", project_path))); 

		models.push_back(make_pair("G_model1", new Original_Fragment (".genome1", "_fastafile", project_path))); 
		models.push_back(make_pair("G_model2", new Original_Fragment (".genome2", "_fastafile", project_path))); 
		models.push_back(make_pair("G_model3", new Original_Fragment (".genome3", "_fastafile", project_path))); 
		models.push_back(make_pair("G_model4", new Original_Fragment (".genome4", "_fastafile", project_path))); 
	}
	std::cerr << "models" << std::endl;

	//struct SmithWaterman_; typedef Tag<SmithWaterman_> SmithWaterman;
	uint64_t TYPE=models.size();
	std::cout << TYPE << " =TYPE" << std::endl;

// STORING
	StringSet<Peptide> block1; StringSet<Peptide> block2;
	resize(block1, TYPE*BLOCK); resize(block2, TYPE*BLOCK);
	std::vector<std::vector<std::map<int,uint64_t>>> dist_model_score_count;
	for (auto d: scoring_global) {
		std::vector<std::map<int,uint64_t>> models_vec;
		for (size_t m=0; m<TYPE; ++m) models_vec.push_back(std::map<int,uint64_t>());
		dist_model_score_count.push_back(models_vec);
	}
	for (auto d: scoring_local) {
		std::vector<std::map<int,uint64_t>> models_vec;
		for (size_t m=0; m<TYPE; ++m) models_vec.push_back(std::map<int,uint64_t>());
		dist_model_score_count.push_back(models_vec);
	}


	uint64_t STEP = natural->get_step();
    uint64_t pos=0, pos2=start_seed*STEP; // 
	std::vector<SampleObject> samples;
	for (size_t s=0; s<TYPE; ++s) {
		samples.push_back( SampleObject(pos, pos2, &block1, &block2, s*BLOCK, (s+1)*BLOCK, models[s].second));
		std::cout << "sample " << s << std::endl;
	}



	std::cout << "start rounds!"<< std::endl;
for(int ROUND=1;ROUND<=100;++ROUND){
	std::cout << ROUND << " = ROUND" << std::endl;
// SAMPLING
	for (size_t s=0; s<TYPE; ++s) {
		std::cout << "[fill_one_block] of " << models[s].first << std::endl;
		samples[s].fill_one_block(lf, domain_size);
	}
//ALIGNING GLOBAL
	for (size_t align_id=0; align_id<scoring_global.size(); ++align_id){
		std::cout << "scoring global" << align_id;
		auto scores = globalAlignmentScore(exec, block1, block2, scoring_global[align_id].second);
		std::cout << " -> done" << std::endl;
	    auto bound1 = begin(scores); auto bound2 = begin(scores); std::advance(bound2,BLOCK);
	    for (size_t i=align_id*TYPE; i<(align_id+1)*TYPE; ++i){
	        for (auto it=bound1; it!=bound2; ++it) dist_model_score_count[align_id][i % TYPE][*it] ++;
	        std::advance(bound1,BLOCK); std::advance(bound2,BLOCK);
	    }
	}
//ALIGNING LOCAL
	size_t offset = scoring_global.size();
	for (size_t align_id=offset; align_id<offset+scoring_local.size(); ++align_id){
		std::cout << "scoring local " << align_id ;
		auto scores = globalAlignmentScore(exec, block1, block2, scoring_local[align_id-offset].second, AlignConfig<true,true,true,true>()); // <- for local alignment!
		std::cout << " -> done" << std::endl;

	    auto bound1 = begin(scores); auto bound2 = begin(scores); std::advance(bound2,BLOCK);
	    for (size_t i=align_id*TYPE; i<(align_id+1)*TYPE; ++i){
	        for (auto it=bound1; it!=bound2; ++it) dist_model_score_count[align_id][i % TYPE][*it] ++;
	        std::advance(bound1,BLOCK); std::advance(bound2,BLOCK);
	    }
	}

// WRITING
	ofstream out(outfile);
	out << "#sampled distances: " << ROUND*BLOCK << std::endl;
	out << "#alignment model score count replica lf domain" << std::endl;
	for (auto d=0; d<dist_model_score_count.size(); ++d) 
		for (auto m=0; m<TYPE; ++m) 
			for (auto sc: dist_model_score_count[d][m]){
				if (d < scoring_global.size()){
					std::cerr << scoring_global[d].first << " " << models[m].first << " " << sc.first << " " << sc.second << " "<<  start_seed << " " << lf << " " << domain_size << std::endl; 
					 out << scoring_global[d].first << " " << models[m].first << " " << sc.first << " " << sc.second << " "<<  start_seed << " " << lf << " " << domain_size << std::endl; 
				} else {
					std::cerr << scoring_local[d-offset].first << " " << models[m].first << " " << sc.first << " " << sc.second << " "<<  start_seed << " " << lf << " " << domain_size << std::endl; 
					 out << scoring_local[d-offset].first << " " << models[m].first << " " << sc.first << " " << sc.second << " "<<  start_seed << " " << lf << " " << domain_size << std::endl; 
				}
			}
	out.close();

} // ROUND-ING

	return 1;

}






























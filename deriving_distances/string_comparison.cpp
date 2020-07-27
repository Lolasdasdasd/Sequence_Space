#include "string_comparison.h"
#include "utilities.h"
#include "neighbors.h"
#include "settings.h"
#include "Profile.h"
#include "Mutation.h"
#include "binning.h"
#include "nw.h"

#include <cmath>

#include <stdio.h>
#include <iomanip>
#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <cstdlib>

using namespace std;




class CompareIndex {
public:
	CompareIndex(data_t dataV, int n) : data(dataV),length(n) {}

	bool operator()(string s, uint32_t index) {return strncmp(s.c_str(), &(*data)[index], length)<0; }
	bool operator()(uint32_t index, string s) {return strncmp(&(*data)[index], s.c_str(), length)<0; }
	bool operator()(uint32_t index, uint32_t index2) { return strncmp(&(*data)[index], &(*data)[index2], length)<0; }

private:
	data_t data;
	int length;
};

void gen_frag(size_t f, string* frag){
	for (size_t i=0; i<frag->size(); ++i){
		size_t mod = f % setting::A;
		frag->at(frag->size()-1-i) = setting::valid_chars[mod];
		f = (f-mod) / setting::A;;
	}
}

indices_t gen_fragpos(string frag){
	indices_t ret = 0;
	const size_t F = frag.size();
	for (size_t i=0; i<F; ++i){
		ret += setting::char2int[frag[F-1-i]] * pow(setting::A,i);
	}
	return ret;
}


int produce_maximal_minimal_unique_string_length(string fasta){
	map<uint32_t,uint64_t> m = maximal_minimal_unique_string_length(fasta);
	ofstream out;
	string dir = fasta.substr(0, fasta.find_last_of("\\/")) + "/";
	umask(0);
	mkdir((dir+"results").c_str(), 0700);
	/*out.open(dir+"results/mmusl", ios::out);
	out << "#" << fasta << endl;
	out << "Maximal_unique_string_length count" << endl;
	for (auto v:m) out << v.first << " " << v.second << endl;
	out.close();*/
	return m.size();
	
}

void existing_neighbors(string fasta, int lf, string frag, string prof_file, float overlap){
	Pointmutation pm = Pointmutation();
	set<string> pm1nei = pm.expand_set(pm.mutate(frag));
	auto data = NamedVector<char>::getNamedVector(fasta, "data");
	auto sorted = NamedVector<uint32_t>::getNamedVector(fasta, "sorted");
	Profile p(lf);
	for (string nei: pm1nei){
		int h = hits(data, sorted, nei);
		if (h > 0) p.insert(nei,h);
	}
	for (string nei: pm1nei){
		if (p.overlap(nei,0) > overlap) cout << nei << endl;
	}
	
	string prof = p.get_profile();
	if (prof_file.size() > 1){
		ofstream out;
		out.open(prof_file, ios::out);
		out << prof;
		out.close();
	}
}

void print_hits_exp_hits(string fasta, string frag){
	auto data = NamedVector<char>::getNamedVector(fasta, "data");
	auto sorted = NamedVector<uint32_t>::getNamedVector(fasta, "sorted");
	auto all_valid_hits = NamedVector<uint32_t>::getNamedVector(fasta, "valid_hits");

	const map<char, double> prob = get_AA_prob_file(NamedVectorBase::cacheFilename(fasta,"_AAcount.map"));
	double exp_hits = frag_prob(&prob, frag)*all_valid_hits->at(frag.size());//data->size();



	int h = hits(data, sorted, frag);
	cout << frag << " " << h << " "<< exp_hits;	
}

int produce_hit_counts(string fasta, int lf){
	map<uint32_t,uint32_t> m = hit_counts_per_protein(fasta, lf);
	ofstream out;
	string dir = fasta.substr(0, fasta.find_last_of("\\/")) + "/";
	umask(0);
	mkdir((dir+"results").c_str(), 0700);
	mkdir((dir+"results/hit_counts").c_str(), 0700);

	out.open((dir+"results/hit_counts/hc."+to_string(lf)), ios::out);
	out << "#" << fasta << " lf=" << lf << endl;
	//vector<string> most = most_often_occ(fasta, lf);
	//for (auto s:most) out << "# " << s << endl;
	out << "hits count" << endl;
	for (auto v:m) {
		out << v.first << " " << v.second << endl;
	}
	out.close();
	return m.size();
}

map<char, double> get_AA_prob_file(string filename){
	ifstream file(filename); // declare file stream: http://www.cplusplus.com/reference/iostream/ifstream/
	string value;
	char c;
	map<char, double> m, m2;
	getline(file, value);
	while(file.good()){
		getline(file, value,','); // read a string until next comma: http://www.cplusplus.com/reference/string/getline/
		if (!file.good()) break;
		c = (value.c_str())[0];
		getline(file, value); // read a string until next comma: http://www.cplusplus.com/reference/string/getline/
		if (!file.good()) break;
		m[c] = stoi(value);
	}

	m.erase('X');
	double sum = 0.0;
	for (size_t i=0; i<setting::all_chars.size(); ++i){
		m2[setting::exchange_chars[i]] += m[setting::all_chars[i]];
		sum += m[setting::all_chars[i]];
	}
	for (auto p:m2)	m2[p.first] = m2[p.first]/sum;
	return m2;
}

map<char, double> get_AA_prob(data_t data){
	map<char,uint32_t> counter;
	for (indices_t i=0; i<data->size(); ++i){
		if (data->at(i)) counter[data->at(i)]++;
	}
	double sum=0; for (auto val: counter) sum += val.second;
	map<char, double> ret;
	for (auto val: counter) ret[val.first] = double(val.second)/sum;
	return ret;
}

double frag_prob(const map<char, double>* prob, string frag){
	double ret = 1.0;
	for (auto ch:frag) ret *= prob->at(ch);
	return ret;
}

// checked 12.12.17 with 3 and 4mers
indices_t valid_hits_of_length(data_t data, size_t lf){
	indices_t ret = 0;
	size_t current_length = 0;
	for (indices_t di=0; di<data->size(); ++di){
		if (data->at(di) == 'X' || data->at(di) == '\0') current_length = 0;
		else current_length ++;
		if (current_length >= lf) ret ++;
	}
	return ret;
}

void word_count_analysis(string fasta, size_t lf, double factor){
	auto data = NamedVector<char>::getNamedVector(fasta, "data");
	auto sorted = NamedVector<uint32_t>::getNamedVector(fasta, "sorted");
	auto diverse_frags = NamedVector<indices_t>::getNamedVector(fasta, "diverse_frags", lf);
	string dir = fasta.substr(0, fasta.find_last_of("\\/")) + "/"; umask(0); mkdir((dir+"results").c_str(), 0700);
	ofstream f_over_under, f_prob_occ;
	f_over_under.open(dir+"results/over_unter_"+to_string(lf)+"mer_" + to_string(factor) + ".csv", ios::out);
	f_prob_occ.open(dir+"results/prob_occ_"+to_string(lf)+"mer.csv", ios::out);
	f_over_under << "# word_count_analysis / " << fasta << " / lf=" << lf << " / factor=" << factor << endl;
	f_over_under << "fragment,hits,exp_hits(data->size()),prob,OU" << endl;
	f_prob_occ << "# word_count_analysis : expected occurrence and observed" << endl;

	const map<char, double> prob = get_AA_prob(data);
	double min_prob = 1; for (auto val:prob) if (val.second < min_prob) min_prob = val.second;
	min_prob = pow(min_prob, lf);
	map<pair<double,uint32_t>, uint32_t> exp2obs;
	indices_t total_hits = valid_hits_of_length(data, lf);
	indices_t hits_without_overreps = total_hits;
	double abs_error = 0.0;
	double sum_exp_hits = 0.0;
	double N = total_hits;//pow(20,lf);
	for (indices_t df=0; df<diverse_frags->size(); ++df){
		string frag = string(&data->at(diverse_frags->at(df))).substr(0,lf);
		int h = hits(data, sorted, frag);
		double f_prob = frag_prob(&prob, frag);
		double exp_hits = f_prob*total_hits;
		sum_exp_hits += exp_hits;
		if (h > exp_hits) abs_error += (h-exp_hits)/N;
		else abs_error -= (h-exp_hits)/N;
		if (factor >= 1){
			if (double(h)*factor < exp_hits){
				//cerr << frag << "," << h << "," << round(exp_hits) << "," << f_prob << ",U" << endl;
				f_over_under << frag << "," << h << "," << round(exp_hits) << "," << f_prob << ",U" << endl;
			} else if(double(h)/factor > exp_hits) {
				//cerr << frag << "," << h << "," << round(exp_hits) << "," << f_prob << ",O" << endl;
				f_over_under << frag << "," << h << "," << round(exp_hits) << "," << f_prob << ",O" << endl;
				hits_without_overreps = hits_without_overreps - h + exp_hits;
			}
		}
		exp2obs[make_pair(round((log2(f_prob)-log2(min_prob))/log2(1.01) ),h)]++;
	}
	cout << "Missing exp_hits = " << total_hits-sum_exp_hits << " / abs_error = " << abs_error << endl;
	cout << "[CHECK] " << hits_without_overreps << " = hits_without_overreps / total_hits = " << total_hits << endl;
	f_prob_occ << "x[min_prob*1.01^x=prob],prob,expected_hits,hits,count" << endl;
	for (auto val:exp2obs) {
		f_prob_occ << val.first.first;
		f_prob_occ << "," << min_prob*pow(1.01,val.first.first);
		f_prob_occ << "," << min_prob*pow(1.01,val.first.first)*total_hits;
		f_prob_occ << "," << val.first.second;
		f_prob_occ << "," << val.second << endl;
	}
	f_prob_occ.close(); f_over_under.close();
}

void all_conditional_entropy(string fasta, size_t hits, bool drop_X){
	auto data = NamedVector<char>::getNamedVector(fasta, "data");
	auto sorted = NamedVector<uint32_t>::getNamedVector(fasta, "sorted");
	auto common_prefix = NamedVector<uint16_t>::getNamedVector(fasta, "prelen");
	auto string_length = NamedVector<uint16_t>::getNamedVector(fasta, "strlen");
	for (size_t lf=1; lf<30; ++lf){
		conditional_entropy(fasta, lf,data, sorted, common_prefix, string_length, hits,hits, drop_X);
	}
}

void subsequent_AAs(data_t data, sorted_t sorted, string fragment, map<char, uint32_t>* entropical_end){
	entropical_end->clear();
	indices_t start = string2sorted(fragment ,data, sorted);
	for (indices_t si=start; ; ++si){
		if (si >= sorted->size()) break;
		if (sorted->at(si)+fragment.size() > data->size()) break;
		if (fragment.compare(string(&data->at(sorted->at(si))).substr(0,fragment.size())) != 0) break;
		char ch = data->at(sorted->at(si)+fragment.size()); // character ch following context c
		if (ch != '\0' && ch != 'X') {
			if (entropical_end->count(ch)) entropical_end->at(ch)++;
			else entropical_end->insert(make_pair(ch,1));
		}
	}
}

char predict_subsequent_AA(data_t data, sorted_t sorted, string fragment, int h){
	map<char, uint32_t> entropical_end;
	subsequent_AAs(data, sorted, fragment, &entropical_end);
	cerr << "Predict follower of : " << fragment << " with #hits=" << h << endl;
	char max = 'X';
	for (auto val:entropical_end) {
		cerr << val.first << " : " << val.second*100.0/float(h) << "%" << endl;
		if (val.second > entropical_end[max]) max = val.first;
	}
	return max;
}






void heuristic_inter_distance_twilight(string fasta, size_t lf, int min_twilight, int max_twilight){
	ofstream f_twilight;
	string dir = fasta.substr(0, fasta.find_last_of("\\/")) + "/";
	umask(0);
	mkdir((dir+"results").c_str(), 0700);
	string out_dir = dir+"results/twilight"+to_string(lf)+"_"+to_string(min_twilight)+"_"+to_string(max_twilight) + "/";
	mkdir((out_dir).c_str(), 0700);
	ofstream f_out;
	auto data = NamedVector<char>::getNamedVector(fasta, "data");
	auto sorted = NamedVector<uint32_t>::getNamedVector(fasta, "sorted");
	auto protstart = NamedVector<uint32_t>::getNamedVector(fasta, "protstart");
	auto protcount = NamedVector<uint32_t>::getNamedVector(fasta, "fasta_prot_count");

	uint64_t M = sorted->size();
	uint64_t m1 = next_greater_prime(uint64_t(sqrt(M))); 
	while (M % m1 == 0) m1 = next_greater_prime(m1);
	uint64_t m2 = next_greater_prime(m1*2);  //19.4.2018 (m1+uint64_t(sqrt(m1)))
	while (M % m2 == 0) m2 = next_greater_prime(m2);

	double random_seq_id = 0.0;
	const map<char, double> prob = get_AA_prob_file(NamedVectorBase::cacheFilename(fasta,"_AAcount.map"));
	for (auto p:prob) random_seq_id += p.second*p.second;
	random_seq_id *= double(lf);
	auto protein_names = NamedVector<uint32_t>::getNamedVector(fasta, "protein_names");
	ifstream fprotein; fprotein.open(NamedVector<uint32_t>::cacheFilename(fasta, "_protein_names.list"));

	cout << "[twilight zone] m1 = " << m1 <<  " m2 = " << m2 << endl;
	cerr << "# m1=" << m1 << " m2 = " << m2 << " lf=" << lf << endl;
	cerr << "iteration deviation" << endl;
	indices_t di1, di2, si1=0, si2=0, h=0;
	for (uint64_t i=0; ; ){
		// pick two random points
		si1 = indices_t((si1+m1) % M); 
		di1 = si1; 
		di1 = sorted->at(si1);
		if (string(&data->at(di1)).size() < lf) continue;
		si2 = indices_t((si2+m2) % M);  
		di2 = si2;
		di2 = sorted->at(si2);
		if (string(&data->at(di2)).size() < lf) continue;
		if (di2>=di1){
			if (di2-di1 < lf) continue;
		} else if (di1-di2 < lf) continue;
		++i;
		//cerr << di1 << " " << di2 << " " << double(di1)-double(di2) << endl;
		int PM = PM_neighbor(&(data->at(di1)), &(data->at(di2)), lf);
		if (PM < 0 || PM>lf) continue;
		if (PM>= min_twilight && PM<= max_twilight){
		    ++h;
			string f1 = string(&(data->at(di1))).substr(0,lf);
			string f2 = string(&(data->at(di2))).substr(0,lf);

			f_out.open(out_dir+to_string(h)+"_1.fasta", ios::out);
			auto p = data2protpos(protstart, di1, protcount);
			string pp = pretty_protein_name(fasta,p.first, &fprotein, protein_names);
			f_out << ">" << pp << " PM=" << PM << endl << f1 << endl;
			f_out.close();

			f_out.open(out_dir+to_string(h)+"_2.fasta", ios::out);
			p = data2protpos(protstart, di2, protcount);
			pp = pretty_protein_name(fasta,p.first, &fprotein, protein_names);
			f_out << ">" << pp << " PM=" << PM << endl << f2 << endl;
			f_out.close();
			/*for (size_t j=0; j<lf; ++j){
				if (data->at(di1+j)==data->at(di2+j)) cout << "O";
				else cout << ".";
			} cout << endl;*/
			//print_protein_with_fragment(fasta, f2, true, false);
		}
	}
}


void cluster_string_max_PM(string fasta, string frag, int max_PM){
	auto data = NamedVector<char>::getNamedVector(fasta, "data");
	auto sorted = NamedVector<uint32_t>::getNamedVector(fasta, "sorted");
	uint32_t si  = string2sorted(frag, data, sorted);
	size_t lf = frag.size();
	string frag_friend = ""; int pm = lf;
	if (si >= sorted->size()){
		for (indices_t di=0; di<sorted->size(); ++di){
			int PM = PM_neighbor(frag, &(data->at(di)), lf);
			if (PM >= 0){
				if (PM < pm){
					pm = PM; frag_friend=string(&(data->at(di))).substr(0,lf);
				}
			}
		}
		si  = string2sorted(frag_friend, data, sorted);
		cout << " fragment " << frag << " no identical match, take closest neighbor:" << endl;
		cout << " neighbor " << frag_friend << " PM=" <<  pm << endl;
	}


	indices_t di = sorted->at(si);
	connected_component_max_PM(fasta, di, di, frag.size(), max_PM);
}

void cluster_max_PM(string fasta, size_t lf, int max_PM){
	auto data = NamedVector<char>::getNamedVector(fasta, "data");
	auto sorted = NamedVector<uint32_t>::getNamedVector(fasta, "sorted");
	for (indices_t di = 0; ; di += 10003){
		di = di % data->size();
		string frag = string(&data->at(di)).substr(0,lf);
		if (frag.size() != lf) continue;
		connected_component_max_PM(fasta, di, di, frag.size(), max_PM);
	}
}


void connected_component_max_PM(string fasta, indices_t di1, indices_t di2, size_t lf, int max_PM){
	auto data = NamedVector<char>::getNamedVector(fasta, "data");
	auto sorted = NamedVector<uint32_t>::getNamedVector(fasta, "sorted");
	set<indices_t> cc1,cc2, cc;
	set<indices_t>* p_cc1 = &cc1;
	set<indices_t>* p_cc2 = &cc2;
	p_cc1->insert(di1); cc.insert(di1);

	int round = 0;
	while (p_cc1->size() > 0){
		++ round; cerr << "[ROUND] " << round << endl;
		for (auto comp_di:*p_cc1){
			cerr << "[COMP]" << string(&(data->at(comp_di))).substr(0,lf) << " "<<  comp_di << endl;
			for (indices_t di=0; di<sorted->size(); ++di){
				if (di == comp_di) { di+=(lf-1); continue; }
				int PM = PM_neighbor(&(data->at(comp_di)), &(data->at(di)), lf);
				if (PM >= 0 && PM <= max_PM){
					if (di == di2) cerr << "[DEEP HOMOLOGY] round=" << round <<  " currently found homologs:" << cc.size()+p_cc1->size()+p_cc2->size() << endl;
					if (p_cc1->find(di) == p_cc1->end() && cc.find(di) == cc.end()){ 
						p_cc2->insert(di);
						cc.insert(di);
						cerr <<  string(&(data->at(di))).substr(0,lf) <<  " " << di << " " << PM << endl;
					} //else cout << "back-find" << endl;
				}
			}
			map<int,int> distance_counter;
			cerr << "cc.size=" << cc.size() << endl;
			for (auto val1:cc){
				for (auto val2:cc){
					if (val1>=val2) continue;
					distance_counter[PM_neighbor(&(data->at(val1)), &(data->at(val2)), lf)]++;
				}
			}
			for (auto val:distance_counter) cerr << val.first << " " << val.second << endl;
		}

		p_cc1->clear();
		swap(p_cc1, p_cc2);
	}

	map<int,int> distance_counter;
	cout << "# cc.size=" << cc.size() <<endl;
	cout << "#" << string(&(data->at(di1))).substr(0,lf) << endl;
	if (di1 != di2 ){
		cout << "#" << string(&(data->at(di2))).substr(0,lf) << endl;
		if (cc.find(di2) != cc.end()) cout << "#[DEEP HOMOLOGY]" << endl;
	}
	for (auto val1:cc){
		for (auto val2:cc){
			if (val1>=val2) continue;
			distance_counter[PM_neighbor(&(data->at(val1)), &(data->at(val2)), lf)]++;
		}
	}
	for (auto val:distance_counter) cout << val.first << " " << val.second << endl;

}

int compositional_overlap(indices_t di1, indices_t di2, data_t data, size_t lf){
	map<char, size_t> comp1;
	map<char, size_t> comp2;
	for (indices_t i=0; i<lf; ++i){
		comp1[data->at(di1+i)]++;
		comp2[data->at(di2+i)]++;
	}
	int overlap = 0;
	for (auto val:comp1){
		overlap += min(val.second, comp2[val.first]);
	}
	return overlap;
}

// TODO
//void frag_shuf_para_comp(string fasta, size_t lf, int start, int org, int locality_PM){
void frag_shuf_para_comp(string fasta, size_t lf, int org){
	ofstream f_frag_shuf;
	string dir = fasta.substr(0, fasta.find_last_of("\\/")) + "/";
	umask(0); 
	mkdir((dir+"results").c_str(), 0700);
	mkdir((dir+"results/inter_dist").c_str(), 0700);
	mkdir((dir+"results/frag_pairs").c_str(), 0700);
	mkdir((dir+"results/locality").c_str(), 0700);
	string frag_pair = dir+"results/frag_pairs/" +to_string(lf)+"mer_para"+ to_string(org) +"_frag_pairs.csv";
	string frag_pair_prot = dir+"results/frag_pairs/" +to_string(lf)+"mer_prot"+ to_string(org) +"_frag_pairs.csv";
	f_frag_shuf.open(frag_pair, ios::out); f_frag_shuf.close();
	f_frag_shuf.open(frag_pair_prot, ios::out); f_frag_shuf.close();

	auto data = NamedVector<char>::getNamedVector(fasta, "data");
	auto protstart = NamedVector<uint32_t>::getNamedVector(fasta, "protstart");
	auto protcount = NamedVector<uint32_t>::getNamedVector(fasta, "fasta_prot_count");
	auto orgstart = NamedVector<uint32_t>::getNamedVector(fasta, "orgstart");

	uint64_t range = (orgstart->at(org+1)-orgstart->at(org)-1);
	const uint64_t org_start = orgstart->at(org);
	cout << " # orgstart = " << org_start << "  til " << org_start + range << endl; 
	uint64_t M = range;
	uint64_t m1 = next_greater_prime(uint64_t(sqrt(M)));  while (M % m1 == 0) m1 = next_greater_prime(m1);
	uint64_t m2 = next_greater_prime(m1*2); while (M % m2 == 0) m2 = next_greater_prime(m2);

	indices_t di1 = org_start;
	indices_t di2 = org_start+1;

	vector<map<size_t, uint64_t>> storage; // nat, rand, prot_nat, prot_rand 
	for (size_t i=0; i<4; ++i) storage.push_back(map<size_t, uint64_t>());
	vector<size_t> frag_pair_count(lf+1, 0);
	vector<size_t> frag_pair_count_prot(lf+1, 0);
	size_t max_frag_pair_count = 500;
	uint64_t round = 0;
	uint64_t N = 0;
	
	size_t transition = cal_transition(lf, 0.0000001, 0.0636);
	pair<size_t, size_t> twi_zone = extend_to_twi_zone(lf, transition, lf*0.2);
	cout << "Twilight zone: " << twi_zone.first << " " << twi_zone.second << endl;
	map<size_t,map<string,uint64_t>> locality;
	map<size_t,map<string,uint64_t>> locality_rand;
	map<size_t,map<string,uint64_t>> locality_comp;
	map<size_t,map<string,uint64_t>> locality_prot;
	map<size_t,map<string,uint64_t>> locality_rand_prot;
	map<size_t,map<string,uint64_t>> locality_comp_prot;
	bool same_prot;

	for (uint64_t i=0; i<pow(10.0, 8.0); ){
		++N;
		// pick two random points
		di1 = org_start + indices_t((di1-org_start+m1) % M);
		di2 = org_start + indices_t((di2-org_start+m2) % M);
		if (di2<=di1) continue;

		if (di1 == org_start) di2 = org_start + indices_t((di2-org_start+m2) % M);

		// Ultimate end
		if (di1 == org_start && di2 == org_start+1){
			cerr << "[INFO] made a round, valid pairs=" << i << " seen pairs=" << N << " range=" << range << ", possible pairs=" << range*range << endl;
			exit(0);
		}

		// not long enough
		if (string(&data->at(di1)).size() < lf) continue;
		if (string(&data->at(di2)).size() < lf) continue;

		// different organism
		pair<uint32_t, uint32_t> p1 = data2protpos(protstart, di1, protcount);
		pair<uint32_t, uint32_t> p2 = data2protpos(protstart, di2, protcount);
		if (protpos2prot(p1) == protpos2prot(p2)) same_prot = true;
		else same_prot = false;
		if (di1-di2 < lf) continue; //overlapping fragments

		int PM = PM_neighbor(&(data->at(di1)), &(data->at(di2)), lf);
		if (PM < 0 || PM>lf) continue;

		// ###### valid fragment pair ###########################

		++i;
		if (same_prot) storage[2][PM]++;
		else storage[0][PM]++;

		string str1 = string(&(data->at(di1))).substr(0,lf);
		string str2 = string(&(data->at(di2))).substr(0,lf);


		if (PM >= twi_zone.first && PM<= twi_zone.second && PM<=lf-2){
			size_t d=0; char ch = '\0';
			vector<size_t> loc_pos;
			for (size_t ii=0; ii<lf; ++ii) if (data->at(di1+ii) == data->at(di2+ii)) loc_pos.push_back(ii);
			for (int jj=0; jj<loc_pos.size(); jj++){
				for (int jjj=jj+1; jjj<loc_pos.size(); jjj++){
					if (same_prot) locality_prot[(PM<<8) + (loc_pos[jjj]-loc_pos[jj])][string(1,data->at(di1+loc_pos[jj]))+string(1,data->at(di1+loc_pos[jjj]))]++;
					else locality[(PM<<8) + (loc_pos[jjj]-loc_pos[jj])][string(1,data->at(di1+loc_pos[jj]))+string(1,data->at(di1+loc_pos[jjj]))]++;
				}
			}
		}

		if (frag_pair_count[PM] < max_frag_pair_count) {
			frag_pair_count[PM]++;
			f_frag_shuf.open(frag_pair, ios::app);      f_frag_shuf << lf << " " << PM << " " << int(PM)-int(transition) << " " << str1 << " " << str2 << endl; f_frag_shuf.close();
		}
		if (frag_pair_count_prot[PM] < max_frag_pair_count) {
			frag_pair_count[PM]++;
			f_frag_shuf.open(frag_pair_prot, ios::app); f_frag_shuf << lf << " " << PM << " " << int(PM)-int(transition) << " " << str1 << " " << str2 << endl; f_frag_shuf.close();
		}

		size_t PMr = 0;
		for (size_t ttt=0; ttt<10; ++ttt){
			random_shuffle(str1.begin(), str1.end());
			random_shuffle(str2.begin(), str2.end());
			PMr = 0; for (size_t jj=0; jj<lf; ++jj) if (str1[jj] != str2[jj]) ++PMr;
			if (same_prot) storage[3][PMr]++;
			else storage[1][PMr]++;
			size_t d=0; char ch = '\0';
			vector<size_t> loc_pos;
			for (size_t ii=0; ii<lf; ++ii) if (str1[ii] == str2[ii]) loc_pos.push_back(ii);
			for (int jj=0; jj<loc_pos.size(); jj++){
				for (int jjj=jj+1; jjj<loc_pos.size(); jjj++){
					if (PM >= twi_zone.first && PM<= twi_zone.second && PM<=lf-2){
						if (same_prot) locality_comp_prot[(PM<<8) + (loc_pos[jjj]-loc_pos[jj])][string(1, str1[loc_pos[jj]])+string(1, str2[loc_pos[jjj]])]++;
						else           locality_comp[(PM<<8) + (loc_pos[jjj]-loc_pos[jj])][string(1, str1[loc_pos[jj]])+string(1, str2[loc_pos[jjj]])]++;
					} else {
						if (same_prot) locality_comp_prot[(PM<<8)][string(1, str1[loc_pos[jj]])+string(1, str2[loc_pos[jjj]])]++;
						else           locality_comp[(PM<<8)][string(1, str1[loc_pos[jj]])+string(1, str2[loc_pos[jjj]])]++;
					}

					if (PMr >= twi_zone.first && PMr<= twi_zone.second && PMr<=lf-2 && same_prot){
						locality_rand_prot[(PMr<<8) + (loc_pos[jjj]-loc_pos[jj])][string(1, str1[loc_pos[jj]])+string(1, str2[loc_pos[jjj]])]++;
					}
					if (PMr >= twi_zone.first && PMr<= twi_zone.second && PMr<=lf-2 && !(same_prot)){
						locality_rand[(PMr<<8) + (loc_pos[jjj]-loc_pos[jj])][string(1, str1[loc_pos[jj]])+string(1, str2[loc_pos[jjj]])]++;
					}
				}
			}
		}

		// evaluate!
		if (i % 1000000 == 0){
			cerr << "eval! " << i << " N=" << N << endl;
			string name = dir+"results/inter_dist/" +to_string(lf)+"mer_para_"+ to_string(org) +".csv";
			f_frag_shuf.open(name, ios::out);	
			f_frag_shuf << "# Replacement for binominal model, whole fragment shuffling" << endl;
			f_frag_shuf << "lf PM nat_count rand_count prot_count prot_rand_count" << endl;
			for (size_t pm=0; pm<=lf; ++pm) {
				f_frag_shuf << lf << " " << pm << " " << storage[0][pm] << " "  << storage[1][pm] << " "  << storage[2][pm] << " "  << storage[3][pm] << endl;
			}
			f_frag_shuf.close();
			// FILE: Locality
			name = dir+"results/locality/" +to_string(lf)+"mer_para_"+ to_string(org) +"_locality.csv";
			f_frag_shuf.open(name, ios::out);	
			f_frag_shuf << "# AA Patterns in twilight zone " << twi_zone.first << "-" << twi_zone.second << endl;
			f_frag_shuf << "# transition at " << transition << endl;
			f_frag_shuf << "lf PM twi_offset d AA1 AA2 nat_count rand_count shuffle_count prot_count prot_rand_count prot_shuffle_count" << endl;
			for (int pm=twi_zone.first; pm<=twi_zone.second; ++pm){
				for (uint32_t d=1;d<lf;++d){
					for (auto AA1: setting::valid_chars){
						for (auto AA2: setting::valid_chars){
							string key = string(1, AA1)+string(1, AA2);
							f_frag_shuf << lf << " " << pm << " " << int(pm)-int(transition) << " " << d << " "<< AA1 << " " << AA2 << " ";
							f_frag_shuf << locality[(pm<<8)+d][key] <<" " << locality_rand[(pm<<8)+d][key] << " " << locality_comp[(pm<<8)+d][key] << " ";
							f_frag_shuf << locality_prot[(pm<<8)+d][key] <<" " << locality_rand_prot[(pm<<8)+d][key] << " " << locality_comp_prot[(pm<<8)+d][key] << endl;
						}
					}
				}
			}
			f_frag_shuf.close();
		}
	}
}



//XXX OLD
void frag_shuf_para_comp(string fasta, size_t lf, int start, int org, int locality_PM){
	ofstream f_frag_shuf;
	string dir = fasta.substr(0, fasta.find_last_of("\\/")) + "/";
	umask(0); 
	mkdir((dir+"results").c_str(), 0700);
	mkdir((dir+"results/frag_shuf_op_comp").c_str(), 0700);
	auto data = NamedVector<char>::getNamedVector(fasta, "data");
	auto protstart = NamedVector<uint32_t>::getNamedVector(fasta, "protstart");
	auto protcount = NamedVector<uint32_t>::getNamedVector(fasta, "fasta_prot_count");
	auto orgstart = NamedVector<uint32_t>::getNamedVector(fasta, "orgstart");

	uint64_t range = (orgstart->at(org+1)-orgstart->at(org)-1);
	const uint64_t org_start = orgstart->at(org);
	cout << " # orgstart = " << org_start << "  til " << org_start + range << endl; 

	vector<map<string,uint64_t>> locality;
	vector<map<string,uint64_t>> locality_rand;
	for (size_t ii=0; ii<=lf; ++ii) locality.push_back(map<string,uint64_t>());
	for (size_t ii=0; ii<=lf; ++ii) locality_rand.push_back(map<string,uint64_t>());


	uint64_t M = range;
	uint64_t m1 = next_greater_prime(uint64_t(sqrt(M))); 
	while (M % m1 == 0) m1 = next_greater_prime(m1);
	uint64_t m2 = next_greater_prime(m1*2);  
	while (M % m2 == 0) m2 = next_greater_prime(m2);

	indices_t di1 = org_start;
	indices_t di2 = org_start+1;

	cerr << "START IS INVALID!" << endl;
	cerr << "m1 = " << m1 << " m2 = " << m2 << "  M = " << M << endl;
	cerr << "di1 = " << di1 << " di2 = " << di2 << endl;

	
	vector<vector<map<size_t, uint64_t>>> storage; // 10compo -> 3 choices : nat,rand,same_protein
	vector<map<size_t, uint64_t>> choice; for (int i=0; i<3; i++) choice.push_back(map<size_t, uint64_t>());
	for (int i=0; i<10; i++) storage.push_back(choice);

	uint64_t N =0;
	bool same_prot = 0;
	for (uint64_t i=0; i<range*(range-1)/2; ){
		++N;
		di1 = org_start + indices_t((di1-org_start+m1) % M);
		di2 = org_start + indices_t((di2-org_start+m2) % M);

		if (di1 == org_start) {
			di2 = org_start + indices_t((di2-org_start+m2) % M);
			cerr << "semi round " << di1-org_start << " " << di2-org_start << endl;
		}

		// Ultimate end
		if (di1 == org_start && di2 == org_start+1){
			cerr << "[INFO] made a round, valid pairs=" << i << " seen pairs=" << N << " range=" << range << ", possible pairs=" << range*range << endl;
			exit(0);
		}

		if (di2<=di1) continue;

		//cerr << "di1=" << di1-org_start << " di2=" << di2-org_start << " diff=" << int(di1)-int(di2) << endl;

		// not long enough
		if (string(&data->at(di1)).size() < lf) continue;
		if (string(&data->at(di2)).size() < lf) continue;

		pair<uint32_t, uint32_t> p1 = data2protpos(protstart, di1, protcount);
		pair<uint32_t, uint32_t> p2 = data2protpos(protstart, di2, protcount);
		if (protpos2organism(p1) != protpos2organism(p2)) {
			cerr << "[ERROR] protpos2organism(p1) != protpos2organism(p2)" << endl;
			cerr << i << "=i " << di1 << " " << di2 << endl;
			cerr << " org1=" << protpos2organism(p1) << " org2=" << protpos2organism(p2) << endl;
			exit(0);
		}
		if (protpos2prot(p1) == protpos2prot(p2)) same_prot = true;
		else same_prot = false;

		//overlapping fragments
		if (di2>=di1){
			if (di2-di1 < lf) continue;
		} else if (di1-di2 < lf) continue;


		// compositional overlap
		int compo = compositional_overlap(di1,di2,data,lf);
		//if (double(compo)/double(lf) <= 0.5) continue;
		int compo_index = min(floor((double(compo)/double(lf)*10.0)), 9.0);

		int PM = PM_neighbor(&(data->at(di1)), &(data->at(di2)), lf);
		if (PM < 0 || PM>lf) continue;
		// ********** valid fragment pair ****************
		++i;
		storage[compo_index][0][PM]++;
		if (same_prot) storage[compo_index][2][PM]++;

		if (PM == locality_PM){
			size_t d=0; char ch = '\0';
			vector<size_t> loc_pos;
			for (size_t ii=0; ii<lf; ++ii) if (data->at(di1+ii) == data->at(di2+ii)) loc_pos.push_back(ii);
			for (int jj=0; jj<loc_pos.size(); jj++){
				for (int jjj=jj+1; jjj<loc_pos.size(); jjj++){
					locality[jjj-jj][string(1,data->at(di1+jj))+string(1,data->at(di1+jjj))]++;
				}
			}
		}

		string str1 = string(&(data->at(di1))).substr(0,lf);
		string str2 = string(&(data->at(di2))).substr(0,lf);

		random_shuffle(str1.begin(), str1.end());
		random_shuffle(str2.begin(), str2.end());

		PM = 0;
		for (size_t ii=0; ii<lf; ++ii) if (str1[ii] != str2[ii]) ++PM;
		// XXX store shuffle
		storage[compo_index][1][PM]++;

		// evaluate!
		if (i % 10000 == 0){
			cerr << "eval! " << i << " N=" << N << " / " << range*range << endl;
			string name = dir+"results/frag_shuf_op_comp/" +to_string(lf)+"mer_"+ to_string(start) +"_para_"+to_string(org)+".csv";
			f_frag_shuf.open(name, ios::out);	
			f_frag_shuf << "# Replacement for binominal model, whole fragment shuffling" << endl;
			f_frag_shuf << "# compo are 10 bins 0:0- 1:10- 2:20- 3:30- 4:40- 5:50-60% 6:60-70% 7:70-80% 8:80-90% 9:90-100%" << endl;
			f_frag_shuf << "lf compo PM nat_count rand_count same_prot_count" << endl;
			for (size_t c=0; c<10; ++c){
				for (size_t pm=0; pm<=lf; ++pm) {
					f_frag_shuf << lf << " " << c << " " << pm << " ";
					f_frag_shuf << storage[c][0][pm] << " "  << storage[c][1][pm] <<  " "  << storage[c][2][pm] << endl;
				}
			}
			f_frag_shuf.close();

			name = dir+"results/frag_shuf_op_comp/" +to_string(lf)+"mer_"+ to_string(start) +"_para_"+to_string(org)+"_locality_"+to_string(locality_PM)+"PM.csv";
			f_frag_shuf.open(name, ios::out);	
			f_frag_shuf << "# Locality of lf=" << to_string(lf) << " PM=" << to_string(locality_PM) << endl;
			f_frag_shuf << "distance A1 A2 nat_count rand_count" << endl;
			for (size_t ii=0; ii<locality.size(); ++ii){
					for (auto val:locality[ii]){
						f_frag_shuf << ii << " " << val.first[0] << " " << val.first[1] << " " << val.second << " "<< locality_rand[ii][val.first] << endl;
					}
			}
			f_frag_shuf.close();

		}

	}
	cerr << "m1 = " << m1 << " m2 = " << m2 << "  M = " << M << endl;
}

size_t cal_transition(size_t lf, double precision, double match){
	for (size_t pm=lf; pm>=0; --pm){
		double exp = pow(1.0-match, pm) * pow(match, lf-pm) * binom(lf, pm);
		cout << "pm=" << pm << " exp=" << exp << endl;
		if (exp <= precision) return pm;
	}
	return 0;
}

pair<size_t, size_t> extend_to_twi_zone(size_t lf, size_t transition, size_t radius){
	size_t min = (transition>radius)?transition-radius:0;
	size_t max = (transition+radius<lf)?transition+radius:lf;
	return(make_pair(min, max));
}


// ORTHOGONAL
void frag_shuf_ortho_comp(string fasta, size_t lf, int start){
	ofstream f_frag_shuf;
	string dir = fasta.substr(0, fasta.find_last_of("\\/")) + "/";
	umask(0); 
	mkdir((dir+"results").c_str(), 0700);
	mkdir((dir+"results/inter_dist").c_str(), 0700);
	mkdir((dir+"results/frag_pairs").c_str(), 0700);
	mkdir((dir+"results/locality").c_str(), 0700);
	string frag_pair = dir+"results/frag_pairs/" +to_string(lf)+"mer_"+ to_string(start) +"_frag_pairs.csv";
	f_frag_shuf.open(frag_pair, ios::out); f_frag_shuf.close();

	auto data = NamedVector<char>::getNamedVector(fasta, "data");
	auto protstart = NamedVector<uint32_t>::getNamedVector(fasta, "protstart");
	auto protcount = NamedVector<uint32_t>::getNamedVector(fasta, "fasta_prot_count");

	uint64_t M = data->size();
	uint64_t m1 = next_greater_prime(uint64_t(sqrt(M))); 
	while (M % m1 == 0) m1 = next_greater_prime(m1);
	uint64_t m2 = next_greater_prime(m1*2);  //19.4.2018 (m1+uint64_t(sqrt(m1)))
	while (M % m2 == 0) m2 = next_greater_prime(m2);

	indices_t di1 = 0;
	indices_t di2 = (M*start*m2/1000) % M; // 1000 valid start values!

	vector<map<size_t, uint64_t>> storage; // nat, rand, nat_ali, rand_ali
	for (size_t i=0; i<2; ++i) storage.push_back(map<size_t, uint64_t>());
	vector<size_t> frag_pair_count(lf+1, 0);
	size_t max_frag_pair_count = 100;
	uint64_t round = 0;
	uint64_t N = 0;

	map<size_t,uint64_t> para_hits;
	size_t transition = cal_transition(lf, 0.0000001, 0.0636);
	pair<size_t, size_t> twi_zone = extend_to_twi_zone(lf, transition, lf*0.2);

	cout << "Twilight zone: " << twi_zone.first << " " << twi_zone.second << " transition = " << transition << endl;
	map<size_t,map<string,uint64_t>> locality;
	map<size_t,map<string,uint64_t>> locality_rand;
	map<size_t,map<string,uint64_t>> locality_comp;
	//map<size_t,map<size_t,<map<string,uint64_t>>>> locality_rand;

	for (uint64_t i=0; i<pow(10.0,8.0); ){
		++N;
		// pick two random points
		di1 = indices_t((di1+m1) % M);
		di2 = indices_t((di2+m2) % M);  
		if (di2<=di1) continue;

		if (di1 == 0) { // semi-round
			di2 = indices_t((di2+m2) % M);
			cerr << "semi round " << di1 << " " << di2 << endl;
		}
		if (di1 == 0 && di2 == (M*(start+1)*m2/1000) % M){ // Ultimate end
			cerr << "[INFO] made a round, valid pairs=" << i << " seen pairs=" << N << " range=" << M << ", possible pairs=" << double(M*M) << endl;
			exit(0);
		}

		// not long enough
		if (string(&data->at(di1)).size() < lf) continue;
		if (string(&data->at(di2)).size() < lf) continue;

		int PM = PM_neighbor(&(data->at(di1)), &(data->at(di2)), lf);
		if (PM < 0 || PM>lf) continue;

		// ###### valid fragment pair ###########################

		// different organism
		pair<uint32_t, uint32_t> p1 = data2protpos(protstart, di1, protcount);
		pair<uint32_t, uint32_t> p2 = data2protpos(protstart, di2, protcount);
		if (protpos2organism(p1) == protpos2organism(p2)){
			para_hits[PM]++;
			continue;
		}


		++i;
		storage[0][PM]++;

		string str1 = string(&(data->at(di1))).substr(0,lf);
		string str2 = string(&(data->at(di2))).substr(0,lf);


		if (PM >= twi_zone.first && PM<= twi_zone.second && PM<=lf-2){
			size_t d=0; char ch = '\0';
			vector<size_t> loc_pos;
			for (size_t ii=0; ii<lf; ++ii) if (data->at(di1+ii) == data->at(di2+ii)) loc_pos.push_back(ii);
			for (int jj=0; jj<loc_pos.size(); jj++){
				for (int jjj=jj+1; jjj<loc_pos.size(); jjj++){
					locality[(PM<<8) + (loc_pos[jjj]-loc_pos[jj])][string(1,data->at(di1+loc_pos[jj]))+string(1,data->at(di1+loc_pos[jjj]))]++;
				}
			}
		}

		if (frag_pair_count[PM] < max_frag_pair_count) {
			frag_pair_count[PM]++;
			f_frag_shuf.open(frag_pair, ios::app);
			f_frag_shuf << lf << " " << PM << " " << int(PM)-int(transition) << " " << str1 << " " << str2 << endl;
			f_frag_shuf.close();
		}

		size_t PMr = 0;
		for (size_t ttt=0; ttt<10; ++ttt){
			random_shuffle(str1.begin(), str1.end());
			random_shuffle(str2.begin(), str2.end());
			PMr = 0; for (size_t jj=0; jj<lf; ++jj) if (str1[jj] != str2[jj]) ++PMr;
			storage[1][PMr]++;
			size_t d=0; char ch = '\0';
			vector<size_t> loc_pos;
			for (size_t ii=0; ii<lf; ++ii) if (str1[ii] == str2[ii]) loc_pos.push_back(ii);
			for (int jj=0; jj<loc_pos.size(); jj++){
				for (int jjj=jj+1; jjj<loc_pos.size(); jjj++){
					if (PM >= twi_zone.first && PM<= twi_zone.second && PM<=lf-2){
						locality_comp[(PM<<8) + (loc_pos[jjj]-loc_pos[jj])][string(1, str1[loc_pos[jj]])+string(1, str2[loc_pos[jjj]])]++;
					} else {
						locality_comp[(PM<<8)][string(1, str1[loc_pos[jj]])+string(1, str2[loc_pos[jjj]])]++;
					}

					if (PMr >= twi_zone.first && PMr<= twi_zone.second && PMr<=lf-2){
						locality_rand[(PMr<<8) + (loc_pos[jjj]-loc_pos[jj])][string(1, str1[loc_pos[jj]])+string(1, str2[loc_pos[jjj]])]++;
					} else locality_rand[(PMr<<8)][string(1, str1[loc_pos[jj]])+string(1, str2[loc_pos[jjj]])]++;
				}
			}
		}

		// evaluate!
		if (i % 1000000 == 0){
			cerr << "eval! " << i << " N=" << N << endl;
			string name = dir+"results/inter_dist/" +to_string(lf)+"mer_"+ to_string(start) +"_ortho.csv";
			f_frag_shuf.open(name, ios::out);	
			f_frag_shuf << "# Replacement for binominal model, whole fragment shuffling" << endl;
			f_frag_shuf << "lf PM nat_count para_count rand_count" << endl;
			for (size_t pm=0; pm<=lf; ++pm) {
				f_frag_shuf << lf << " " << pm << " " << storage[0][pm] << " "  << para_hits[pm] << " " << storage[1][pm] << endl;
			}
			f_frag_shuf.close();
			// FILE: Locality
			name = dir+"results/locality/" +to_string(lf)+"mer_"+ to_string(start) +"_locality.csv";
			f_frag_shuf.open(name, ios::out);	
			f_frag_shuf << "# AA Patterns in twilight zone " << twi_zone.first << "-" << twi_zone.second << endl;
			f_frag_shuf << "# transition at " << transition << endl;
			f_frag_shuf << "lf PM twi_offset d AA1 AA2 nat_count rand_count shuffle_count" << endl;
			for (int pm=twi_zone.first; pm<=twi_zone.second; ++pm){
				for (uint32_t d=1;d<lf;++d){
					for (auto AA1: setting::valid_chars){
						for (auto AA2: setting::valid_chars){
							string key = string(1, AA1)+string(1, AA2);
							f_frag_shuf << lf << " " << pm << " " << int(pm)-int(transition) << " "<< d << " "<< AA1 << " " << AA2 << " ";
							f_frag_shuf << locality[(pm<<8)+d][key] <<" " << locality_rand[(pm<<8)+d][key] << " " << locality_comp[(pm<<8)+d][key] << endl;
						}
					}
				}
			}
			f_frag_shuf.close();
		}
	}
}


















double compo_exp(string str1, string str2){
	map<char, size_t> compo1;
	map<char, size_t> compo2;
	for (int i=0; i<str1.size(); ++i) compo1[str1[i]]++;
	for (int i=0; i<str2.size(); ++i) compo2[str2[i]]++;
	double ret = 0.0;
	for (auto val: compo1){
		ret += compo1[val.first] * compo2[val.first];
	}
	return ret/double(str1.size() + str2.size() / 2.0);
}




void frag_shuf_model(string fasta, size_t max_lf, uint64_t start){
	ofstream f_frag_shuf;
	string dir = fasta.substr(0, fasta.find_last_of("\\/")) + "/";
	umask(0);
	mkdir((dir+"results").c_str(), 0700);
	auto data = NamedVector<char>::getNamedVector(fasta, "data");
	auto sorted = NamedVector<uint32_t>::getNamedVector(fasta, "sorted");

	uint64_t M = sorted->size();
	uint64_t m1 = next_greater_prime(uint64_t(sqrt(M))); 
	while (M % m1 == 0) m1 = next_greater_prime(m1);
	uint64_t m2 = next_greater_prime(m1+uint64_t(sqrt(m1)));  //19.4.2018 (m1+uint64_t(sqrt(m1)))
	while (M % m2 == 0) m2 = next_greater_prime(m2);

	vector<map<int, uint64_t>> context_shuf(max_lf+1);
	vector<map<int, uint64_t>> frag_shuf(max_lf+1);
	vector<map<int, uint64_t>> nat_frag(max_lf+1);

	cerr << "[frag_shuf_model] m1 = " << m1 <<  " m2 = " << m2 << endl;
	cerr << "# m1=" << m1 << " m2 = " << m2 << " max_lf=" << max_lf << endl;
	cerr << "iteration deviation" << endl;
	indices_t di1, di2;
	indices_t si1 = indices_t(start*sqrt(m1)) % M;
	indices_t si2 = indices_t(2*start*sqrt(m1)) % M;
	cerr << "si1 = " << si1 << " si2 = " << si2 << endl;
	for (uint64_t i=0; i<sorted->size()*(sorted->size()-1)/2; ){
		// pick two random points
		si1 = indices_t((si1+m1) % M); 
		di1 = si1;

		si2 = indices_t((si2+m2) % M);  
		di2 = si2;

		for (size_t lf=1; lf<=max_lf; ++lf){
			// not long enough
			if (string(&data->at(di1)).size() < lf) break;
			if (string(&data->at(di2)).size() < lf) break;

			//overlapping fragments
			/*if (di2>=di1){
				if (di2-di1 < lf) break;
			} else if (di1-di2 < lf) break;
			*/

			// valid fragment pair
			int PM = PM_neighbor(&(data->at(di1)), &(data->at(di2)), lf);
			nat_frag[lf][PM]++;
			if (PM < 0 || PM>lf) break;
			string str1 = string(&(data->at(di1))).substr(0,lf);
			string str2 = string(&(data->at(di2))).substr(0,lf);

			auto it2 = str1.begin(); 
			auto itt2 = str2.begin(); 

			advance(it2, lf/2);
			random_shuffle(str1.begin(), it2);
			random_shuffle(it2, str1.end());

			advance(itt2, lf/2);
			random_shuffle(str2.begin(), itt2);
			random_shuffle(itt2, str2.end());

			PM = 0;
			for (size_t i=0; i<lf; ++i) if (str1[i] != str2[i]) ++PM;
			context_shuf[lf][PM]++;



			// all shuf
			random_shuffle(str1.begin(), str1.end());random_shuffle(str2.begin(), str2.end());
			PM = 0;
			for (size_t i=0; i<lf; ++i) if (str1[i] != str2[i]) ++PM;
			frag_shuf[lf][PM]++;
		}
		++i;
		// evaluate!
		if (i % 1000000 == 0){
			cerr << "eval! " << i << endl;
		    f_frag_shuf.open(dir+"results/intra_dis_demi_frag_shuf_" + to_string(start) +".csv", ios::out);	
			f_frag_shuf << "# Replacement for binominal model, whole fragment shuffling context= half of fragment shuffle" << endl;
			f_frag_shuf << "# ALLOW OVERLAP" << endl;
			f_frag_shuf << "lf PM exp_count exp_fraction context_count context_fraction nat_count nat_fraction" << endl;
			for (size_t lf=1; lf<=max_lf; ++lf){
				double s = 0.0; for (auto val:frag_shuf[lf]) s += val.second;
				for (size_t pm=0; pm<=lf; ++pm) {
					f_frag_shuf << lf << " " << pm << " ";
					f_frag_shuf << frag_shuf[lf][pm]<< " " << double(frag_shuf[lf][pm])/s << " ";
					f_frag_shuf << context_shuf[lf][pm]<< " " << double(context_shuf[lf][pm])/s << " ";
					f_frag_shuf << nat_frag[lf][pm] << " " << nat_frag[lf][pm]/s << endl;
				}
			}
			f_frag_shuf.close();
		}
	}
}



//#############
void heuristic_inter_distance(string fasta, size_t lf, size_t rep, size_t precision){
	if (rep<2){
		cerr << "[ERROR] in heuristic_inter_distance: please enter at least 2 replica" << endl;
		return;
	}
	ofstream f_intra_dis, f_intra_dis_pos, f_very_similar, f_exchange, f_2PM_distance, f_frag_shuf;
	string dir = fasta.substr(0, fasta.find_last_of("\\/")) + "/";
	umask(0);
	mkdir((dir+"results").c_str(), 0700);
	f_very_similar.open(dir+"results/intra_dis_" + to_string(lf) + "mer_similars.txt", ios::out);
	f_very_similar << "di1 di2 PM string1 string2" << endl;
	auto data = NamedVector<char>::getNamedVector(fasta, "data");
	auto sorted = NamedVector<uint32_t>::getNamedVector(fasta, "sorted");
	map<size_t, uint64_t> PM_distances;
	vector<map<int,uint64_t>> distances(rep);
	vector<map<int,double>> dis_prob(rep);
	vector<map<uint16_t,uint64_t>> exchanges; // PM-> (char,char) -> count
	for (int i=0; i<=lf; ++i) exchanges.push_back(map<uint16_t,uint64_t>());

	vector<uint64_t> PM_pos(lf);
	uint64_t PM_pos_count = 0;
	uint64_t M = sorted->size();
	uint64_t m1 = next_greater_prime(uint64_t(sqrt(M))); 
	while (M % m1 == 0) m1 = next_greater_prime(m1);
	uint64_t m2 = next_greater_prime(m1+uint64_t(sqrt(m1)));  //19.4.2018 (m1+uint64_t(sqrt(m1)))
	while (M % m2 == 0) m2 = next_greater_prime(m2);

	double random_seq_id = 0.0;
	const map<char, double> prob = get_AA_prob_file(NamedVectorBase::cacheFilename(fasta,"_AAcount.map"));
	for (auto p:prob) random_seq_id += p.second*p.second;
	random_seq_id *= double(lf);

	map<int, uint64_t> frag_shuf;

	cout << "[heuristic_inter_distance] m1 = " << m1 <<  " m2 = " << m2 << endl;
	cerr << "# m1=" << m1 << " m2 = " << m2 << " rep=" << rep << " lf=" << lf << " precision=" << precision << endl;
	cerr << "iteration deviation" << endl;
	indices_t di1, di2, si1=0, si2=0;
	double thresh = 1.0/double(pow(10,precision));
	for (uint64_t i=0; ; ){
		// pick two random points
		si1 = indices_t((si1+m1) % M); 
		di1 = si1; 
		di1 = sorted->at(si1);
		if (string(&data->at(di1)).size() < lf) continue;
		si2 = indices_t((si2+m2) % M);  
		di2 = si2;
		di2 = sorted->at(si2);
		if (string(&data->at(di2)).size() < lf) continue;
		if (di2>=di1){
			if (di2-di1 < lf) continue;
		} else if (di1-di2 < lf) continue;
		++i;
		//cerr << di1 << " " << di2 << " " << double(di1)-double(di2) << endl;
		int PM = PM_neighbor(&(data->at(di1)), &(data->at(di2)), lf);
		if (PM < 0 || PM>lf) continue;
		map<char, double> aa1, aa2;
		for (size_t i =0; i<lf; ++i) { aa1[data->at(di1+i)]++; aa2[data->at(di2+i)]++; }
		double E = 0.0;
		for (char c='A'; c<='Y'; ++c) E += aa1[c]*aa2[c];
		E /= pow(double(lf),2.0);
		frag_shuf[int(E*lf*10)]++;

		for (size_t bla=0; bla<lf; ++bla){
			uint16_t key = (uint16_t(data->at(di1+bla))<<8) + uint16_t(data->at(di2+bla));
			exchanges[PM][key]++;
			//cout << (data->at(di1+bla)) << "->" << data->at(di2+bla) << " " << key << " | " << char(key>>8) << "->" << uint8_t(key) << endl;
		}
		if (PM == 2){
			int start = -1;
			for (size_t bla=0; bla<lf; ++bla) {
				if (data->at(di1+bla) != data->at(di2+bla)) { 
					if (start == -1) { start = bla;
					} else {
						PM_distances[bla-start]++;
						f_2PM_distance.open(dir+"results/intra_dis_" + to_string(lf) + "mer_2PM_distance.txt", ios::out);
						f_2PM_distance << "2PM_distance count" << endl;
						for (size_t bh=1;bh<lf;++bh) f_2PM_distance << bh << " " << PM_distances[bh] << endl;
						f_2PM_distance.close();
						break;
					}
				}
			}
		}
		if (PM <= lf/4 && PM>0) {
			for (size_t bla=0; bla<lf; ++bla) if (data->at(di1+bla) != data->at(di2+bla)) { PM_pos[bla]++; }
			f_very_similar << di1 << " "<< di2 << " " << PM << " ";
			f_very_similar << string(&data->at(di1)).substr(0,lf) << " "<<  string(&data->at(di2)).substr(0,lf) << endl;
			PM_pos_count++;
		}
		distances[i % rep][PM]++;
		// evaluate!
		if (i % 10000000 == 0){
		    f_frag_shuf.open(dir+"results/intra_dis_" + to_string(lf) + "mer_frag_shuf_heuri.csv", ios::out);	
			double s = 0.0; for (auto val:frag_shuf) s += val.second;
			f_frag_shuf << "# Replacement for binominal model" << endl << "PM count fraction" << endl;
			for (auto val:frag_shuf) f_frag_shuf << double(val.first)/10.0 << " " << val.second << " " <<  double(val.second)/s;
			f_frag_shuf.close();
			// [FILE] exchange matrices *************
			f_exchange.open(dir+"results/intra_dis_" + to_string(lf) + "mer_exchange.csv", ios::out);
			f_exchange << "#PM-distance-wise exchange matrix, " << i << " comparisons\nlf PM-distance";
			for (auto ex:exchanges[lf-1]) f_exchange << " " << char((ex.first)>>8) << char(ex.first);
			f_exchange << endl;
			for (int bla=0; bla<=lf; ++bla){
				f_exchange << lf << " " << bla;
				for (auto ex:exchanges[lf-1]) f_exchange << " " << exchanges[bla][ex.first];
				f_exchange << endl;
			}
			f_exchange.close(); // ******************

			// [FILE] position specific mutation **********
			f_intra_dis_pos.open(dir+"results/intra_dis_" + to_string(lf) + "mer_pos.csv", ios::out);
			f_intra_dis_pos << "# total frags with PM_dis <= " << lf/4 << " is " << PM_pos_count<< endl;
			f_intra_dis_pos << "lf position mutation" << endl;
			for (size_t bla=0; bla<lf; ++bla) f_intra_dis_pos << lf << " " <<bla << " " << PM_pos[bla] << endl;
			f_intra_dis_pos.close(); // *******************

			double avg_diff = 0.0, sum = i/rep;
			for (size_t r=0; r<rep; ++r) for (auto val: distances[r]) dis_prob[r][val.first] = double(val.second)/sum;
			size_t zero_sum = 0; for (size_t r1=0; r1<rep; ++r1) zero_sum += dis_prob[r1][0];
			for (size_t r1=0; r1<rep; ++r1){
				for (size_t r2=r1+1; r2<rep; ++r2){
					double diff = 0.0;
					for (size_t j=0; j<=3; ++j) diff += pow(distances[r1][j]-distances[r2][j], 2.0);
					avg_diff += diff;
				}
			}
			avg_diff = avg_diff/double(rep*(rep-1));
			vector<uint64_t> counter(lf+1);
			for (size_t r=0; r<rep; ++r) for (size_t j=0; j<=lf; ++j) counter[j] += distances[r][j];
			uint64_t the_sum = 0;
			for (auto val: counter) the_sum += val;
			// [FILE] intra distances **********
			f_intra_dis.open(dir+"results/intra_dis_" + to_string(lf) + "mer.csv", ios::out);
			f_intra_dis << "# m1=" << m1 << " m2 = " << m2 << " rep=" << rep << " lf=" << lf << " precision=" << precision << " i="<<i<< endl;
			f_intra_dis << "# avg_diff = " << avg_diff << endl;
			f_intra_dis << "lf PM_distance count fraction" << endl;
			for (size_t i=0; i<=lf; ++i) f_intra_dis << lf << " " << i << " "<< counter[i] << " "<< double(counter[i])/double(the_sum) << endl;
			f_intra_dis.close(); // *************
			if (avg_diff < thresh && zero_sum >= 30){
				return;
			}
		}
	}
}

void exhaustive_inter_distance(string fasta, size_t lf){
	auto data = NamedVector<char>::getNamedVector(fasta, "data");
	auto sorted = NamedVector<uint32_t>::getNamedVector(fasta, "sorted");
	auto df = NamedVector<indices_t>::getNamedVector(fasta, "diverse_frags", lf);
	auto common_prefix = NamedVector<uint16_t>::getNamedVector(fasta, "prelen");
	auto string_length = NamedVector<uint16_t>::getNamedVector(fasta, "strlen");
	vector<uint64_t> all_hits;
	map<int,uint64_t> distances;
	uint64_t same=0, n_hits=0, bla=0;
	for (uint32_t dfi=0; dfi<df->size(); ++dfi){
		string frag = string(&data->at(df->at(dfi))).substr(0,lf);
		uint32_t si  = string2sorted(frag, data, sorted);
		uint64_t h = hits(si, lf, common_prefix, string_length);
		all_hits.push_back(h);
		n_hits += h;
	}
	cout << n_hits << " hits captured -> "<< n_hits*n_hits << " 'picks' possible" << endl;
	if (lf == 1) for (uint32_t i=0; i<all_hits.size(); ++i) cout << i << " " << all_hits[i] << endl;
	for (uint32_t dfi=0; dfi<df->size(); ++dfi){
		uint64_t h = all_hits[dfi];
		auto frag = &data->at(df->at(dfi));
		//distances[0] += h*h;
		same += h*h;
		for (uint32_t dfi2=0; dfi2<df->size(); ++dfi2){
			//if (dfi2 == dfi) { distances[0] += h*(h-1); continue; }
			distances[PM_neighbor(frag, &(data->at(df->at(dfi2))), lf)] += h*all_hits[dfi2];
			//distances[1] += h*all_hits[dfi2];
			bla += h*all_hits[dfi2];
		}
	}
	cout << "BLA " << bla << endl;
	uint64_t sum = 0;
	for (auto val:distances) {
		sum += val.second;
		cout << sum << endl; 
	}
	cout << "same " << same << " div/sum^2: " << double(same)/double(n_hits*n_hits) << endl;
	for (auto val:distances) cout << val.first << " " << double(val.second) << " " << double(val.second)/double(sum) <<" " << sum << endl;
}


void predict_subsequent_AAs(string fasta, string protein, size_t lf){
	auto data = NamedVector<char>::getNamedVector(fasta, "data");
	auto sorted = NamedVector<uint32_t>::getNamedVector(fasta, "sorted");
	map<char, uint32_t> entropical_end;
	cerr << "***************" << endl << " * Protein = " << protein << endl;
	float total_certainty = 0.0;
	float real_prob = 0.0;
	size_t not_pred = 0;
	for (size_t i=0; i<protein.size()-lf; ++i){
		string fragment = protein.substr(i, lf);
		char real_AA = protein[i+lf];
		cerr << "Predict follower of : " << fragment << endl;
		subsequent_AAs(data, sorted, fragment, &entropical_end);
		size_t h = 0; for (auto val:entropical_end) h += val.second;
		char pred_AA = 'X';
		for (auto val:entropical_end) {
			cerr << val.first << " : " << val.second*100.0/float(h) << "%" << endl;
			if (val.second > entropical_end[pred_AA]) pred_AA = val.first;
		}
		if (h > 0){
			float certainty = entropical_end[pred_AA]*100.0/float(h);
			total_certainty += certainty;
			real_prob += entropical_end[real_AA]*100.0/float(h);
		}
		if (pred_AA == 'X') not_pred++;

		if (pred_AA == real_AA) cerr << "CORRECT";
		else cerr << "WRONG";
		cerr << " predict " << pred_AA << " is " << real_AA << endl;
	}
	float norm = float(protein.size()-lf-not_pred);
	cerr << protein.size()-lf << " " << not_pred << " " << norm << endl;
	cerr << "Certainty = " << total_certainty/norm << " norm over all: " << total_certainty/float(protein.size()-lf) << " // real prob = " << real_prob/norm << " norm over all: " << real_prob/float(protein.size()-lf) << endl;
}

void conditional_entropy(string fasta, size_t lf, data_t data, sorted_t sorted, prelen_t common_prefix, strlen_t string_length, size_t min_hits, size_t max_hits, bool drop_X){
	//cout << "Hits >= " << min_hits << " and  <= " << max_hits << endl;
	uint32_t h = 0; char ch;
	map<char, uint32_t> entropical_end;
	double H = 0.0;

	indices_t total = 0, total_sum=0;
	for (indices_t sl=0; sl<=string_length->size(); ++sl){
		if (drop_X && contains_char(data, sorted->at(sl), lf, 'X')) continue;
		if (string_length->at(sl) >= lf) ++total;
	}
	vector<pair<double,double>> entropies;
	// find prefix blocks
	size_t context_count = 0;
	for (indices_t si=0; si<sorted->size(); si += h){
		//if (sorted->at(si) > sorted->size()-lf) continue; // should not be a problem
		for (;;){
			h = hits(si,lf-1, common_prefix, string_length); // hits of context c
			if (h == 0) ++si;
			else break;
		}
		// check frag of prefix block
		if (drop_X && contains_char(data, sorted->at(si), lf-1, 'X')) continue; // invalid prefix block
		entropical_end.clear();
		for (indices_t si2=si; si2<si+h; ++si2){
			if (string_length->at(si2) < lf-1) { cerr << "[ERROR] should not happen" << endl; continue;} // should not happen
			ch = data->at(sorted->at(si2)+lf-1); // character ch following context c
			if     (!drop_X && ch != '\0') entropical_end[ch]++;
			else if (drop_X && ch != '\0' && ch != 'X') entropical_end[ch]++;
		}
		indices_t sum = 0;
		for (auto val:entropical_end) sum += val.second;

		context_count++;
		total_sum += sum;
		double ent = entropy(&entropical_end, 2);
		entropies.push_back(make_pair(sum,ent));
		H += sum * ent;
	}
	H /= double(total_sum);

	double variance = 0.0;
	for (auto val:entropies) {
		variance += val.first * (val.second-H) * (val.second-H);
	}
	variance /= double(total_sum);
	//cout << "entropy variance total_sum context_count lf min_hits max_hits" << endl;
	cout << H << " " << variance << " " << total_sum << " " << context_count << " " << lf << " " << min_hits << " " << max_hits << endl;
}






map<uint32_t, uint32_t> frags_per_org(string fasta, data_t data, sorted_t sorted, string fragment, orgstart_t orgstart, bool print_orgs){
	auto interval = std::equal_range(sorted->begin(), sorted->end(), fragment, CompareIndex(data,fragment.size()));
	map<uint32_t, uint32_t> org_hits;
	for (auto si=interval.first; si!=interval.second; ++si){
		//!!! *si is a index in data !
		indices_t di = *si;
		auto low = lower_bound(orgstart->begin(), orgstart->end(), di);
		if (*low != di) --low;
		org_hits[low-orgstart->begin()]++;
	}
	if (print_orgs) for (auto o: org_hits) cout << organism_name(fasta, o.first) << " " << o.second << endl;
	return org_hits;
}


set<uint32_t> orgs_with_frag(string fasta, data_t data, sorted_t sorted, string fragment, orgstart_t orgstart, bool print_orgs){
	auto interval = std::equal_range(sorted->begin(), sorted->end(), fragment, CompareIndex(data,fragment.size()));
	set<uint32_t> orgs;
	for (auto si=interval.first; si!=interval.second; ++si){
		indices_t di = sorted->at(*si);
		auto low = lower_bound(orgstart->begin(), orgstart->end(), di);
		if (*low != di) --low;
		orgs.insert(low-orgstart->begin());
	}
	if (print_orgs) for (auto o: orgs) cout << organism_name(fasta, o) << endl;
	return orgs;
}

int hits(data_t data, sorted_t sorted, string fragment){
	auto interval = std::equal_range(sorted->begin(), sorted->end(), fragment, CompareIndex(data,fragment.size()));
	return interval.second - interval.first; 
}

int hits(const string fasta, string fragment){
	data_t data = NamedVector<char>::getNamedVector(fasta, "data");
	sorted_t sorted = NamedVector<uint32_t>::getNamedVector(fasta, "sorted");

	auto interval = std::equal_range(sorted->begin(), sorted->end(), fragment, CompareIndex(data,fragment.size()));
	//cout << interval.second << "-" << interval.first << " hits of string, interval" << interval.second - interval.first << endl;
	return interval.second - interval.first; 
}

pair<uint32_t, uint32_t> hits_per_prot(uint32_t si, size_t fraglen, sorted_t sorted, prelen_t common_prefix, strlen_t string_length, protstart_t protstart, protcount_t protcount){
	uint32_t h = 0;
	set<uint32_t> prots;
	if (string_length->at(si) >= fraglen) {
		h=1;
		prots.insert(protpos2prot(data2protpos(protstart, sorted->at(si), protcount)));
	}
	else return make_pair(0,0);
	for (size_t i=si; i<common_prefix->size(); ++i){
		//if (i % 100000000 == 0) cout << "  " << i << " cp: "<<  common_prefix->at(i)<<endl;
		//if (si % (sorted->size()/10000) == 0) progress_bar_cerr(si, sorted->size());
		if (string_length->at(i) >= fraglen){
			if (common_prefix->at(i) >= fraglen) {
				++h;
				prots.insert(protpos2prot(data2protpos(protstart, sorted->at(i+1), protcount))); // i should never exeed sorted->size() because last common_prefix == 0;
			} else break;
		}
	}
	return make_pair(prots.size(),h);
}

uint32_t hits(uint32_t si, size_t fraglen, prelen_t common_prefix, strlen_t string_length){
	uint32_t h = 0;
	if (string_length->at(si) >= fraglen) h=1;
	else return 0;
	indices_t i;
	for (i=si; i<common_prefix->size(); ++i){
		if (string_length->at(i) >= fraglen){
			if (common_prefix->at(i) >= fraglen) {
				++h;
			} else {
				//cerr << "BREAK of hits at si=" << i << endl;
				break;
			}
		}
	}
	if (fraglen == 0) --h;
	//cout << "i=" << i << " h=" << h << endl;
	return h;
}


void locality_in_protein(string original_prot, bool randomize){
	map<char,int> waiting;
	map<char,map<int,int>> distances; // char : distance, count
	string prot = original_prot;
	if (randomize) random_shuffle(prot.begin(), prot.end());
	int i=0;
	for (auto ch:prot){
		++i;
		if (!waiting.count(ch)) {
			// ### 1 ### distance to beginning ###
			distances[ch][i]++;
			cout << ch << " "<<i << endl;
		} else {
			// ### 2 ### intermediate distances! ###
			distances[ch][waiting[ch]]++;
			cout << ch << " " << waiting[ch] << endl;
		}
		waiting[ch] = 0;
		for (auto val:waiting) waiting[val.first]++; // all wait one longer
	}
	// ### 3 ### distance to end ###
	for (auto val:waiting) {
		distances[val.first][waiting[val.first]-1]++;
		cout << val.first << " " << waiting[val.first]-1 << endl;
	}

	for(auto ch_count: distances){
		//for (auto val:ch_count.second) cout << ch_count.first << " " << val.first << " " << val.second << endl;
	}

}

void locality_in_proteins(string fasta){
	auto data = NamedVector<char>::getNamedVector(fasta, "data");
	auto sorted = NamedVector<uint32_t>::getNamedVector(fasta, "sorted");
	auto common_prefix = NamedVector<uint16_t>::getNamedVector(fasta, "prelen");
	auto string_length = NamedVector<uint16_t>::getNamedVector(fasta, "strlen");
	auto protstart = NamedVector<uint32_t>::getNamedVector(fasta, "protstart");
	auto protcount = NamedVector<uint32_t>::getNamedVector(fasta, "fasta_prot_count");

	map<char,int> waiting; // waiting time for each character
	map<char,map<int,int>> all_distances; // char : distance, count

	map<char,map<int,int>> distances; // char : distance, count
	char ch;
	map<char,vector<double>> entropies;
	cout << protstart->size() << " = #prots" << endl;
	for (uint32_t ps=0; ps<protstart->size()-1; ps++){
		//cout << " ### PROTEIN " << ps << " ###" << endl;
		//print_protein(fasta, ps);
		waiting.clear();
		distances.clear();
		for (auto di=protstart->at(ps); di < protstart->at(ps+1); ++di){
			ch = data->at(di);
			// ### 1 ### distance to beginning ###
			if (!waiting.count(ch)) {
				//distances[ch] = map<int,int>();
				//distances[ch][di-start+1]++;
				//all_distances[ch][di-start+1]++;
			} else {
			// ### 2 ### intermediate distances! ###
				distances[ch][waiting[ch]]++;
				all_distances[ch][waiting[ch]]++;
			}
			waiting[ch] = 0;
			for (auto val:waiting) waiting[val.first]++;
		}
		// ### 3 ### distance to end ###
		//for (auto val:waiting) {
			//distances[val.first][waiting[val.first]-1]++;
			//all_distances[val.first][waiting[val.first]-1]++;
		//}
		for (auto val:distances){
			entropies[val.first].push_back(entropy(&val.second, 2));
		}
		if (ps % 100000 == 0) { cout << ps << " / " << protstart->size() << endl; }
		//if (ps == 100000) break;
	}
	for(auto ch_count: all_distances){
		for (auto val:ch_count.second) cout << ch_count.first << " " << val.first << " " << val.second << endl;
	}
	for (auto e:entropies) {
		auto v = e.second;
		double sum = std::accumulate(v.begin(), v.end(), 0.0);
		double mean = sum / v.size();
		double sq_sum = std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
		double stdev = sqrt(sq_sum / v.size() - mean * mean);
		cout << e.first << " " << mean << " "<< stdev << endl;

		//for (auto val:e.second) cout << val << " "; cout << endl;
	}
}



map<uint32_t,uint32_t> hit_counts_per_protein(string fasta, size_t lenfrag) {
	data_t data = NamedVector<char>::getNamedVector(fasta, "data");
	sorted_t sorted = NamedVector<uint32_t>::getNamedVector(fasta, "sorted");
	shared_ptr<CachedVector<uint16_t>> common_prefix = NamedVector<uint16_t>::getNamedVector(fasta, "prelen");
	shared_ptr<CachedVector<uint16_t>> string_length = NamedVector<uint16_t>::getNamedVector(fasta, "strlen");
	shared_ptr<CachedVector<uint32_t>> protstart = NamedVector<uint32_t>::getNamedVector(fasta, "protstart");
	shared_ptr<CachedVector<uint32_t>> protcount = NamedVector<uint32_t>::getNamedVector(fasta, "fasta_prot_count");
	map<uint32_t,uint32_t> ret;

	pair<uint32_t, uint32_t> hpp;
	uint32_t si;
	cout << endl;
	for (si=0;si<sorted->size();){
		
		hpp = hits_per_prot(si, lenfrag, sorted, common_prefix, string_length, protstart, protcount);
		if (hpp.first == 0) {++si; continue;}
		if (!(contains_char(data, sorted->at(si), lenfrag, 'X'))) ret[hpp.first]++;
		si += hpp.second;
	}
	return ret;
}

void error_to_exp(string fasta, size_t lf){
	data_t data = NamedVector<char>::getNamedVector(fasta, "data");
	sorted_t sorted = NamedVector<uint32_t>::getNamedVector(fasta, "sorted");
	auto all_valid_hits = NamedVector<uint32_t>::getNamedVector(fasta, "valid_hits");
	const map<char, double> prob = get_AA_prob_file(NamedVectorBase::cacheFilename(fasta,"_AAcount.map"));

	uint64_t si_start = 0;
	string frag = string(&data->at(sorted->at(0))).substr(0,lf); 
	while (frag.length() != lf) {frag = string(&data->at(sorted->at(si_start))).substr(0,lf); ++si_start;}
	double count = 1.0;
	double total_P = 0.0, error_under = 0.0, error_over = 0.0;
	double N = all_valid_hits->at(lf);
	for (uint64_t si=si_start; si<sorted->size(); ++si){
		if (string(&data->at(sorted->at(si))).substr(0,lf).compare(frag) == 0){
			count += 1.0;
		} else {
			double P = frag_prob(&prob, frag);
			double exp_hits = P*N;
			total_P += P;
			if (count > exp_hits) error_over += (count-exp_hits);
			else if (count < exp_hits) error_under += (exp_hits-count);
			frag = string(&data->at(sorted->at(si))).substr(0,lf);
			while (frag.length() != lf) { ++si; frag = string(&data->at(sorted->at(si))).substr(0,lf);}
			count=1;
		}
 		//if (si % 1000000 == 0){
		//	cout << si << ") " << total_P << " " << error_over+error_under << endl;
		//}
	}
	double P = frag_prob(&prob, frag);
	double exp_hits = P*N;
	total_P += P;
	if (count > exp_hits) error_over += (count-exp_hits);
	else if (count < exp_hits) error_under += (exp_hits-count);

	error_under += (1.0-total_P)*N;
	cout << " lf error_under error_over total_prob_visited valid_hits" << endl;
	cout << lf << " " << error_under << " " << error_over << " "<< total_P << " " << N << endl;
}


map<uint32_t,uint32_t> hit_counts(string fasta, size_t lenfrag) {
	cout << "[WARING] not per protein count!" << endl;
	map<uint32_t,uint32_t> ret;
	auto common_prefix = NamedVector<uint16_t>::getNamedVector(fasta, "prelen");
	auto string_length = NamedVector<uint16_t>::getNamedVector(fasta, "strlen");

	uint32_t count = 0;
	if (string_length->at(0) >= lenfrag) count = 1;
	for (size_t i=0; i<common_prefix->size(); ++i){
		auto l = string_length->at(i);
		if (l >= lenfrag){
			if (common_prefix->at(i) >= lenfrag) {
				count++;
			} else {
				ret[count]++;
				count = 1;
			}
		} else count = 1;
	}
	//ret[count]++;
	return ret;
}


bool contains_char(data_t data, indices_t di, size_t lf, char ch){
	for (size_t i=0; i<lf; ++i){
		if (data->at(di+i) == ch) return true;
	}
	return false;
}


vector<string> most_often_occ(string fasta, int fraglen, int min_max_count){
	auto data = NamedVector<char>::getNamedVector(fasta, "data");
	auto sorted = NamedVector<uint32_t>::getNamedVector(fasta, "sorted");
	auto common_prefix = NamedVector<uint16_t>::getNamedVector(fasta, "prelen");
	auto string_length = NamedVector<uint16_t>::getNamedVector(fasta, "strlen");
	auto protstart = NamedVector<uint32_t>::getNamedVector(fasta, "protstart");
	auto protcount = NamedVector<uint32_t>::getNamedVector(fasta, "fasta_prot_count");

	uint32_t max_count = min_max_count;
	vector<string> ret;
	string s;
	uint32_t si = 0;
	pair<uint32_t, uint32_t> hpp;
	for (;si<sorted->size();){
		hpp = hits_per_prot(si, fraglen, sorted, common_prefix, string_length, protstart, protcount);
		if (hpp.first == 0) {++si; continue;}
		if (hpp.first > max_count){
			max_count = hpp.first;
			ret.clear();
		}
		if (hpp.first == max_count){
			s = (&(data->at(sorted->at(si))));
			s = s.substr(0,fraglen);
			ret.push_back(s);
		}
		si += hpp.second;
		if (si % (sorted->size()/1000) == 0) progress_bar_cerr(si, sorted->size());
	}

	if (max_count == 5 && ret.size() == 0) cout << "[WARNING] mosoften_occ_t was set to be min "<< min_max_count << "!" << endl;
	return ret;
}





// mmusl
map<uint32_t,uint64_t> maximal_minimal_unique_string_length(string fasta) {
	map<uint32_t,uint64_t> ret;
	auto data = NamedVector<char>::getNamedVector(fasta, "data");
	auto sorted = NamedVector<uint32_t>::getNamedVector(fasta, "sorted");
	auto common_prefix = NamedVector<uint16_t>::getNamedVector(fasta, "prelen");
	auto protstart = NamedVector<uint32_t>::getNamedVector(fasta, "protstart");
	auto protcount = NamedVector<uint32_t>::getNamedVector(fasta, "fasta_prot_count");
	auto orgstart = NamedVector<uint32_t>::getNamedVector(fasta, "orgstart");
	auto string_length = NamedVector<uint16_t>::getNamedVector(fasta, "strlen");
	cout << "preparations ok" << endl;
	map<int, uint64_t> ortho; // dup_size -> o->count / p->count
	map<int, uint64_t> para;
	map<uint32_t, uint64_t> counter;
	uint32_t con = 0, zero_left_enlonged = 0, left_max = 0, look_back = 0;

	for (indices_t i=0; i<sorted->size(); ++i){
		//if (i == 100000000) break;
		if (i % (sorted->size()/1000) == 0) {
			cerr << "(" << i << ")" << endl;
		} //progress_bar(i,sorted->size());
		const indices_t s = sorted->at(i); // index in data
		const indices_t dup_string_size = common_prefix->at(i);

		if (dup_string_size > 0){ // if unique_suffix -> kommt nur von einem string! <-enlongment works
			// X CHECK
			if (contains_char(data, s, dup_string_size, 'X')) { continue; }
			// FIRST DUP CHECK -> look back
			if (i > 0 && common_prefix->at(i-1) >= dup_string_size) {
				look_back++;
				continue;
			}
			const indices_t first_dup_index = i; // first duplicate

			// OVERLAP CHECK
			uint32_t diff;
			indices_t cp = first_dup_index+1;
			for (;;){ // for all duplications

				if (sorted->at(cp) > s) diff = sorted->at(cp)-s;
				else diff = s-sorted->at(cp);

				if (diff >= dup_string_size){
					cp = indices_t(-1); // non-overlapping duplicate found!!!
					break; // overlap check
				}
				if (common_prefix->at(cp) < dup_string_size) break; // last duplicate
				++cp;
			}

			if (cp != indices_t(-1)) { // only duplicate is overlapping 
				con++; //cerr << string(&(data->at(s))).substr(0,dup_string_size) << " "<< i << endl; 
				continue; //!!!
			}
			// LEFT MAXIMAL CHECK
			if (s>0 && data->at(s-1) != '\0' && data->at(s-1) != 'X'){
				string f = string(&(data->at(s-1)));
				uint32_t longet_dup_sindex = string2sorted(f, data, sorted);
				indices_t enlonged_dup_string_size = common_prefix->at(longet_dup_sindex);

				if (dup_string_size+1 == enlonged_dup_string_size) {
					left_max++;
					continue; // only take if maximal to left
				}
			} else { // cannot be left-enlonged
				indices_t cp = first_dup_index;
				for (;;){ //-> look of any duplicate that can be left enlonged
					if (cp >= sorted->size()) break;
					
					if (sorted->at(cp) > 0){
						if (data->at(sorted->at(cp)-1) && data->at(sorted->at(cp)-1) != 'X'){ // other duplicate found that can be left enlonged
							string f = string(&(data->at(sorted->at(cp)-1)));
							uint32_t longet_dup_sindex = string2sorted(f, data, sorted);
							indices_t enlonged_dup_string_size = common_prefix->at(longet_dup_sindex);
							//indices_t enlonged_unique_string_size = is_unique_string_of_size(string2sorted(f ,data, sorted), data, sorted, common_prefix);
							if (dup_string_size+1 == enlonged_dup_string_size){
								cp = indices_t(-1); break;
							}
						}
					}
					if (!(common_prefix->at(cp) >= dup_string_size)) break; // last duplicate
					++cp;
				}
				if (cp == indices_t(-1)) {
					zero_left_enlonged++;
					continue;
				}
			}

			// TAKE !!!
			string dup = string(&(data->at(s))).substr(0,dup_string_size);
			pair<double,double> ortho_para = ortholog_paralog(fasta,dup,data,sorted,orgstart);

			ortho[dup_string_size] += ortho_para.first;
			para[dup_string_size] += ortho_para.second;
			counter[dup_string_size] ++;
		}
	}
	cout << "look back: " << look_back << endl;
	cout << "left_max: " << left_max << endl;
	cout << "con because of overlap: " << con << endl; 
	cout << "zero left enlonged: " << zero_left_enlonged << endl; 
	uint64_t mask = (uint64_t(1) << 32)-1;
	cout << "length n_strings ortho para" << endl;
	for (auto s:ortho){ // dup_string_size->#orthos
		int strlength = s.first;
		cout << strlength << " " << counter[strlength] << " " << s.second << " " << para[strlength] << endl;
	}

	//ortho_para.print();
	return counter;
}










pair<double,double> ortholog_paralog(string fasta, string frag, data_t data, sorted_t sorted, sorted_t orgstart){
	map<uint32_t, uint32_t> org_count = frags_per_org(fasta, data, sorted, frag, orgstart, false); // limiting factor
	uint32_t T = 0; for (auto o:org_count) T += o.second;
	double ort = 0.0, par = 0.0;
	for (auto o:org_count) {
		ort += double(o.second*(T-o.second));
		if (o.second > 1) par += o.second*(o.second-1.0);
	}
	ort = ort/2.0;
	par = par/2.0;
	return make_pair(ort,par);
}




// sorted->at(si) = index in data
size_t is_unique_string_of_size(indices_t si, data_t data, sorted_t sorted, shared_ptr<CachedVector<uint16_t>> common_prefix){
	auto d = sorted->at(si);
	auto o = max( (si>0)?common_prefix->at(si-1):common_prefix->at(si), common_prefix->at(si)); // check for overlapp with previous and following sequence
	return (d+o < data->size() && data->at(d+o)) ? o+1 : 0 ;
}


//occurs in only one prot, min length
// smalles form of 
size_t is_unique_string_per_proof_size_t(uint32_t si, data_t data, sorted_t sorted, shared_ptr<CachedVector<uint16_t>> common_prefix){
	//auto s = sorted->at(si);
	auto o = max((si>0)?common_prefix->at(si-1):common_prefix->at(si), common_prefix->at(si));

	size_t length = o-1;
	set<uint32_t> pos;

	for( ;; --length){
		// go doen and up from si, as long as overlap>= length XXX
	}
	return o;
}

indices_t string2df(string& f, data_t data, diverse_frags_t diverse_frags, CachedVector<indices_t>::iterator b, CachedVector<indices_t>::iterator e){
	auto res = lower_bound(b, e, f, [&](auto it, string val){
		return string(&data->at(it)).substr(0,f.size()) < val;
	});
	if (string(&data->at(*res)).substr(0,f.size()).compare(f) == 0) return res-diverse_frags->begin();
	else return indices_t(-1);
}


// return si
// find index where fragment matches data->at(sorted->at(INDEX))
// corrected 10.8.17 XXX
indices_t string2sorted(string& fragment, data_t data, sorted_t sorted){
	auto res = lower_bound(sorted->begin(), sorted->end(), fragment, [&](auto it, string val){
		return string(&data->at(it)).substr(0,fragment.size()).compare(val) < 0; // changed to <= 5.12.17
	});
	if (string(&data->at(*res)).substr(0,fragment.size()).compare(fragment) == 0) return res-sorted->begin();
	else return indices_t(-1);
}


void inter_distance(string fasta, string fasta_rand, int lf, indices_t start, indices_t stop, int step, string name){
	auto data_n = NamedVector<char>::getNamedVector(fasta, "data");

	auto data_r = NamedVector<char>::getNamedVector(fasta_rand, "data");
	auto sorted = NamedVector<uint32_t>::getNamedVector(fasta_rand, "sorted");
	auto string_length = NamedVector<uint16_t>::getNamedVector(fasta_rand, "strlen");

	map<int, uint64_t> all_counter;
	map<int, uint64_t> tmp_counter;
	map<uint64_t, uint64_t> dis_lumpsize_counter;

	int pm;
	
	uint64_t m = next_greater_prime(uint64_t(sqrt(double(data_n->size())))); // modified 25.10.

	cout << "[inter_dictance_to_random] m = " << m << endl;
	indices_t di;
	for (uint64_t di_h=start; di_h<stop; ++di_h){
		di = indices_t((di_h*m) % data_n->size()); // modified 17.10.
		if (string(&data_n->at(di)).size() < size_t(lf) || contains_char(data_n, di, lf, 'X')) continue;
		if (di % step == 0) cerr << "di= " << di << endl;
		tmp_counter.clear();
		// natural data point determined
		indices_t si=0; if (step > 1) si = (di % step); // modified 24.10. from "> 0" to "> 1" XXX
		for (; si<sorted->size(); si += step){
			if (string_length->at(si) < lf || sorted->at(si) == di) continue;
			pm = PM_neighbor(&data_n->at(di), &(data_r->at(sorted->at(si))), lf); if (pm < 0) continue;
			tmp_counter[pm]++;
		}
		uint64_t sum = 0;
		int ind_size = tmp_counter.begin()->first;
		for (auto val:tmp_counter) {
			//distance = val.first // count = val.second
			all_counter[val.first] += val.second;
			sum += val.second; // current lump size
			while (ind_size <= val.first){
				dis_lumpsize_counter[(uint64_t(ind_size)<<56)+sum]++;
				++ind_size;
			}
		}
	}
	string fn = NamedVector<char>::cacheFilename(fasta, "/results/inter" + name + "_" + to_string(lf) + "_" + to_string(start) + "_" + to_string(stop) +"_" + to_string(step) +".txt");
	ofstream fs; fs.open(fn);
	fs << "# " << fasta << " vs " << fasta_rand << " m = "<<  m << " start = " << start << " stop = " << stop << " lf = " << lf << " step = " << step << endl;
	fs << "distance count" << endl;
	for (int i=0; i<=lf; ++i) fs << i << " " << all_counter[i] << endl;
	fs << "distance lumpsize count" << endl;
	for (auto val:dis_lumpsize_counter) fs << (val.first>>56) << " " << (val.first&((uint64_t(1)<<56)-1)) << " " << val.second << endl; 
	fs.close();
}



//XXX ##########################################
void intra_distance(string fasta, int lf, indices_t start, indices_t stop, int step){
	auto data = NamedVector<char>::getNamedVector(fasta, "data");
	auto sorted = NamedVector<uint32_t>::getNamedVector(fasta, "sorted");
	auto string_length = NamedVector<uint16_t>::getNamedVector(fasta, "strlen");

	map<int, uint64_t> all_counter;
	map<int, uint64_t> tmp_counter;
	map<uint64_t, uint64_t> dis_lumpsize_counter;

	int pm;
	uint64_t m = next_greater_prime(uint64_t(sqrt(double(data->size()))));

	cout << "[intra_dictance] m = " << m << endl;
	//for (indices_t di=start; di<stop; ++di){
	indices_t di;
	for (uint64_t di_h=start; di_h<stop; ++di_h){
		di = indices_t((di_h*m) % data->size()); // modified 17.10.
		if (string(&data->at(di)).size() < size_t(lf) || contains_char(data, di, lf, 'X')) continue;
		if (di % step == 0) cerr << "di= " << di << endl;
		tmp_counter.clear();
		indices_t si=0; if (step > 0) si = (di % step);
		for (; si<sorted->size(); si += step){
			if (string_length->at(si) < lf || sorted->at(si) == di) continue;
			pm = PM_neighbor(&data->at(di), &(data->at(sorted->at(si))), lf); if (pm < 0) continue;
			tmp_counter[pm]++;
		}
		uint64_t sum = 0;
		int ind_size = tmp_counter.begin()->first;
		for (auto val:tmp_counter) {
			//distance = val.first // count = val.second
			all_counter[val.first] += val.second;
			sum += val.second; // current lump size
			while (ind_size <= val.first){
				dis_lumpsize_counter[(uint64_t(ind_size)<<56)+sum]++;
				++ind_size;
			}
		}
	}
	string fn = NamedVector<char>::cacheFilename(fasta, "/results/intra" + to_string(lf) + "_" + to_string(start) + "_" + to_string(stop) +"_" + to_string(step) +".txt");
	ofstream fs; fs.open(fn);
	//fs << "# " << fasta << endl;
	fs << "# " << fasta << " m = " << m << " start = " << start << " stop = " << stop << " lf = " << lf << " step = " << step << endl;
	fs << "distance count" << endl;
	for (int i=0; i<=lf; ++i) fs << i << " " << all_counter[i] << endl;
	fs << "distance lumpsize count" << endl;
	for (auto val:dis_lumpsize_counter) fs << (val.first>>56) << " " << (val.first&((uint64_t(1)<<56)-1)) << " " << val.second << endl; 
	fs.close();
}



map<uint32_t,uint64_t> minimal_unique_strings(string fasta) {
	map<uint32_t,uint64_t> ret;
	auto data = NamedVector<char>::getNamedVector(fasta, "data");
	auto sorted = NamedVector<uint32_t>::getNamedVector(fasta, "sorted");
	auto common_prefix = NamedVector<uint16_t>::getNamedVector(fasta, "prelen");

	if (data->at(sorted->at(0)+1)) ret[common_prefix->at(0)+1]++;
	cout << &data->at(sorted->at(0)) << " " << common_prefix->at(0)+1<< endl;
	for (size_t i=1; i<sorted->size(); ++i){
		auto s = sorted->at(i);
		auto o = max(common_prefix->at(i-1), common_prefix->at(i));
		if (data->at(s+o)) {// unique following character
			ret[o+1]++;
			if (o<=3) cout << &data->at(sorted->at(i)) << " " << o+1<< endl;
		}
	}
	return ret;
}


void print_vectors(string fasta){
	cout << "*** void prinvectors_t(string fasta) *** start" << endl;
	data_t data = NamedVector<char>::getNamedVector(fasta, "data");
	sorted_t sorted = NamedVector<uint32_t>::getNamedVector(fasta, "sorted");
	shared_ptr<CachedVector<uint16_t>> common_prefix = NamedVector<uint16_t>::getNamedVector(fasta, "prelen");
	shared_ptr<CachedVector<uint16_t>> string_length = NamedVector<uint16_t>::getNamedVector(fasta, "strlen");
	shared_ptr<CachedVector<uint32_t>> protstart = NamedVector<uint32_t>::getNamedVector(fasta, "protstart");
	shared_ptr<CachedVector<uint32_t>> protcount = NamedVector<uint32_t>::getNamedVector(fasta, "fasta_prot_count");

	size_t add = 0;
	size_t start = 100000000;
	for (size_t i=start; i<sorted->size() && i<start+100000+add; ++i){
		if (common_prefix->at(i) >=10){
		cout << right << std::setw(5) << i<<") ";
		cout << left << std::setw(30) << string(&(data->at(sorted->at(i)))).substr(0,30) <<" ";
		//cout << string(&(data->at(sorted->at(i)))).substr(0,max_prelen+1) << " ";
		cout << " data:" <<  right << std::setw(10) << sorted->at(i);
		cout << " prefix:" <<  std::setw(3) <<(int) common_prefix->at(i);
		pair<uint32_t, uint32_t> protpos = data2protpos(protstart, sorted->at(i), protcount);
		cout << " organism/prot/pos:" << std::setw(4) << (protpos.second&((1<<12)-1)) << " / " << std::setw(6) <<protpos.first << "/" << std::setw(5) << (protpos.second>>12);
		cout << " strlen:" << std::setw(6) << string_length->at(i) << endl;
		} else add++;
	}
	cout << "*** void prinvectors_t(string fasta) *** end" << endl;
}


void print_protein_with_fragment(string fasta, string frag, bool print_prot, bool pdb){
	auto data = NamedVector<char>::getNamedVector(fasta, "data");
	sorted_t sorted = NamedVector<uint32_t>::getNamedVector(fasta, "sorted");
	auto protstart = NamedVector<uint32_t>::getNamedVector(fasta, "protstart");
	auto protcount = NamedVector<uint32_t>::getNamedVector(fasta, "fasta_prot_count");
	auto protein_names = NamedVector<uint32_t>::getNamedVector(fasta, "protein_names");
	indices_t start = string2sorted(frag, data, sorted);
	ifstream fprotein; fprotein.open(NamedVector<uint32_t>::cacheFilename(fasta, "_protein_names.list"));
	if (start >=sorted->size()) return;
	for (indices_t si = start; si<sorted->size(); ++si){
		if (frag.compare(string(&data->at(sorted->at(si))).substr(0,frag.size())) != 0) break;
		pair<uint32_t, uint32_t> protpos = data2protpos(protstart, sorted->at(si), protcount);
		uint32_t protein_number = protpos2prot(protpos);
		if (pdb) cout << full_protein_name(fasta, protpos.first, &fprotein, protein_names).substr(1,4) << endl;
		else cout << full_protein_name(fasta, protpos.first, &fprotein, protein_names) << " : " << protpos.first << " pos="<< protpos2pos(protpos)<< endl;
		string the_protein = protein(fasta, protein_number);
		if (print_prot) cout << the_protein.substr(0,protpos2pos(protpos)) << endl << frag << endl << the_protein.substr(protpos2pos(protpos)+frag.size(), -1) << endl;
	}
}


/*
 --- AAAAADWVGGALKGEFGQQQTTGQIVFDATISMFPVLGEGTAARDTIAIMMHMADNRQAVEDRWTWIKLVLCLIAVVPIFGGVLKGVGKLVIRAFEKSEDLAKLAEELVLFVNRMCHGNAYEWLRKLDFTTYQGKVIDGLTDALDRFSGACQYIVRNMGNVLPSHVVAYLSGMPAKLQPIRNAANRMIPQALRDMNNCLA (2/2) orgs=1 :2 --- 
 Paraburkholderia /  Uncharacterized protein(#3098410) / 26 AAAAADWVGGALKGEFGQQQTTGQI
 Paraburkholderia /  Uncharacterized protein(#3098508) / 26 AAAAADWVGGALKGEFGQQQTTGQI
*/

void filter_and_where_from(string fasta, Bounds b_fraglen, Bounds b_hits, Bounds b_orgs, bool print_all){
	auto data = NamedVector<char>::getNamedVector(fasta, "data");
	sorted_t sorted = NamedVector<uint32_t>::getNamedVector(fasta, "sorted");
	auto common_prefix = NamedVector<uint16_t>::getNamedVector(fasta, "prelen");
	auto string_length = NamedVector<uint16_t>::getNamedVector(fasta, "strlen");
	auto protstart = NamedVector<uint32_t>::getNamedVector(fasta, "protstart");
	auto protcount = NamedVector<uint32_t>::getNamedVector(fasta, "fasta_prot_count");
	auto protein_names = NamedVector<uint32_t>::getNamedVector(fasta, "protein_names");
	ifstream fprotein; fprotein.open(NamedVector<uint32_t>::cacheFilename(fasta, "_protein_names.list"));

	cout << "filter_and_where_from " << b_fraglen.first << "-" << b_fraglen.second << endl;
	cout << "organism / protein_name (#protein_number) / position_in_protein" << endl;
	cerr << "hits_per_prot total_hits n_organisms" << endl;
	cerr << "min/max hits = " << b_hits.first << "/" << b_hits.second << endl;
	cerr << "min/max orgs = " << b_orgs.first << "/" << b_orgs.second << endl;
	pair<uint32_t, uint32_t> hpp;
	// *** Bounds b_fraglen loop ***
	for (size_t f=b_fraglen.first; f<=b_fraglen.second; ++f){
		cout << "### fraglen = " << f <<" / of max = " << b_fraglen.second << " ###"<< endl;
		map<uint32_t, map<uint32_t,uint32_t>> counter;
		for (size_t i=0; i<sorted->size(); ){
			//if (i % (sorted->size()/1000) == 0) progress_bar_cerr(i,sorted->size());
			hpp = hits_per_prot(i, f, sorted, common_prefix, string_length, protstart, protcount);
			// *** Bounds b_hits check ***
			if (hpp.second >= b_hits.first && hpp.second <= b_hits.second) {
				set<uint32_t> orgs;
				vector<pair<uint32_t, uint32_t>> propos_stack_t;
				for (size_t j=i; j<i+hpp.second; ++j){
					pair<uint32_t, uint32_t> protpos = data2protpos(protstart, sorted->at(j), protcount);
					uint32_t organism = protpos2organism(protpos);
					orgs.insert(organism);
					propos_stack_t.push_back(protpos);
				}
				// *** Bounds b_orgs check ***
				if (orgs.size()>=b_orgs.first && orgs.size()<=b_orgs.second){
					if (counter.count(hpp.first) == 0) counter[hpp.first] = map<uint32_t,uint32_t>();
					counter[hpp.first][orgs.size()] ++;
					if (print_all) {
						cerr << hpp.first << " " << hpp.second << " " << orgs.size() << endl;
						cout << " --- "<< string(&data->at(sorted->at(i))).substr(0,f) << " (" << hpp.first<<"/"<<hpp.second << ") orgs=" << orgs.size()<< " ("<< hpp.first<< "/" << hpp.second <<") --- " << endl;
						for (size_t j=i; j<i+hpp.second; ++j){
							pair<uint32_t, uint32_t> protpos = propos_stack_t[j-i];
							uint32_t organism = protpos2organism(protpos);
							uint32_t protein = protpos2prot(protpos);
							uint32_t pos = protpos2pos(protpos);
							cout << organism << ":" << organism_name(fasta, organism) << " / " << pretty_protein_name(fasta,protein, &fprotein, protein_names) << " / " << pos <<  " " << endl;
						}
					}
				}
			}
			i += hpp.second;
			if (hpp.second == 0) ++i;
		}
		for (auto hits_map:counter){
			for (auto org_count:hits_map.second) cout << f << " " << hits_map.first << " " << org_count.first << " " << org_count.second << endl;
		}
	}
	fprotein.close();
}

string organism_name(string fasta, uint32_t org){
	return file_content(NamedVector<int>::cacheFilename(fasta, "/_fasta_labels/" + to_string(org) + ".txt"));
}



string full_protein_name(string fasta, uint32_t protein, ifstream* fprotein, protein_names_t pn){
	string ugly_protein;
	fprotein->seekg(pn->at(protein));
	getline(*fprotein, ugly_protein);
	return ugly_protein;
}

string pretty_protein_name(string fasta, uint32_t protein, ifstream* fprotein, protein_names_t pn){
	string ugly_protein;
	fprotein->seekg(pn->at(protein));
	getline(*fprotein, ugly_protein);
	int head = ugly_protein.find_first_of(" ");
	int tail = ugly_protein.find_first_of("=")-3; // " OS="
	return ugly_protein.substr(head, tail-head) + "(#"+to_string(protein)+")";
}


uint32_t protpos2organism(pair<uint32_t, uint32_t> p){
	return p.second&((1<<12)-1);
}

uint32_t protpos2prot(pair<uint32_t, uint32_t> p){
	return p.first;
}

uint32_t protpos2pos(pair<uint32_t, uint32_t> p){
	return p.second>>12;
}


// corrected 16.5.2018
pair<uint32_t, uint32_t> data2protpos(protstart_t protstart, uint32_t pos, protcount_t protcount){
	auto prot_start = upper_bound(protstart->begin(), protstart->end(), pos); prot_start--;
	uint32_t prot = prot_start-protstart->begin();
	auto org_start = upper_bound(protcount->begin(), protcount->end(), prot); // corrected 20.8.2018 org_start-- not necessary! protcount holds the number of protein until there 
	//cerr << "pos = " << pos << " *prot_start=" << *prot_start << " *org_start=" << *org_start << " org=" << org_start-protcount->begin() << endl;
	return make_pair(prot, ((pos-*prot_start)<<12)+ org_start-protcount->begin());
}

void intra_prot_aa_dist(uint32_t org, string fasta, string filename){
	auto data = NamedVector<char>::getNamedVector(fasta, "data");
	auto protstart = NamedVector<uint32_t>::getNamedVector(fasta, "protstart");
	auto protcount = NamedVector<uint32_t>::getNamedVector(fasta, "fasta_prot_count");
	ofstream out; out.open(filename);
	out << "di_start,di_stop,";	for (char c='A';c<='Z';++c) out << c << ","; out << "sum" << endl;
	uint32_t p=0; if (org != 0) p=protstart->at(protcount->at(org-1));
	auto m = amino_acid_distribution(fasta,p, protstart->at(protcount->at(org)), data);
	uint32_t sum = 0;

	uint32_t startprot = 0; if (org != 0) startprot=protcount->at(org-1);
	uint32_t stopprot = protcount->at(org);
	cerr << startprot << "=startprot/stopprot=" << stopprot << endl;
	for (uint32_t i=startprot; i<stopprot; ++i){
		auto m = amino_acid_distribution(fasta, protstart->at(i), protstart->at(i+1), data);
		out << protstart->at(i) << "," << protstart->at(i+1) << ",";
		sum = 0;
		for (auto e:m) {
			out << e.second << ",";
			sum += e.second;
		}
		out << sum << ","<< organism_name(fasta, org) << endl;
	}
	out.close();
}

void intra_prot_aa_dist(string fasta, string filename){
	auto data = NamedVector<char>::getNamedVector(fasta, "data");
	auto protstart = NamedVector<uint32_t>::getNamedVector(fasta, "protstart");
	ofstream out; out.open(filename);
	uint32_t sum;
	out << "di_start,di_stop,";
	for (char c='A';c<='Z';++c) out << c << ",";
	out << "sum" << endl;
	for (uint32_t i=0; i<protstart->size()-1; ++i){
		auto m = amino_acid_distribution(fasta, protstart->at(i), protstart->at(i+1), data);
		out << protstart->at(i) << "," << protstart->at(i+1) << ",";
		sum = 0;
		for (auto e:m) {
			out << e.second << ",";
			sum += e.second;
		}
		out << sum << endl;
	}
	out.close();
}

void intra_organism_aa_dist(string fasta, string filename){
	auto data = NamedVector<char>::getNamedVector(fasta, "data");
	auto protstart = NamedVector<uint32_t>::getNamedVector(fasta, "protstart");
	auto protcount = NamedVector<uint32_t>::getNamedVector(fasta, "fasta_prot_count");
	ofstream out; out.open(filename);
	uint32_t sum;
	out << "di_start,di_stop,";
	for (char c='A';c<='Z';++c) out << c << ",";
	out << "sum,organism" << endl;
	size_t p=0;
	for(size_t org=0; org<protcount->size()-1; ++org){
		out << p << "," << protstart->at(protcount->at(org)) << ",";
		auto m = amino_acid_distribution(fasta,p, protstart->at(protcount->at(org)), data);
		sum = 0;
		for (auto e:m) {
			out << e.second << ",";
			sum += e.second;
		}
		out << sum << ","<< organism_name(fasta, org) << endl;
		p=protstart->at(protcount->at(org));
	}
	out.close();
}

map<char, uint32_t> amino_acid_distribution(string fasta, indices_t di_start, indices_t di_stop, data_t data){
	map<char,uint32_t> ret;
	for (char c='A';c<='Z';++c) ret[c] = 0;
	char ch;
	for (indices_t i=di_start; i<di_stop; ++i){
		ch = data->at(i);
		if (ch) ret[ch]++;
	}
	return ret;
}


string protein_at(uint32_t di, data_t data, sorted_t sorted, protstart_t protstart, protcount_t protcount){
	pair<uint32_t, uint32_t> protpos = data2protpos(protstart, di, protcount);
	string protein = string(&data->at(protstart->at(protpos2prot(protpos))));
	return protein;
}


string protein(string fasta, uint32_t protein_number){
	auto data = NamedVector<char>::getNamedVector(fasta, "data");
	auto protstart = NamedVector<uint32_t>::getNamedVector(fasta, "protstart");
	return string(&data->at(protstart->at(protein_number)));
}




void print_protein(string fasta, uint32_t protein_number){
	auto data = NamedVector<char>::getNamedVector(fasta, "data");
	auto protstart = NamedVector<uint32_t>::getNamedVector(fasta, "protstart");
	auto protcount = NamedVector<uint32_t>::getNamedVector(fasta, "fasta_prot_count");
	auto protein_names = NamedVector<uint32_t>::getNamedVector(fasta, "protein_names");

	pair<uint32_t, uint32_t> protpos = data2protpos(protstart, protstart->at(protein_number), protcount);
	ifstream fprotein; fprotein.open(NamedVector<uint32_t>::cacheFilename(fasta, "_protein_names.list"));
	cout << "> Protein id=" << protein_number << " / " << pretty_protein_name(fasta, protein_number, &fprotein, protein_names) << " / "<< organism_name(fasta, protpos2organism(protpos)) << endl;
	cout << protein(fasta, protein_number) << endl;
}

void ANA_lf_recur_dup(string fasta){
	auto data = NamedVector<char>::getNamedVector(fasta, "data");
	auto sorted = NamedVector<uint32_t>::getNamedVector(fasta, "sorted");
	auto all_valid_hits = NamedVector<uint32_t>::getNamedVector(fasta, "valid_hits");
	cout << "###\nlf hits\n";
	for (size_t i=0; i<all_valid_hits->size(); ++i) cout << i << " " << all_valid_hits->at(i) << endl;
	
	auto common_prefix = NamedVector<uint16_t>::getNamedVector(fasta, "prelen");
	auto string_length = NamedVector<uint16_t>::getNamedVector(fasta, "strlen"); // used for frags occurring once
	map<int, map<int,int>> lf_abu_count;
	map<int,int> lf_hits;

	for (indices_t i=0; i<common_prefix->size(); ++i){
		uint16_t cp = common_prefix->at(i);
		if (lf_hits.size() >= cp){
			auto it = lf_hits.begin(); advance(it, cp);
			if (it != lf_hits.end()){
				for ( ;it != lf_hits.end(); ++it){
					if (lf_abu_count[it->first].size() == 0) lf_abu_count[it->first] = map<int,int>();
					// WRITE
					lf_abu_count[it->first][it->second+1]++;
				}
				it=lf_hits.begin(); advance(it, cp);
				lf_hits.erase(it, lf_hits.end());
			}
		}
		for (size_t j=max(size_t(cp+1),size_t(((i>0)?common_prefix->at(i-1):0)+1)); j<=string_length->at(i); ++j) lf_abu_count[j][1]++; //XXX fixed? 3.4.
		for (size_t j=1; j<=cp; ++j) lf_hits[j]++;
	}
	// end of clean ###
	for (auto it = lf_hits.begin() ;it != lf_hits.end(); ++it){
		if (lf_abu_count[it->first].size() == 0) lf_abu_count[it->first] = map<int,int>();
		// WRITE
		lf_abu_count[it->first][it->second+1]++;
	}
	cout << "#######\nlf abu count" << endl;
	for (auto val:lf_abu_count) {
		for (auto v:val.second){
			cout << val.first << " " << v.first << " " << v.second << endl;
		}
	}

	/*
	// start from longest fragments
	for (auto lf_map=lf_abu_count.rbegin(); lf_map!=lf_abu_count.rend(); ++lf_map){
		// start from greatest abundance
		for (auto abu_count=(lf_map->second).rbegin();  abu_count!=(lf_map->second).rend(); ++abu_count){
			// if abundace >= 2 -> duplicate!
			if (abu_count->first >= 2){
				int c = abu_count->second; // number of diverse duplicates
				for (int l=1;l<lf_map->first;++l){
					for (auto p:lf_abu_count[l]){
						lf_abu_count[l][p->first] -= c;
					}
				}
			}
		}
	}
	*/
}

void ANA_lf_abu_frags(string fasta, size_t lf){
	auto data = NamedVector<char>::getNamedVector(fasta, "data");
	auto sorted = NamedVector<uint32_t>::getNamedVector(fasta, "sorted");


}


void ANA_frag_intra_inter(string fasta, size_t lf, int count){

}

void ANA_frag_exp_abu(string fasta, size_t lf, int N, int take_every){
	auto data = NamedVector<char>::getNamedVector(fasta, "data");
	auto sorted = NamedVector<uint32_t>::getNamedVector(fasta, "sorted");
	auto common_prefix = NamedVector<uint16_t>::getNamedVector(fasta, "prelen");
	auto string_length = NamedVector<uint16_t>::getNamedVector(fasta, "strlen");
	auto protstart = NamedVector<uint32_t>::getNamedVector(fasta, "protstart");
	auto protcount = NamedVector<uint32_t>::getNamedVector(fasta, "fasta_prot_count");
	auto all_valid_hits = NamedVector<uint32_t>::getNamedVector(fasta, "valid_hits");
	const map<char, double> prob = get_AA_prob_file(NamedVectorBase::cacheFilename(fasta,"_AAcount.map"));
	for (auto p:prob) cout << p.first << " " << p.second << endl;
	double DB = all_valid_hits->at(lf);
	for (indices_t i=0; i<sorted->size(); ){
		while (string_length->at(i) < lf) ++i; 
		string frag = string(&data->at(sorted->at(i))).substr(0,lf);

		auto hpp = hits_per_prot(i, lf, sorted, common_prefix, string_length, protstart, protcount);
		double exp_hits = frag_prob(&prob, frag)*DB;
		//if (hpp.second > exp_hits*N && i%take_every==0)
		cout << frag << " " << exp_hits << " " << hpp.second << " " << double(hpp.second)/exp_hits << endl;
		if (hpp.second==0) ++i;
		else i += hpp.second;
		//cerr << i << endl;
	}

}


void print_strings_longer_then(string fasta, size_t lf){
	auto sorted = NamedVector<uint32_t>::getNamedVector(fasta, "sorted");
	auto protstart = NamedVector<uint32_t>::getNamedVector(fasta, "protstart");
	auto protcount = NamedVector<uint32_t>::getNamedVector(fasta, "fasta_prot_count");
	auto string_length = NamedVector<uint16_t>::getNamedVector(fasta, "strlen");
	set<uint32_t> all_prots;
	for (indices_t i=0; i<string_length->size(); ++i){
		if (string_length->at(i) >= lf ){
			auto p = data2protpos(protstart, sorted->at(i), protcount);
			auto prot = protpos2prot(p);
			if (all_prots.find(prot) == all_prots.end()){
				all_prots.insert(prot);
				cout << prot << "=protein / length=" << string_length->at(i) << endl;
				print_protein(fasta, prot);
			}
		}
	}

}

void exp_ortho_para_ratio(string fasta){
	auto data = NamedVector<char>::getNamedVector(fasta, "data");
	auto sorted = NamedVector<uint32_t>::getNamedVector(fasta, "sorted");
	auto orgstart = NamedVector<uint32_t>::getNamedVector(fasta, "orgstart");
	auto string_length = NamedVector<uint16_t>::getNamedVector(fasta, "strlen");

	vector<map<int,int>> total_length_count;
	ofstream f_org_str, f_para;
	string dir = fasta.substr(0, fasta.find_last_of("\\/")) + "/";
	umask(0);
	mkdir((dir+"results").c_str(), 0700);
	f_org_str.open(dir+"results/org_length_count.csv", ios::out);
	f_para.open(dir+"results/lf_count_exp_para.csv", ios::out);

	uint32_t i,stop;
	size_t string_size;
	size_t m=0;
	uint32_t maxO = orgstart->size()-1;
	cout << "org length count" << endl;
	f_org_str << "org length count" << endl;
	for (uint32_t o=0; o<orgstart->size(); ++o){
		i = orgstart->at(o);
		stop = (o<orgstart->size()-1)?(orgstart->at(o+1)):data->size();
		cerr << "o=" << o << " start=" << i << " stop=" << stop << endl;
		string_size = 0;
		size_t mi=0;
		map<int,int> length_count;
		for (;;){
			if (data->at(i) == '\0'){ 
				length_count[string_size] ++;
				mi = max(mi, string_size);
				string_size = 0;
			}
			if (i == stop) {
				cout << "STOP !" << string(&data->at(i)) << endl;
				if (stop == data->size()){
					length_count[string_size] ++;
					mi = max(mi, string_size);
					string_size = 0;
				}
				map<int,int> length_count_new;
				m = max(m,mi);
				total_length_count.push_back(map<int,int>());
				for (int j=mi; j>0; --j){
					for (int k=j; k>0; --k){
						length_count_new[k] += length_count[j]*(j-k+1);
						total_length_count[o][k] += length_count[j]*(j-k+1);
					}
				}
				for (auto val:length_count_new) {
					//cout << o << " " << val.first << " " << val.second << endl;
					f_org_str << o << " " << val.first << " " << val.second << endl;
				}
				break;
			}
			++string_size;
			++i;
		}
		if (o == maxO) break;
	}
	cerr << "lf count p_para p_ortho" << endl;
	f_para << "lf count p_para p_ortho ortho_para" << endl;
	for (uint32_t l=0; l<m; ++l){
		long double P_para = 0.0;
		long double T = 0.0;
		for (uint32_t o=0; o<maxO; ++o)	T += (total_length_count[o][l]);
		for (uint32_t o=0; o<maxO; ++o){
			if (total_length_count[o][l] < 0) cerr << "[ERROR] o=" << o << " l=" << l <<" total_length_count[o][l]=" << total_length_count[o][l] << endl;
			P_para += (total_length_count[o][l]/T) * (total_length_count[o][l]/T);
		}
		long double ortho = 1.0;
		ortho -= P_para;
		f_para << l << " " << T << " "  << P_para << " " << ortho << " " << ortho/P_para << endl;
	}

	cout << "sorted to: " << dir+"results/lf_count_exp_para.csv" << endl;
}


vector<map<char, double>> protwise_AA_prob(data_t data, indices_t di1, indices_t di2, size_t max_lf){
	vector<map<char, double>> aa_count(max_lf+1);
	indices_t start = di1;
	indices_t s = di2-start;
	for (uint16_t p=0; p<s-1; ++p){
		uint16_t q = s-p;
		char c = data->at(start+p);
		if (c == 'X') cerr << "[ERROR] encountered a wild X" << endl;
		char state = 0;
		size_t lf=1;
		while (lf <= min(uint16_t(p+1),q) && lf <= max_lf){
			aa_count[lf][c] += lf;
			++lf;
		}
		uint16_t stat = lf-1;
		while (lf <= max(uint16_t(p+1), q) && lf <= max_lf){
			aa_count[lf][c] += stat;
			++lf;
		}
		while (lf <= max_lf){
			stat--;
			aa_count[lf][c] += stat;
			++lf;
		}
	}
	for (int lf = 1; lf<=max_lf; ++lf) {
		double sum = 0.0;
		for (auto v:aa_count[lf]) sum += double(v.second);
		for (auto v:aa_count[lf]) aa_count[lf][v.first] /= sum;
	}
	return aa_count;
}

void calc_lf_AAprob_S(string fasta){
	auto data = NamedVector<char>::getNamedVector(fasta, "data");
	uint16_t max_lf = 300;
	vector<map<char, uint64_t>> aa_count;
	for (uint16_t j=0; j<max_lf+1; ++j)  aa_count.push_back(map<char, uint64_t>());
	indices_t start = 0;
	for (indices_t i=0; i<data->size(); ++i){
		start = i;
		while (data->at(i) != '\0') ++i;
		uint16_t s = i-start;
		//cout << "s=" << s << " start=" << start << " i=" << i << endl; 
		//cout << string(&data->at(start)) << endl;
		for (uint16_t p=0; p<s-1; ++p){
			uint16_t q = s-p;
			char c = data->at(start+p);
			if (c == 'X') cerr << "[ERROR] encountered a wild X" << endl;
			char state = 0;
			//cout << "  p=" << p << " q=" << q << " c=" << c << endl;
			uint16_t lf=1;
			while (lf <= min(uint16_t(p+1),q) && lf <= max_lf){
				aa_count[lf][c] += lf;
				++lf;
			}
			uint16_t stat = lf-1;
			while (lf <= max(uint16_t(p+1), q) && lf <= max_lf){
				aa_count[lf][c] += stat;
				++lf;
			}
			while (lf <= max_lf){
				stat--;
				aa_count[lf][c] += stat;
				++lf;
			}
		}
	}
	cout << "lf"; for (auto v:aa_count[1]) cout << " " << v.first; cout << " repick_prob" << endl;
	for (int lf = 1; lf<aa_count.size(); ++lf) {
		cout << lf << ")";
		double sum = 0.0;
		for (auto v:aa_count[lf]) sum += double(v.second);

		double mu = 0.0;
		for (auto v:aa_count[lf]) {
			cout << " "  << double(v.second)/sum;
			mu += pow(double(double(v.second)/sum), 2.0);
		} cout << " " << mu << endl;
	}

}



void calc_prot_AA(string fasta){
	auto data = NamedVector<char>::getNamedVector(fasta, "data");
	auto protstart = NamedVector<uint32_t>::getNamedVector(fasta, "protstart");
	auto protcount = NamedVector<uint32_t>::getNamedVector(fasta, "fasta_prot_count");

	//setting::valid_chars[mod];
	//setting::A;;
	
	size_t max_lf = 300;
	vector<vector<long double>> PTc; // PTC[lf][char]
	vector<long double> T;
	vector<long double> v; for (char c='A'; c<='Y'; ++c) v.push_back(0.0);
	for (size_t lf = 0; lf <= max_lf; ++lf) {
		PTc.push_back(v);
		T.push_back(0.0);
	}

	// calc PTc and T
	for (uint32_t pi=1; pi<protstart->size(); ++pi){
		long double lp = protstart->at(pi)-protstart->at(pi-1)-1;
		vector<map<char, double>> aa_count = protwise_AA_prob(data, protstart->at(pi-1), protstart->at(pi), max_lf);
		if (aa_count[1][0] != 0 || aa_count[1]['X'] != 0) continue;

		for (size_t lf = 1; lf <= max_lf && lf <=lp; ++lf){
			long double w = lp-lf+1;
			T[lf] += w;
			for (char c='A'; c<='Y'; ++c) PTc[lf][c-'A'] += aa_count[lf][c] * w;
		}
		if (pi % (protstart->size()/1000) == 0) cerr << pi << "  calc" << endl;
	}

	// print PTc
	for (size_t lf = 0; lf <= max_lf; ++lf) {
		cout << "#" << lf << ")";
		long double sum=0.0;
		for (char c='A'; c<='Y'; ++c) {
			PTc[lf][c-'A'] /= T[lf];
			sum += PTc[lf][c-'A'];
			cout << " " << PTc[lf][c-'A'];
			if (PTc[lf][c-'A']>1 || PTc[lf][c-'A']<0) return;
		} cout << endl; //" | "<< sum << endl;
	}

	map<uint32_t, double> bino = calc_binom(max_lf);

	vector<vector<long double>> Plf_d; //Plf_d[lf][d<=lf]
	for (size_t lf = 0; lf <= max_lf; ++lf) { Plf_d.push_back(vector<long double>(lf+1)); }

	cerr << " Now calc Plf_d" << endl;
	for (uint32_t pi=1; pi<protstart->size(); ++pi){
		vector<map<char, double>> aa_count = protwise_AA_prob(data, protstart->at(pi-1), protstart->at(pi), max_lf);
		if (aa_count[1][0] != 0 || aa_count[1]['X'] != 0) continue;
		long double lp = protstart->at(pi)-protstart->at(pi-1)-1;

		for (size_t lf = 0; lf <= max_lf && lf <=lp; ++lf) {
			long double w = double(lp-lf+1)/T[lf];
			long double mu  = 0.0;
			for (char c='A'; c<='Y'; ++c) mu += aa_count[lf][c] * PTc[lf][c-'A'];
			for (size_t d=0.0; d<=lf; ++d){
				Plf_d[lf][d] += bino[(lf<<16)+(lf-d)] * pow(mu, lf-d) * pow(1.0-mu, d) * w;
			}
		}
		if (pi % (protstart->size()/1000) == 0) cerr << pi << endl;
	}

	cout << "# plf_d" << endl << "lf";
	for (size_t lf = 0; lf <= max_lf; ++lf) { cout << " " << lf; } cout << endl;
	for (size_t lf = 1; lf <= max_lf; ++lf) {
		long double sum = 0.0;
		cout << lf;
		for (size_t d=0.0; d<=lf; ++d) {cout << " " << Plf_d[lf][d]; sum += Plf_d[lf][d]; }
		for (size_t d=lf+1; d<=max_lf; ++d) cout << " 0";
		cout << endl;
	}
}



map<uint32_t, double> calc_binom(size_t max_lf){
	map<uint32_t, double> ret;
	for (size_t i=0; i<=max_lf; ++i){
		for (size_t j=i; j<=max_lf; ++j){
			double b = binom(j,i);
			ret[(i<<16)+j] = b;
			ret[(j<<16)+i] = b;
		}
	}
	return ret;
}

















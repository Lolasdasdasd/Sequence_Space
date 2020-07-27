#include <iostream> 
#include <string>
#include <vector>
#include <map>
#include <set>
#include <stdlib.h> 
#include <sstream>
#include <iostream>
#include <fstream>
#include <algorithm>

using namespace std;

set<string> assi_list_to_set(string s1);
bool has_family_overlap(string s1, string s2);
bool has_potential_homology_overlap(string s1, string s2);

int main(int argc, char** argv){
	if (argc != 3){
		cout << "please provide <.assign> <.non_homologs>" << endl;
		return -1;
	}
	string assign_filename = string(argv[1]);
	string nonhomo_filename = string(argv[2]);

	vector<map<int,string>> assign_ecod;
	for (int i=0; i<=9; ++i) assign_ecod.push_back(map<int,string>());


//  assign
//  0 0 314.1.1.2
//  0 1000000000 208.1.1.9
//  0 10000000
//  0 1001000000 2003.1.9.16
//  0 1002000000 12.3.1.18

    ifstream f_assign, f_nonhomo;
	string line, ECOD;
	int VAR, POS;
	f_assign.open(assign_filename);
	while (std::getline(f_assign, line)) {
			VAR=0; POS=0; ECOD="";
			istringstream(line) >> VAR >> POS >> ECOD;
			//cout << "var=" << VAR << " pos=" << POS << " ecod=" << ECOD << endl;
			if (assign_ecod[VAR].find(POS) != assign_ecod[VAR].end()) {
				assign_ecod[VAR][POS] = assign_ecod[VAR][POS] + "-" + ECOD;
			} else {
				assign_ecod[VAR][POS] = ECOD;
				//cout << "else" << endl;
			}
	}
	f_assign.close();

	f_nonhomo.open(nonhomo_filename);

// # non_homologs
//	0 1 0 28010000 0.11


	bool PRINT=true;

	double overlaps=0, non_overlaps=0, question_mark=0;

	int VAR1, VAR2, POS1, POS2;
	double PROB;
	while (std::getline(f_nonhomo, line)) {
		istringstream(line) >> VAR1 >> VAR2 >> POS1 >> POS2 >> PROB;
		//cout << "var1=" << VAR1 << " var2=" << VAR2 << " pos1=" << POS1 << " pos2=" << POS2 << endl;
		string assi1="", assi2="";
		if (assign_ecod[VAR1].find(POS1) != assign_ecod[VAR1].end()) {
			assi1 = assign_ecod[VAR1][POS1];
			if (assi1.size() == 0) {
				question_mark++; // frag1 not assigned
				if (PRINT) cout << line << " Q" << endl;
				continue;
			}
			if (assign_ecod[VAR2].find(POS2) != assign_ecod[VAR2].end()) {
				// both fragments have at least 1 assigned domain
				assi2 = assign_ecod[VAR2][POS2];
				if (assi2.size() == 0) {
					question_mark++; // frag2 not assigned
					if (PRINT) cout << line << " Q" << endl;
					continue;
				}
				bool overlap = has_potential_homology_overlap(assi1, assi2);
				if (overlap){
					overlaps++;
					if (PRINT) cout << line << " P" << endl;
				} else {
					if (PRINT) cout << line << " C" << endl;
					//cout << line << "  ph=" << overlaps<< " c=" << non_overlaps << " ?=" << question_mark << " " << assi1 << " / " << assi2 << endl;
					non_overlaps++;
				}
        	}
		}
	}
	double T=overlaps+non_overlaps+question_mark;
	cerr << "Questionmark=" << 100.0*question_mark/T << " PotentialHomology=" << 100.0*overlaps/T << " Convergent=" << 100.0*non_overlaps/T << endl;
	f_nonhomo.close();

}

// on class X
bool has_potential_homology_overlap(string s1, string s2){
	set<string> assis1 = assi_list_to_set(s1);
	set<string> assis2 = assi_list_to_set(s2);

	for (auto a1:assis1){
		a1 = a1.substr(0, a1.find(".") + 1);
		for (auto a2:assis2){
			a2 = a2.substr(0, a2.find(".") + 1);
			// 2003.1.2.282 4021.1.1.2
			if (a1.compare(a2) == 0){
				//cout << s1 << " OVERLAP WITH " << s2 << endl;
				return true;
			}
		}
	}
	return false;
}


bool has_family_overlap(string s1, string s2){
	set<string> assis1 = assi_list_to_set(s1);
	set<string> assis2 = assi_list_to_set(s2);

	for (auto a1:assis1){
		for (auto a2:assis2){
			// 2003.1.2.282 4021.1.1.2
			if (a1.compare(a2) == 0){
				//cout << s1 << " OVERLAP WITH " << s2 << endl;
				return true;
			}
		}
	}
	return false;
}


set<string> assi_list_to_set(string s1){
	set<string> ret;
	string delimiter = "-";
	size_t pos = 0;
	string token;
	while ((pos = s1.find(delimiter)) != std::string::npos) {
	    token = s1.substr(0, pos);
	    //std::cout << token << std::endl;
		ret.insert(s1);
	    s1.erase(0, pos + delimiter.length());
	}
	//std::cout << s1 << std::endl;
	ret.insert(s1);
	return ret;

}





















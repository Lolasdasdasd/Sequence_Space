#include "../CachedVector.h"
#include "../string_comparison.h"
#include "../settings.h"

#include <algorithm>
#include <cstdlib>

using namespace std;
char equal_char();

int main(int argc, char** argv){
	if (argc != 4 && argc != 3){
		cout << "please provide: \n<full path! Original project> \n<protein:P / fasta:F / all:A / equal:E> \n[random seed default=0]" << endl;
		exit(0);
	}
	string old_project = string(argv[1]);
	char shuffle_mode = argv[2][0];

	string new_project = old_project;
	switch (shuffle_mode){
		case 'P': { new_project += ".protshuf"; break; }
		case 'F': { new_project += ".fastashuf"; break; }
		case 'A': { new_project += ".allshufX"; break; }
		case 'E': { new_project += ".equal"; break; }
		cout << "[ERROR] wrong shuffle mode: eighter P or F" << endl;
		exit(0);
	}
	int seed;
	if (argc==4) seed = atoi(argv[3]);
	else seed = 0;
	srand(seed);
	new_project += "_" + to_string(seed);

	string old_fasta = old_project + "/_fastafile";
	string new_fasta = new_project + "/_fastafile";

	int ret = 0;
	ret += system(("rm -rf " + new_project).c_str());
	ret += system(("mkdir " + new_project).c_str());
	cout << ("ln -s " + old_fasta + " " + new_fasta) << endl;
	ret += system(("ln -s " + old_fasta + " " + new_fasta).c_str());
	cout << ("ln -s " + old_project + "/_fasta_prot_count " + new_fasta + "/_fasta_prot_count") << endl;
	ret += system(("ln -s " + old_project + "/_fasta_prot_count " + new_project + "/_fasta_prot_count").c_str());
	cout << ("ln -s " + old_project + "/protstart " + new_fasta + "/protstart") << endl;
	ret += system(("ln -s " + old_project + "/protstart " + new_project + "/protstart").c_str());

	
	{
		auto data = NamedVector<char>::getNamedVector(old_fasta, "data");
		auto protstart = NamedVector<uint32_t>::getNamedVector(old_fasta, "protstart");
		auto protcount = NamedVector<uint32_t>::getNamedVector(old_fasta, "fasta_prot_count");
		// indices in data with characters
		WriteCachedVector<uint32_t> tmp(NamedVectorBase::cacheFilename(new_fasta, "data.tmp"));
		tmp.reserve(data->size());

		size_t protein = 0;
		size_t organism = 0;
		vector<indices_t> boundaries;
		indices_t valids = 0; // used in E-mode
		if (shuffle_mode == 'F' || shuffle_mode == 'P'){
			if (shuffle_mode == 'F') boundaries.push_back(0);
			// #boundaries = #protein because protein contains last boundary as well
			for (indices_t i=0; i<data->size(); ++i){
				if (protstart->at(protein) == i){
					if (shuffle_mode == 'P') boundaries.push_back(tmp.size());
					else if (shuffle_mode == 'F' && protcount->at(organism) == protein){
						boundaries.push_back(tmp.size());
						++organism;
					}
					++protein;
				}
				if (data->at(i) && data->at(i) != 'X') tmp.push_back(i); // added data->at(i) != 'X' on 20.10.17
			}
		} else if (shuffle_mode == 'A'){ // shuffle case A = all or E = equal
			for (indices_t i=0; i<data->size(); ++i){
				if (data->at(i) && data->at(i) != 'X')
					tmp.push_back(i);
			}
			boundaries.push_back(0);
			boundaries.push_back(tmp.size());
		} else if (shuffle_mode == 'E'){ // shuffle case A = all or E = equal
			// find all characters in data
			vector<indices_t> pos;
			set<char> found_char;
			indices_t V = setting::valid_chars.size();
			for (indices_t i=0; found_char.size() != V; ++i){
				if (data->at(i) && found_char.find(data->at(i)) == found_char.end()){
					found_char.insert(data->at(i));
					pos.push_back(i);
				}
			}
			// fill tmp with equally many letters of each type
			for (indices_t i=0; i<data->size(); ++i){
				if (data->at(i) && data->at(i) != 'X') {
					tmp.push_back(pos[valids%V]);
					valids++;
				}
			}
			boundaries.push_back(0);
			boundaries.push_back(tmp.size());
		}

		// MODIFY TMP
		cout << "# protein = " << protein << " protstart->size()="<< protstart->size() <<" / # boundaries="<< boundaries.size() <<endl;
		cout << "tmp filled -> now shuffle!..." << endl;
		for (size_t i=0; i<boundaries.size()-1; ++i){ // protstart contains end of last prot as well
			if (shuffle_mode == 'F' && i % (organism/100) == 0) cout << "shuffle fasta " << i << endl;
			if (shuffle_mode == 'P' && i % (protein/100) == 0) cout << "shuffle protein " << i << endl;
			if (shuffle_mode == 'A') cout << "shuffle all tmp.size()" << tmp.size() << " =? " << boundaries[1]<< " (shuffle boundary)"<< endl;

			if (shuffle_mode == 'F' || shuffle_mode == 'P')	random_shuffle(tmp.begin()+boundaries[i], tmp.begin()+boundaries[i+1]);
			else if (shuffle_mode == 'A' || shuffle_mode == 'E'){
				indices_t save, randi, N=tmp.size();
				cout << "RAND_MAX = " << RAND_MAX << " >? " << N << " div="<< float(RAND_MAX)/float(N)<<  " is great number or close to an int? "<< endl;
				for (indices_t i=0; i<N; ++i){
					randi= rand() % N;
					save = tmp[i];
					tmp[i] = tmp[randi];
					tmp[randi] = save;
					if (i % (N/1000) == 0) cout << "shuffle index " << i << endl;
				}
			}
		}
		cout << "all shuffled it is (:" << endl;

		// ### write shuffled data file! ###
		WriteCachedVector<char> new_data(NamedVectorBase::cacheFilename(new_fasta, "data"));
		new_data.reserve(data->size());

		cout << "now start writing to random data..." << endl;
		uint32_t pos = 0;
		for (uint32_t i=0; i<data->size(); ++i){
			if (data->at(i)) {
				if (data->at(i) != 'X') {
					new_data.push_back(data->at(tmp.at(pos)));
					++pos; // change from 14.11.
				} else new_data.push_back('X');
			} else {
				new_data.push_back(0);
			}
		}
		cout << "### old data to randomize ###" << endl;
		for (int i=0; i<1000; ++i) {
			cout << data->at(i);
			if (!data->at(i)) cout << endl;
		}
		cout << "### new_data start like this in WCV ###" << endl;
		for (int i=0; i<1000; ++i) {
			cout << new_data.at(i);
			if (!new_data.at(i)) cout << endl;
		}

	}
	ret += system(("rm -f " + NamedVectorBase::cacheFilename(new_fasta, "data.tmp")).c_str());

	auto read_data = NamedVector<char>::getNamedVector(new_fasta, "data");

	cout << "### open writen data as sharted_ptr ###" << endl;
	for (int i=0; i<1000; ++i) {
		cout << read_data->at(i);
		if (!read_data->at(i)) cout << endl;
	}

	return 0;
}



char equal_char(){
	return setting::valid_chars[rand() % setting::valid_chars.size()];
}

#include "../CachedVector.h"
#include "../utilities.h"
#include "../settings.h"
#include <fstream>
#include <memory> // make_shared

using namespace std;


class CompressedData : public NamedVector<char> {
public:
	CompressedData() : NamedVector<char>("data"){};

protected:
	shared_ptr<CachedVector<char>> generateVector(string fasta, int param) {
		string filename = cacheFilename(fasta,"data");
		string tmp = filename + ".tmp";

		std::ifstream infile(filename);
		if (infile.good()) { // changed on 22.3.2019
			infile.close();
			cerr << "  ### existing data ###" << endl;
			return make_shared<CachedVector<char>>(filename);	
		} 
		{
			std::cerr << filename << std::endl;

			WriteCachedVector<char> v(tmp);

			ifstream fin(fasta, fstream::in);
			ofstream flog;
			flog.open(cacheFilename(fasta,"_data.log"));
			flog << "setting::valid_chars=" << setting::valid_chars << endl;
			flog << "setting::exchange_chars=" << setting::exchange_chars << endl;

			char ch;
			int n_protein = -1, state = 0, length = 0, invalids = 0, valids = 0;
			map<int,int> protlength;
			map<char,int> AAcount;
			map<char,int> inva;
			string protein;
			// ### READ ###	   IN
			while (fin >> noskipws >> ch) {
				//if (valids+invalids+n_protein+1 == 10*1024*1024){ v.push_back('\0'); cerr << "prebreak 10M" << endl; break; }
				if (state == 0 && ch == '>') { // starting header
					state = 1; n_protein++;
					if (n_protein > 0) { 
						protlength[length]++; 
						v.push_back('\0'); 
					}
					length = 0;
				} else if (state == 1 && ch != '>'){ // reading header
					protein += ch;
					if (ch == '\n') {
						state = 0;
					}
				} else if (state == 0 && ch != '\n'){ // reading sequence
					++length; ++AAcount[ch];
					size_t AA = setting::valid_chars.find(ch);
					if (AA != string::npos) {
						v.push_back(ch); AAcount[ch]++; valids++;
					} else {
						if (ch >= 'A' && ch<='Z') v.push_back(setting::exchange_chars[ch-'A']);
						else v.push_back('X'); 
						invalids++; inva[ch]++;
					}
				}
			}
			v.push_back('\0');
			++n_protein;
			flog << "n_protein=" << n_protein << endl;
			flog << "n_chars=" << v.size()-n_protein << endl;
			flog.close();
			fin.close();
			map2file(&inva, "invalid,count", cacheFilename(fasta,"_invalids.map"));
			map2file(&AAcount, "AA,count", cacheFilename(fasta,"_AAcount.map"));
			map2file(&protlength, "protein_length,count", cacheFilename(fasta,"_protein_length.map"));
		}
			rename(tmp.c_str(),filename.c_str());
			return make_shared<CachedVector<char>>(filename);	
		}

} registerClass_data;



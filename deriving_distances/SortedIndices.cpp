#include "../CachedVector.h"
#include <fstream>
#include <algorithm>

using namespace std;

#define TYPE uint32_t
#define NAME "sorted"

class SortedIndices : public NamedVector<TYPE> {
public:
	SortedIndices() : NamedVector<TYPE>("sorted") {};

protected:
	shared_ptr<CachedVector<TYPE>> generateVector(string fasta, int param) {
		return generateVector(fasta);
	}

	shared_ptr<CachedVector<TYPE>> generateVector(string fasta) {

		string filename = cacheFilename(fasta,"sorted");
		string tmp = filename + ".tmp";
		if (fileNewer(filename,cacheFilename(fasta,"data"))) {
			cerr << "  ### existing sorted ###" << endl;
			return make_shared<CachedVector<TYPE>>(filename);
		} 
		// dependencies
		shared_ptr<CachedVector<char>> data = NamedVector<char>::getNamedVector(fasta, "data");
		{
			cerr << "  ### GENERATE sorted ###" << endl;
			WriteCachedVector<TYPE> v(tmp);
			size_t D = data->size();
			v.reserve(D);
			// fill
			for (TYPE i=0; i<data->size(); ++i){
				if (data->at(i)) v.push_back(i);
			}
			// sort
			stable_sort(v.begin(), v.end(), [&](TYPE v1, TYPE v2) {
				int off=0;
				while (data->at(v1+off) == data->at(v2+off)
					&& data->at(v1+off) && data->at(v2+off) 
								&& v1+off < D   && v2+off < D) { off++; }
				return (data->at(v1+off) < data->at(v2+off));
			});
		}
			rename(tmp.c_str(),filename.c_str());
			return make_shared<CachedVector<TYPE>>(filename);
		}

} registerClass_sorted;



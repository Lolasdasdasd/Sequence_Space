#include <string>
#include <map>

namespace setting
{

	size_t A = 20;
	std::string valid_chars   ="ACDEFGHIKLMNPQRSTVWY";
	std::string all_chars     ="ABCDEFGHIJKLMNOPQRSTUVWXYZ";
	std::string exchange_chars="ADCDEFGHIXKLMNKPQRSTCVWXYE";
	std::map<char,size_t> char2int {{'A',0},{'C',1},{'D',2},{'E',3},{'F',4},{'G',5},{'H',6},{'I',7},{'K',8},{'L',9},{'M',10},{'N',11},{'P',12},{'Q',13},{'R',14},{'S',15},{'T',16},{'V',17},{'W',18},{'Y',19}};



}

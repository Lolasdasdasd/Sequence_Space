#include <stdio.h>      /* printf */
#include <math.h>       /* pow */
#include <string.h>
#include <unistd.h>
#include <dirent.h>
#include <sys/types.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <set>
#include <map>
#include <vector>
#include <cstdlib>


#define COLOR_RESET   "\x1b[0m"
#define COLOR_RED     "\x1b[31m"
#define COLOR_GREEN   "\x1b[32m"
#define COLOR_YELLOW  "\x1b[33m"
#define COLOR_BLUE    "\x1b[34m"
#define COLOR_MAGENTA "\x1b[35m"
#define COLOR_CYAN    "\x1b[36m"

#define BG_COLOR_RESET   "\e[0m"
#define BG_COLOR_BLACK   "\e[40m"
#define BG_COLOR_RED     "\e[41m"
#define BG_COLOR_GREEN   "\e[42m"
#define BG_COLOR_YELLOW  "\e[43m"
#define BG_COLOR_BLUE    "\e[44m"
#define BG_COLOR_MAGENTA "\e[45m"
#define BG_COLOR_CYAN    "\e[46m"
#define BG_COLOR_WHITE   "\e[47m"

#define PREV_LINE        "\033[F"

using namespace std;

double binom(size_t n, size_t k){
	if (n < k) cerr << "[ERROR] in binom: n<k   :  "<< n << " <" << k << endl;
	double ret = 1.0;
	if (n == k) return 1.0;
	for (double i=1.0; i<=double(n-k); ++i) ret *= double(n-i+1) / i;
	return ret;
}

double binom(int n, int k){
	double ret = 1.0;
	for (double i=1.0; i<=n-k; ++i) ret *= (n-i+1) / i;
	return ret;
}

vector<string> files_in_dir(string path){
	if (path[path.size()-1] != '/') path += "/";
	DIR *dir;
	struct dirent *ent;
	vector<string> files;
	if ((dir = opendir(path.c_str())) != NULL) {
		while ((ent = readdir(dir)) != NULL) {
			if (ent->d_name[0] != '.'){
				files.push_back(path + ent->d_name);
			}
		}
		closedir(dir);
	} else { 
		cout << "directory " << path << " could not be read..." << endl; 
	}
	return files;
}

string fasta_from_file(string fasta){
	ifstream in(fasta, ios::in);
	string line = "", content = "";
	if (in.is_open()){
		while (in.good()){
			getline(in,line);
			if (line[0]== '>') continue;
			content += line;
		}
	}
	return content;
}

string file_content(string filename, int nrows){
	ifstream in(filename, ios::in);
	string line = "", content = "";
	int n=0;
	if (in.is_open() && n<nrows){
		while (in.good()){
			getline(in,line);
			content += line;
			++n;
		}
	}
	return content;
}

int get_file_size(string filename){
    FILE *p_file = NULL;
    p_file = fopen(filename.c_str(),"rb");
    if (!p_file){
        cout << "[ERROR] utilities::get_file_size fails, no file: " << filename << endl;
        return -1;
    }
    fseek(p_file,0,SEEK_END);
    int size = ftell(p_file);
    fclose(p_file);
    return size;
}

bool file_exists(string file_name){
    ifstream ffile(file_name);
    if (ffile){
        ffile.close();
        return true;
    }
    return false;
}

void timestamp(char* buffer){
    time_t now = time(NULL);
    strftime(buffer, 200, "%Y-%m-%d %H:%M:%S", localtime(&now));
}

void progress_bar_cerr(unsigned long current, unsigned long aim){
	cerr << "\r             \r" << double(current)/double(aim)*100.0 << " % ";
}

void progress_bar(unsigned long current, unsigned long aim){
    printf(COLOR_GREEN "\r%5.1f%% " COLOR_RED "[", current*100.0/aim);

    for (int i=1; i<current*100.0/aim; i++){
        printf("=");
    }
    printf(">");
    for (int i=current*100.0/aim; i<100; i++){
        printf(" ");
    }
    printf("]" COLOR_GREEN" %5.1f%%" COLOR_RESET, 100.0-(current*100.0/aim));
    fflush(stdout);
}

double string_consensus(char const* s1, char const* s2){
	int o = 0;
	for (size_t i=0; i<strlen(s1); ++i){
		for (size_t j=0; j<strlen(s2); ++j){
			//cout << string(s1+i) << " " << string(s2+j) << "    ";
			for (int k=0; ;++k){
				if (s1+i+k && s2+j+k && *(s1+i+k) == *(s2+j+k)) { ++o; } //cout << "+"; }
				else break;
			}
			//cout << endl;
		}
	}
	return double(o/100.0);
}


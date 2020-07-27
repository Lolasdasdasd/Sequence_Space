#include "CachedVector.h"
#include "utilities.h"

#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>

#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>

using namespace std;

typedef pair<uint32_t, uint32_t> Bounds;
void ANA_lf_recur_dup(string fasta);
void ANA_lf_abu_frags(string fasta, size_t lf);
void ANA_frag_intra_inter(string fasta, size_t lf, int count);
void ANA_frag_exp_abu(string fasta, size_t lf, int N, int take_every=1);
void error_to_exp(string fasta, size_t lf);
void print_strings_longer_then(string fasta, size_t lf);
void exp_ortho_para_ratio(string fasta);
void calc_lf_AAprob_S(string fasta);
void calc_prot_AA(string fasta);
int compositional_overlap(indices_t di1, indices_t di2, data_t data, size_t lf);
double compo_exp(string str1, string str2);
size_t cal_transition(size_t lf, double precision, double match);
pair<size_t, size_t> extend_to_twi_zone(size_t lf, size_t transition, size_t radius);

map<char, double> get_AA_prob(data_t data);
map<char, double> get_AA_prob_file(string fasta);
void word_count_analysis(string fasta, size_t lf, double factor);
double frag_prob(const map<char, double>* prob, string frag);
void gen_frag(size_t f, string* frag);
indices_t gen_fragpos(string frag);
indices_t valid_hits_of_length(data_t data, size_t lf);
void existing_neighbors(string fasta, int lf, string frag, string prof_file="", float overlap=0.9);
pair<double,double> ortholog_paralog(string fasta, string frag, data_t data, sorted_t sorted, sorted_t orgstart);
vector<map<char, double>> protwise_AA_prob(data_t data, indices_t di1, indices_t di2, size_t max_lf);
map<uint32_t, double> calc_binom(size_t max_lf);
void frag_shuf_model(string fasta, size_t max_lf, uint64_t start);
void frag_shuf_ortho_comp(string fasta, size_t lf, int start);
void frag_shuf_para_comp(string fasta, size_t lf, int start, int org, int locality_PM=-1); // old
void frag_shuf_para_comp(string fasta, size_t lf, int org); // new no random seed...

int produce_maximal_minimal_unique_string_length(string fasta);
int produce_hit_counts(string fasta, int lf);

void intra_distance(string fasta, int lf, indices_t start, indices_t stop, int step=1);
void inter_distance(string fasta, string fasta_rand, int lf, indices_t start, indices_t stop, int step=1, string name="");
void heuristic_inter_distance(string fasta, size_t lf, size_t rep, size_t precision=9);
void heuristic_inter_distance_twilight(string fasta, size_t lf, int min_twilight, int max_twilight);

void connected_component_max_PM(string fasta, indices_t di1, indices_t di2, size_t lf, int max_PM);
void cluster_string_max_PM(string fasta, string frag, int max_PM);
void cluster_max_PM(string fasta, size_t lf, int max_PM);
void exhaustive_inter_distance(string fasta, size_t lf);
void all_conditional_entropy(string fasta, size_t hits=2, bool drop_X=true);
void conditional_entropy(string fasta, size_t lf, data_t data, sorted_t sorted, prelen_t common_prefix, strlen_t string_length, size_t min_hits=2, size_t max_hits=2, bool drop_X=true);
void locality_in_proteins(string fasta);
void locality_in_protein(string prot, bool randomize=false);

void subsequent_AAs(data_t data, sorted_t sorted, string fragment, map<char, uint32_t>* entropical_end);
void predict_subsequent_AAs(string fasta, string protein, size_t lf);
char predict_subsequent_AA(data_t data, sorted_t sorted, string fragment);

void print_protein_with_fragment(string fasta, string frag, bool print_prot=true, bool pdb=false);

indices_t string2df(string& f, data_t data, diverse_frags_t diverse_frags, CachedVector<indices_t>::iterator b, CachedVector<indices_t>::iterator e);
bool contains_char(data_t data, indices_t di, size_t lf, char ch);

// returns (proteinid, position in protein)

set<uint32_t> orgs_with_frag(string fasta, data_t data, sorted_t sorted, string fragment, orgstart_t orgstart, bool print_orgs);

pair<uint32_t, uint32_t> data2protpos(protstart_t protstart, uint32_t pos, protcount_t protcount);
uint32_t protpos2organism(pair<uint32_t, uint32_t> p);
uint32_t protpos2prot(pair<uint32_t, uint32_t> p);
uint32_t protpos2pos(pair<uint32_t, uint32_t> p);
string pretty_protein_name(string fasta, uint32_t protein, ifstream* fprotein, protein_names_t pn);
string full_protein_name(string fasta, uint32_t protein, ifstream* fprotein, protein_names_t pn);
string organism_name(string fasta, uint32_t org);
string protein_at(uint32_t di, data_t data, sorted_t sorted, protstart_t protstart, protcount_t protcount);
string protein(string fasta, uint32_t protein_number);
void print_protein(string fasta, uint32_t protein_number);

void intra_prot_aa_dist(string fasta, string filename);
void intra_prot_aa_dist(uint32_t organism, string fasta, string filename);
void intra_organism_aa_dist(string fasta, string filename);
map<char, uint32_t> amino_acid_distribution(string fasta, indices_t di_start, indices_t di_stop, data_t data);
// to investigate weird bump
void filter_and_where_from(string fasta, Bounds b_fraglen, Bounds b_hits, Bounds b_orgs, bool print_all=true);

// for a given fragment length calculates how often a fragment of this size occurs
// returns hits->#frags_with_thismany_hits
map<uint32_t,uint32_t> hit_counts(string fasta, size_t lenfrag);
map<uint32_t,uint32_t> hit_counts_per_protein(string fasta, size_t lenfrag);

// ### HITS #########################################
// returns how often the fragment occurs in inputfile
int hits(const string fasta, const string fragment);
int hits(data_t data, sorted_t sorted, const string fragment);
void print_hits_exp_hits(string fasta, string frag);
// given the index in sorted, and fraglen, how many hits?
uint32_t hits(uint32_t si, size_t fraglen, prelen_t common_prefix, strlen_t string_length);
// returns hits_per_protein, total hits
pair<uint32_t,uint32_t> hits_per_prot(uint32_t si, size_t fraglen, sorted_t sorted, prelen_t common_prefix, strlen_t string_length, protstart_t protstart, protcount_t protcount);
map<uint32_t, uint32_t> frags_per_org(string fasta, data_t data, sorted_t sorted, string fragment, orgstart_t orgstart, bool print_orgs);

// ### MAPS, STATISTICS ############################

// returns most often occuring fragments of certain fraglen
vector<string> most_often_occ(string fasta, int fraglen, int min_max_count=5);
// strings with following properties:
// * occur once
// * occur at least 2x if cut off end
map<uint32_t,uint64_t> minimal_unique_strings(string fasta);
// strings with following properties:
// * occur once
// * occur at least 2x if cut off front or end
// * that are not substings of other strings with properties 1+2
map<uint32_t,uint64_t> maximal_minimal_unique_string_length(string fasta);

// ### SEARCH ###

// sorted->at(si) == fragment
// returns si
uint32_t string2sorted(string& fragment, data_t data, sorted_t sorted);

// if a unique string starts at data->at(sorted->at(si))
// returns the length of the minimal unique string, else returns 0
size_t is_unique_string_of_size(uint32_t si, data_t data, sorted_t sorted, prelen_t common_prefix);

// prints all vectors
void print_vectors(string fasta);

double string_consensus(char const* s1, char const* s2);

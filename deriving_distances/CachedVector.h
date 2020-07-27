#ifndef BASICCACHEDVECTOR_H
#define BASICCACHEDVECTOR_H

#include <unistd.h> // getpagesize
#include <sys/mman.h> // mmap
#include <sys/stat.h> // mmap
#include <fcntl.h> // O_RDONLY, lock

#include <iostream>
#include <string>
#include <cerrno>
#include <map>
#include <set>
#include <memory>

using namespace std;



class BasicCachedVector {

  public:
	enum OpenMode { ReadOnly, ReadWrite};
	int openMode[2] = { O_RDONLY, O_RDWR | O_CREAT | O_TRUNC };
	int mmapProt[2] = { PROT_READ, PROT_READ|PROT_WRITE };

	struct FileException : public std::exception { 
		FileException(int e): error_number(e){};
		const char* what() const noexcept override{ return "open on file failed"; }
		int error_number;
	};

	BasicCachedVector() {} // members locally initialized
	BasicCachedVector(string fn, BasicCachedVector::OpenMode m) : m_filename(fn), m_openmode(m) {
		//cerr << "[CONSTRUCTOR] " << fn<< endl;
		m_fd = open(filename().c_str(), BasicCachedVector::openMode[m], 0600);
		m_filesize = 0;
		if (m_fd >= 0){ 
			struct stat buf;
			fstat(m_fd,&buf);
			m_filesize = buf.st_size;
			//cerr << "[INFO] opened " << m_filename << " as cached vector of size in bytes: " << m_filesize <<endl;
		} else {
			cerr << "[ERROR] cannot open " << m_filename << endl;
			throw FileException(errno);
		}
	}
	
	BasicCachedVector(BasicCachedVector&& v) 
		: m_address(v.m_address)
		, m_fd(v.m_fd)
		, m_filename(v.m_filename)
		, m_filesize(v.m_filesize)
		, m_openmode(v.m_openmode) {

		v.m_filename="";
		v.m_filesize=0;
		v.m_address=0;
		v.m_openmode=ReadOnly;
		v.m_fd=-1;
	}

	~BasicCachedVector() {
		if (m_fd >= 0) { munmap_file();	close(m_fd); }
	}

	inline const string& filename() const { return m_filename; }

  protected:
	// getter functions
	inline BasicCachedVector::OpenMode openmode() const { return m_openmode; }
	inline int file_descriptor() const { return m_fd; }
	inline size_t filesize() const { return m_filesize; }

	void* map_resize(size_t fs) {
		cerr << "[map_resize] "<<this->filename() << " " << fs << " bytes"<< endl;
		munmap_file();
		resize_file(fs);
		return map_file();
	}

	void* map_file() {
		void* r;
		r = mmap(m_address, m_filesize, BasicCachedVector::mmapProt[m_openmode], MAP_SHARED, file_descriptor(), 0);
		m_address = 0;
		if (m_filesize == 0) return 0;
		if (r == (void*)-1) throw FileException(errno);
		m_address = r;
		cerr << "[map_file] filename="<< filename() << " of " << m_filesize << " bytes / fd=" << m_fd << endl;
		return r;
	}


  private:
	void munmap_file() { 
		cerr << "[munmap_file] fd=" << m_fd << " / m_filesize=" << m_filesize << " / filename= " << filename() << endl;
		munmap(m_address, m_filesize);
		m_address = 0;
	}

	inline void resize_file(size_t new_size) { // new_size in bytes
		cerr << "[resize_file] " << new_size << " / fd="<< m_fd << endl;
		if (ftruncate(m_fd, new_size)== -1) throw std::invalid_argument("ftruncate in CachedVector.h");
		m_filesize = new_size;
	}

  private:
	void* m_address = 0; // in subclass T*
	int m_fd = -1;
	string m_filename = "";
	size_t m_filesize = 0;
	BasicCachedVector::OpenMode m_openmode = ReadOnly;

};


// ### READ ###
template<class T>
class CachedVector : public BasicCachedVector {
public:
	typedef T value_type;
	typedef size_t size_type;
	typedef ptrdiff_t difference_type;
	typedef value_type& reference;
	typedef const value_type& consreference_t;
	typedef value_type* pointer;
	typedef const value_type* conspointer_t;
	typedef T* iterator;
	typedef const T* consiterator_t;

	CachedVector() {}
	CachedVector(string fn) : CachedVector(fn, ReadOnly){} // public constructor
	CachedVector(CachedVector&& v) // move constructor
		: BasicCachedVector(v)
		, p_start(v.p_start)
		, p_end(v.p_end)
		, p_reserve(v.p_reserve) {
		cerr << "move constructor"<< endl;
		v.p_start=v.p_end=v.p_reserve=0;
	}

protected:
	CachedVector(string fn, BasicCachedVector::OpenMode m) : BasicCachedVector(fn,m) {
		void* mm=map_file();
		p_start = (T*) mm;
		if (p_start) p_end = p_reserve = p_start+(filesize()/sizeof(T));
	}
public:

	// element read access
	inline const T& at(size_type i) const{ return *(p_start+i); }
	inline const T& operator[] (int i) const { return *(p_start+i); } 
	inline const T& front() const { return *p_start; }
	inline const T& back() const { return *(p_end-1); }

	// iterators
	inline iterator begin() { return p_start;}
	inline iterator end() { return p_end;}
	inline consiterator_t begin() const { return p_start; } 
	inline consiterator_t end() const { return p_end; }	 
	inline consiterator_t cbegin() const { return p_start; }
	inline consiterator_t cend() const { return p_end; }	

	// capacity
	//inline bool empty() const { return p_start; }
	inline size_type size() const { return (this->p_end - this->p_start); }
	inline size_type capacity() const { return this->p_reserve - this->p_start; }
	inline size_type max_size() const { return 1ul<<36; } // about 4G elements

	inline bool empty() const { return p_end!=p_start; }

  protected:
	T* p_start=0;
	T* p_end=0;
	T* p_reserve=0;
 };


// ### WRITE ###
template<class T>
class WriteCachedVector : public CachedVector<T> {
  public:
	typedef T value_type;
	typedef size_t size_type;
	typedef ptrdiff_t difference_type;
	typedef value_type& reference;
	typedef const value_type& consreference_t;
	typedef value_type* pointer;
	typedef const value_type* conspointer_t;
	typedef T* iterator;
	typedef const T* consiterator_t;

	WriteCachedVector(string fn) : CachedVector<T>(fn, CachedVector<T>::ReadWrite){
		if (this->p_start )cerr << "[WriteCachedVector(fn)]" << *this->p_start << endl;
	}

	~WriteCachedVector(){ 
		cerr << "[DESTRUCTOR] of WriteChachedVector" << endl;
		this->shrink_to_fit();
	}

	// element read/write access
	T& at(size_type i) { return (*(this->p_start+i)); }
	inline value_type& operator[] (int i) { return (*(this->p_start+i));  } 

	// iterators
	inline iterator begin() { return this->p_start; }
	inline iterator end() { return this->p_end; }	 

  public:
	// reserve space for inserting future elements.
	void reserve(size_t new_size) { // new_size is total map size
		cerr << "[reserve (" << new_size << " elements)]" << endl;
		//if (size_t(this->p_reserve-this->p_end) < new_size) { //shrinks also
			T* s = (T*) this->map_resize(new_size*sizeof(T));
			this->p_end=s+(this->p_end-this->p_start);
			this->p_start=s;
			this->p_reserve=this->p_start+new_size;
		//}
	}
	void shrink_to_size(size_t new_size){
		this->p_end = this->p_start+new_size;
	}

	// modifiers
	inline void push_back(T val){
		intelligenreserve_t(this->size() +1);
		*this->p_end = val;
		++this->p_end;
	}

	inline void pop_back(){ --this->p_end; }

  protected:
	void intelligenreserve_t(size_t s) { // for s elements in total map
		if (size_t(this->p_reserve-this->p_start) > s) { return; }
		//size_t new_size = ((this->p_reserve-this->p_start+getpagesize()-1)/getpagesize())*getpagesize();
		size_t new_size = (this->p_reserve-this->p_start)*2;
		if (new_size == 0) new_size = getpagesize();
		while (new_size < s) {
			new_size*=2;
		}
		reserve(new_size);		
	}

  private:
	void shrink_to_fit(){ // onlz called in destructor of WCV
		cerr << "[shrink_to_fit] calls reserve " << this->size() << endl;
		if (this->p_reserve < this->p_end) this->p_end = this->p_reserve; 
		reserve(this->size());
	}

};

class NamedVectorBase {
protected:
	static set<string>* m_usedNames;

public:
	bool fileNewer(string filename, string dep) {
		struct stat buf,buf2;
		if (stat(filename.c_str(),&buf)) return false;
		if (stat(dep.c_str(),&buf2)) return false;
		return buf2.st_mtime<buf.st_mtime;
	}

	static string cacheFilename(string project, string name) {
		return project.substr(0,project.find_last_of("\\/")) + "/" +name;
	}

	static string cacheFilename(string project, string name, int param) {
		return project.substr(0,project.find_last_of("\\/")) + "/" +name + "_" + to_string(param);
	}

	void makeCacheDir(string project) {
		mkdir((project+".d").c_str(),0700);
	}
};

template<class T>
class NamedVector : public NamedVectorBase {
public:	
	class NameDuplicate {};
	class UnknownName {};

	// stored generator classes, independent of project
	NamedVector(string identifier) {
		string name = identifier;
		if (!m_usedNames) m_usedNames=new set<string>;
		if (!m_usedNames->insert(name).second) throw NameDuplicate();
		if (!m_names) m_names=new map<string, NamedVector*>;
		(*m_names)[name]=this;
	}

	// name -> name_bla
	static shared_ptr<CachedVector<T>> getNamedVector(string project, string name, int param=-1) {
		//m_names has no parameter!
		string cached_vec_name = name;
		if (param != -1) cached_vec_name += "_" + to_string(param);
		weak_ptr<CachedVector<T>> cached = m_vectors[make_pair(project,cached_vec_name)];
		shared_ptr<CachedVector<T>> cached_vec = cached.lock();
		if (cached_vec) return cached_vec;

		auto named_vector = m_names->find(name);
		if (named_vector == m_names->end()){
			cerr << "name " << name << " of vector is unknown..." << endl;
			throw UnknownName();
		}
		cached_vec = named_vector->second->generateVector(project, param);
		m_vectors[make_pair(project,name)] = cached_vec;
		return cached_vec;
	}

protected:
	virtual shared_ptr<CachedVector<T>> generateVector(string project, int param=0)=0;	

private:
	static map<string,NamedVector*>* m_names;	//name -> generator class
	static map<pair<string,string>,weak_ptr<CachedVector<T>>> m_vectors;	// project + name+param -> actual cached vec
	void create_new_name(string name);

};

typedef uint32_t indices_t;

typedef shared_ptr<CachedVector<char>> data_t;
typedef shared_ptr<CachedVector<indices_t>> sorted_t;
typedef shared_ptr<CachedVector<uint16_t>> prelen_t;
typedef shared_ptr<CachedVector<uint16_t>> strlen_t;
typedef shared_ptr<CachedVector<indices_t>> protstart_t;
typedef shared_ptr<CachedVector<uint32_t>> protcount_t;
typedef shared_ptr<CachedVector<uint16_t>> common_prefix_t;
typedef shared_ptr<CachedVector<uint32_t>> protein_names_t;
typedef shared_ptr<CachedVector<uint32_t>> orgstart_t;
typedef shared_ptr<CachedVector<indices_t>> connected_component_mt2_t;
typedef shared_ptr<CachedVector<indices_t>> diverse_frags_t;
typedef shared_ptr<CachedVector<float>> CV_float_t;

#endif

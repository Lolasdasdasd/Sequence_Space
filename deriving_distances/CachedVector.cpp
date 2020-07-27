#include "CachedVector.h"
#include <memory>

set<string>* NamedVectorBase::m_usedNames;

template<class T>
map<string,NamedVector<T>*>* NamedVector<T>::m_names;

template<class T>
map<pair<string,string>,weak_ptr<CachedVector<T>>> NamedVector<T>::m_vectors;

template class NamedVector<char>;
template class NamedVector<uint8_t>;
template class NamedVector<uint16_t>;
template class NamedVector<uint32_t>;
template class NamedVector<float>;

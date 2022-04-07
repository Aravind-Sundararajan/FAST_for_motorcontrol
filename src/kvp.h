#ifndef HEADERFILE_KVP_
#define HEADERFILE_KVP_
template <typename T,typename U>
class kvp //Key Value Pair
{
public:
	T first;
	U second;
	kvp(){
		this->first = T();
		this->second =U();
	}
	kvp(const T& f, const U& s)
	{
		this->first = T(f);
		this->second = U(s);
	}
	kvp(const kvp& k)
	{
		first = k.first;
		second = k.second;
	}
	kvp<T,U>& operator=(const kvp& k)
	{
		(*this).first = k.first;
		(*this).second = k.second;
		return *this;
	}
	virtual ~kvp(){};
	kvp<U,T> swap()
	{
		kvp<U,T> newPair(second, first);
		return newPair;
	}
};
#endif

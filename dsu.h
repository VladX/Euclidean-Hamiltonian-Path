#include <vector>
#include <algorithm>

template<typename T = int>
class DSU {
private:
	T * rk;
	T * p;
public:
	DSU (size_t n) {
		rk = new T[n];
		p = new T[n];
		for (size_t i = 0; i < n; ++i)
			p[i] = i, rk[i] = 1;
	}
	
	~DSU () {
		delete[] rk;
		delete[] p;
	}
	
	inline T get (T v) {
		T cur = v;
		while (cur != p[cur])
			cur = p[cur];
		T ret = cur;
		cur = v;
		while (cur != p[cur])
			v = p[cur], p[cur] = ret, cur = v;
		return ret;
	}

	void unite (T a, T b) {
		a = get(a), b = get(b);
		if (a == b)
			return;
		if (rk[a] < rk[b])
			std::swap(a,b);
		p[b] = a;
		rk[a] += rk[b];
	}
};
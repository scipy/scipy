#ifndef _CYCLIC_CACHE_H_
#define _CYCLIC_CACHE_H_

namespace fft {

class CacheId {
	public:
		CacheId(int n) : m_n(n) {};
                virtual ~CacheId() {};

		virtual bool operator==(const CacheId& other) const
		{
			return is_equal(other);
		};

		virtual bool is_equal(const CacheId& other) const
		{
			return m_n == other.m_n;
		};

	public:
		int m_n;
		
};

template <class T>
class Cache {
	public:
		Cache() {};
		Cache(const T& id) : m_id(id) {};
		virtual ~Cache() {};

		virtual bool operator==(const Cache& other) const 
		{ 
			return other.m_id == m_id;
		};

	public:
		T m_id;
};

template <class T, class U>
class CacheManager {
	public:
		CacheManager(int n) :
			m_n(n),
			m_curn(0),
			m_last(0)
		{
			m_cache = new U*[n];

		};

		virtual ~CacheManager()
		{
			int i;
	
			for (i = 0; i < m_curn; ++i) {
				delete m_cache[i];
			}

			delete[] m_cache;
		}

		virtual U* get_cache(const T& id)
		{
			int i;

			/* Look in the current cache */
			for (i = 0; i < m_curn; ++i) {
				if ( m_cache[i]->m_id == id) {
					m_last = i;
					return m_cache[i];
				}
			}

			/* If still space, create a new cache */
			if (m_curn < m_n) {
				i = m_curn;
				++m_curn;
			} else {
				i = (m_last < m_n - 1) ? m_last + 1 : 0;
				delete m_cache[i];
			}

			m_cache[i] = new U(id);
			m_last = i;
			return m_cache[i];
		};

	private:
		U** m_cache;		
		int m_n;
		int m_curn;
		int m_last;
};

}

#endif

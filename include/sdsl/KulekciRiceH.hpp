#ifndef INCLUDED_RiceH
#define INCLUDED_RiceH

#include "sdsl/int_vector.hpp"
#include "sdsl/util.hpp"
#include <set> // for calculating the alphabet size
#include <map> // for mapping a symbol to its lexicographical index
#include <algorithm> // for std::swap
#include <stdexcept>
#include <vector>
#include <queue>
#include <utility>
#include "sdsl/select_support_mcl.hpp"

//! Namespace for the succinct data structure library.
namespace sdsl
{

template<class bit_vector_type = bit_vector, class select_support_type = select_support_mcl<1,1>>
class KulekciRiceH
{
    public:
        typedef int_vector<>::size_type              size_type;
        typedef int_vector<>::value_type             value_type;
        typedef random_access_const_iterator<KulekciRiceH> const_iterator;
        typedef const_iterator                       iterator;
        typedef wt_tag                               index_category;
        typedef int_alphabet_tag                     alphabet_category;
        typedef typename bit_vector::difference_type difference_type;
        enum 	{lex_ordered=1};

        typedef std::pair<value_type, size_type>     point_type;
        typedef std::vector<point_type>              point_vec_type;
        typedef std::pair<size_type, point_vec_type> r2d_res_type;


    protected:

        //mutable int_vector_type m_data;     // array keeps track of path offset in select-like methods
        bit_vector_type m_unaryBits;
        bit_vector m_binaryBits;
        select_support_type  m_select_1_unaryBits;
        uint64_t m_size;
        const int m_factor = 5;

        void copy(const KulekciRiceH& wt) {
        	m_unaryBits      = wt.m_unaryBits;
        	m_binaryBits      = wt.m_binaryBits;
        	m_select_1_unaryBits      = wt.m_select_1_unaryBits;
        	m_size      = wt.m_size;
        }


    public:
	const size_type&       sigma = 1;

    template<uint8_t int_width>
    KulekciRiceH(int_vector_buffer<int_width>& buf, size_type size)  {

        if (0 == size)
            return;
        size_type n = buf.size();  // set n
        if (n < size) {
            throw std::logic_error("n="+util::to_string(n)+" < "+util::to_string(size)+"=size");
            return;
        }

        bit_vector unaryBits = {1};
        bit_vector binaryBits;
        select_support_type select_1_unaryBits;
        const uint64_t jthBit[64] =  // jthBit[i] = 2^i
        {
            1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8192,16384,32768,65536,131072,262144,524288,1048576,2097152,
            4194304,8388608,16777216,33554432,67108864,134217728,268435456,536870912,1073741824,2147483648,4294967296,
            8589934592,17179869184,34359738368,68719476736,137438953472,274877906944,549755813888,1099511627776,2199023255552,
            4398046511104,8796093022208,17592186044416,35184372088832,70368744177664,140737488355328,281474976710656,
            562949953421312,1125899906842620,2251799813685250,4503599627370500,9007199254740990,18014398509482000,
            36028797018964000,72057594037927900,144115188075856000,288230376151712000,576460752303423000,1152921504606850000,
            2305843009213690000,4611686018427390000
        };

        size_t unaryLength = 1,binaryLength=0,i,j;
        unsigned long a,b,r,one=1;


        m_size = n;

		binaryBits.resize(n*m_factor);
		for(i=0; i<n ; i++)
		{
			a = (buf[i] >> m_factor);
			unaryBits.resize(unaryLength+a+1);
			for(j=0; j<a; j++) unaryBits[unaryLength++] = 0; //if you dont do this initialization, some of the padded bits may not be set to 0, which corrupts the coding then.
			unaryBits[unaryLength++]=1;

			for(j=0; j<m_factor; j++) binaryBits[i*m_factor+j] = buf[i] & jthBit[j];
		}
		unaryBits.resize(unaryBits.size()+64);//pad 64bits
		binaryBits.resize(binaryBits.size()+64);//pad 64bits

        m_unaryBits = bit_vector_type(std::move(unaryBits));
		util::init_support(m_select_1_unaryBits, &m_unaryBits);

		m_binaryBits.swap(binaryBits);

    }

    value_type operator[](size_type i)const {
        assert(i < size());
        unsigned long a,b,r,one=1;
		unsigned long mask=(one<<m_factor)-1;
		a = m_select_1_unaryBits(i+2) - m_select_1_unaryBits(i+1) - 1;
		r = ( a << m_factor ) + (m_binaryBits.get_int(i*m_factor) & mask);
		return r;
    };

        //! Default constructor
        KulekciRiceH() {m_size = 0;
        };



        //! Copy constructor
        KulekciRiceH(const KulekciRiceH& wt) {
            copy(wt);
        }

        //! Copy constructor
        KulekciRiceH(KulekciRiceH&& wt) {
            *this = std::move(wt);
        }

        //! Assignment operator
        KulekciRiceH& operator=(const KulekciRiceH& wt) {
            if (this != &wt) {
                copy(wt);
            }
            return *this;
        }

        //! Assignment move operator
        KulekciRiceH& operator=(KulekciRiceH&& wt) {
            if (this != &wt) {
            	m_binaryBits      = std::move(wt.m_binaryBits);
            	m_select_1_unaryBits      = std::move(wt.m_select_1_unaryBits);
            	m_unaryBits      = std::move(wt.m_unaryBits);
            	m_size      = std::move(wt.m_size);
            }
            return *this;
        }

        //! Swap operator
        void swap(KulekciRiceH& wt) {
            if (this != &wt) {
                m_binaryBits.swap(wt.m_binaryBits);

            	m_unaryBits.swap(wt.m_unaryBits);
            	util::swap_support(m_select_1_unaryBits, wt.m_select_1_unaryBits, &m_unaryBits, &(wt.m_unaryBits));
            	m_size = wt.m_size;
            }
        }

        //! Returns the size of the original vector.
        size_type size()const {
            return m_size;
        }

        //! Returns whether the wavelet tree contains no data.
        bool empty()const {
            return m_binaryBits.empty();
        }

        //! Calculates how many symbols c are in the prefix [0..i-1] of the supported vector.
        /*!
         *  \param i The exclusive index of the prefix range [0..i-1], so \f$i\in[0..size()]\f$.
         *  \param c The symbol to count the occurrences in the prefix.
         *    \returns The number of occurrences of symbol c in the prefix [0..i-1] of the supported vector.
         *  \par Time complexity
         *       \f$ \Order{\log |\Sigma|} \f$
         *  \par Precondition
         *       \f$ i \leq size() \f$
         */
        size_type rank(size_type i, value_type c)const {
            return 0;
        };



        //! Calculates how many occurrences of symbol wt[i] are in the prefix [0..i-1] of the original sequence.
        /*!
         *  \param i The index of the symbol.
         *  \return  Pair (rank(wt[i],i),wt[i])
         *  \par Precondition
         *       \f$ i < size() \f$
         */
        std::pair<size_type, value_type>
        inverse_select(size_type i)const {
            return         std::pair<size_type, value_type>(0,0);
        }

        //! Calculates the i-th occurrence of the symbol c in the supported vector.
        /*!
         *  \param i The i-th occurrence.
         *  \param c The symbol c.
         *  \par Time complexity
         *       \f$ \Order{\log |\Sigma|} \f$
         *  \par Precondition
         *       \f$ 1 \leq i \leq rank(size(), c) \f$
         */
        size_type select(size_type i, value_type c)const {
            return 0;
        };
        void interval_symbols(size_type i, size_type j, size_type& k,
                              std::vector<value_type>& cs,
                              std::vector<size_type>& rank_c_i,
                              std::vector<size_type>& rank_c_j) const {
        }

        template<class t_ret_type = std::tuple<size_type, size_type, size_type>>
        t_ret_type lex_count(size_type i, size_type j, value_type c)const {
            assert(i <= j and j <= size());
                        return t_ret_type {i, 0, 0};
        };
        template<class t_ret_type = std::tuple<size_type, size_type>>
        t_ret_type lex_smaller_count(size_type i, value_type c) const {

            return t_ret_type {i, 0};
        }
        std::pair<size_type, std::vector<std::pair<value_type, size_type>>>
        range_search_2d(size_type lb, size_type rb, value_type vlb, value_type vrb,
                        bool report=true) const {
                        return make_pair(0, std::vector<std::pair<value_type, size_type>>());
        }
        const_iterator begin()const {
            return const_iterator(this, 0);
        }

        //! Returns a const_iterator to the element after the last element.
        const_iterator end()const {
            return const_iterator(this, size());
        }


        //! Serializes the data structure into the given ostream
        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += m_binaryBits.serialize(out, child, "m_binaryBits");
            written_bytes += m_unaryBits.serialize(out, child, "m_unaryBits");
            written_bytes += m_select_1_unaryBits.serialize(out, child, "m_select_1_unaryBits");
            written_bytes += write_member(this->m_size,out,child, "size");

            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Loads the data structure from the given istream.
        void load(std::istream& in) {
            m_binaryBits.load(in);
            m_unaryBits.load(in);
            m_select_1_unaryBits.load(in, &m_unaryBits);
            read_member(this->m_size, in);

        }

        //! return the path to the leaf for a given symbol
        std::pair<uint64_t,uint64_t> path(value_type c) const {
            return {0,c};
        }
};

}// end namespace sdsl
#endif


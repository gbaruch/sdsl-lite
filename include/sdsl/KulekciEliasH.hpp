#ifndef INCLUDED_EliasH
#define INCLUDED_EliasH

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
class KulekciEliasH
{
    public:
        typedef int_vector<>::size_type              size_type;
        typedef int_vector<>::value_type             value_type;
        typedef random_access_const_iterator<KulekciEliasH> const_iterator;
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

        void copy(const KulekciEliasH& wt) {
        	m_unaryBits      = wt.m_unaryBits;
        	m_binaryBits      = wt.m_binaryBits;
        	m_select_1_unaryBits      = wt.m_select_1_unaryBits;
        	m_size      = wt.m_size;
        }


    public:
	const size_type&       sigma = 1;

    template<uint8_t int_width>
    KulekciEliasH(int_vector_buffer<int_width>& buf, size_type size)  {

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

        size_t unaryLength = 1,binaryLength=0,i,j;
        unsigned long a,b,r,one=1;


        m_size = n;
		for(i=0; i<n ; i++)
		{
			a = buf[i]+2;
			r = bits::hi(a);
			unaryLength  += r+1;
			unaryBits.resize(unaryLength);
			for(j=unaryLength-r-1; j<unaryLength; j++) unaryBits[j]=0; //if you dont do this initialization, some of the padded bits may not be set to 0, which corrupts the coding then.
			unaryBits[unaryLength-1]=1;
			binaryBits.resize(binaryLength+r);
			for(j=0; j<r; j++) binaryBits[binaryLength++] = a & (one<<j);
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
        a = m_select_1_unaryBits(i+1);//a-i is the beginning bit address of x[i] on the binaryBits array
		b = m_unaryBits.get_int(a+1);//instead of performing a second select, seek for the next set bit via bitwise operations, thanks to Simon Gog for pointing out this
		r = m_binaryBits.get_int(a-i);
		a = (one << bits::lo(b))-1;
		r = r & a;
		r = r + a + 1;
        return r-2;
    };

        //! Default constructor
        KulekciEliasH() {m_size = 0;
        };



        //! Copy constructor
        KulekciEliasH(const KulekciEliasH& wt) {
            copy(wt);
        }

        //! Copy constructor
        KulekciEliasH(KulekciEliasH&& wt) {
            *this = std::move(wt);
        }

        //! Assignment operator
        KulekciEliasH& operator=(const KulekciEliasH& wt) {
            if (this != &wt) {
                copy(wt);
            }
            return *this;
        }

        //! Assignment move operator
        KulekciEliasH& operator=(KulekciEliasH&& wt) {
            if (this != &wt) {
            	m_binaryBits      = std::move(wt.m_binaryBits);
            	m_select_1_unaryBits      = std::move(wt.m_select_1_unaryBits);
            	m_unaryBits      = std::move(wt.m_unaryBits);
            	m_size      = std::move(wt.m_size);
            }
            return *this;
        }

        //! Swap operator
        void swap(KulekciEliasH& wt) {
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

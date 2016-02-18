#ifndef INCLUDED_KulekciRiceW
#define INCLUDED_KulekciRiceW

#include "sdsl/int_vector.hpp"
#include "sdsl/util.hpp"
#include <set> // for calculating the alphabet size
#include <map> // for mapping a symbol to its lexicographical index
#include <algorithm> // for std::swap
#include <stdexcept>
#include <vector>
#include <queue>
#include <utility>

//! Namespace for the succinct data structure library.
namespace sdsl
{

template<class wt_type>
class KulekciRiceW
{
    public:
        typedef int_vector<>::size_type              size_type;
        typedef int_vector<>::value_type             value_type;
        typedef random_access_const_iterator<KulekciRiceW> const_iterator;
        typedef const_iterator                       iterator;
        typedef wt_tag                               index_category;
        typedef int_alphabet_tag                     alphabet_category;
        typedef typename bit_vector::difference_type difference_type;
        enum 	{lex_ordered=1};

        typedef std::pair<value_type, size_type>     point_type;
        typedef std::vector<point_type>              point_vec_type;
        typedef std::pair<size_type, point_vec_type> r2d_res_type;


    protected:

        wt_type m_wt;
        bit_vector m_rawBits;
        uint64_t m_size;
        const int m_factor = 5;


    public:
	const size_type&       sigma = 1;

    template<uint8_t int_width>
    KulekciRiceW(int_vector_buffer<int_width>& buf, size_type size)  {

        if (0 == size)
            return;
        size_type n = buf.size();  // set n
        if (n < size) {
            throw std::logic_error("n="+util::to_string(n)+" < "+util::to_string(size)+"=size");
            return;
        }

        int_vector<> labels;
        size_t i,j,bitvectorSize;

        labels.resize(n);

	    uint64_t one=1;
	    m_rawBits.resize(n * m_factor);

        m_size = n;
		for(i=0; i<n; i++)
	    {
	        labels[i]= buf[i] >> m_factor; //compute the label
	        //compute the payload, which is the least significant factor bits. Append these bits to the rawBits
	        for(j=0; j<m_factor; j++) m_rawBits[i*m_factor+j] = buf[i] & (one<<j);//insert from least significant to most significant to be able to retrieve in one 64-bit read operation
	    }

m_rawBits.resize(m_rawBits.size() + 64);
	    construct_im(m_wt, labels);//create the huffman shaped wavelet tree of the labels
    }

    value_type operator[](size_type i)const {
        assert(i < size());
        size_t j, label, value;
        
        uint64_t mask,maskPlus1,one=1;
	    mask = (one<<m_factor)-1;
	    maskPlus1 = mask+1;
	    
		label     = m_wt[i];
        value     = m_rawBits.get_int(i*m_factor) & mask;
        value    += label*maskPlus1;
        return value;
	};

void copy(const KulekciRiceW& wt) {
        	m_wt      = wt.m_wt;
        	m_rawBits      = wt.m_rawBits;
        	m_size      = wt.m_size;
        }

        //! Default constructor
        KulekciRiceW() {m_size = 0;
        };



        //! Copy constructor
        KulekciRiceW(const KulekciRiceW& wt) {
            copy(wt);
        }

        //! Copy constructor
        KulekciRiceW(KulekciRiceW&& wt) {
            *this = std::move(wt);
        }

        //! Assignment operator
        KulekciRiceW& operator=(const KulekciRiceW& wt) {
            if (this != &wt) {
                copy(wt);
            }
            return *this;
        }

        //! Assignment move operator
        KulekciRiceW& operator=(KulekciRiceW&& wt) {
            if (this != &wt) {
            	m_wt      = std::move(wt.m_wt);
            	m_rawBits      = std::move(wt.m_rawBits);
            	m_size      = std::move(wt.m_size);
            }
            return *this;
        }

        //! Swap operator
        void swap(KulekciRiceW& wt) {
            if (this != &wt) {
                m_wt.swap(wt.m_wt);

                	m_rawBits.swap(wt.m_rawBits);
            	m_size = wt.m_size;
            }
        }

        //! Returns the size of the original vector.
        size_type size()const {
            return m_size;
        }

        //! Returns whether the wavelet tree contains no data.
        bool empty()const {
            return m_rawBits.empty();
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
            written_bytes += m_wt.serialize(out, child, "m_wt");
            
            written_bytes += m_rawBits.serialize(out, child, "m_rawBits");
            
            written_bytes += write_member(this->m_size,out,child, "size");

            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Loads the data structure from the given istream.
        void load(std::istream& in) {
        	m_wt.load(in);
        	m_rawBits.load(in);
            read_member(this->m_size, in);

        }

        //! return the path to the leaf for a given symbol
        std::pair<uint64_t,uint64_t> path(value_type c) const {
            return {0,c};
        }
};

}// end namespace sdsl
#endif

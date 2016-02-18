/* sdsl - succinct data structures library
    Copyright (C) 2009 Simon Gog

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/ .
*/
/*! \file wt_int.hpp
    \brief wt_int.hpp contains a specialized class for a wavelet tree of a
           sequence of the numbers. This wavelet tree class takes
           less memory than the wt_pc class for large alphabets.
    \author Simon Gog, Shanika Kuruppu
*/
#ifndef INCLUDED_FIXED_LENGTH
#define INCLUDED_FIXED_LENGTH

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

//! A wavelet tree class for integer sequences.
/*!
 *    \par Space complexity
 *        \f$\Order{n\log|\Sigma|}\f$ bits, where \f$n\f$ is the size of the vector the wavelet tree was build for.
 *
 *  \tparam t_bitvector   Type of the bitvector used for representing the wavelet tree.
 *  \tparam t_rank        Type of the support structure for rank on pattern `1`.
 *  \tparam t_select      Type of the support structure for select on pattern `1`.
 *  \tparam t_select_zero Type of the support structure for select on pattern `0`.
 *
 *   @ingroup wt
 */
template<class t_intvector   = int_vector<0>>
class wt_fixedLength
{
    public:
		typedef t_intvector                          int_vector_type;
        typedef int_vector<>::size_type              size_type;
        typedef int_vector<>::value_type             value_type;
        typedef random_access_const_iterator<wt_fixedLength> const_iterator;
        typedef const_iterator                       iterator;
        typedef wt_tag                               index_category;
        typedef int_alphabet_tag                     alphabet_category;
        enum 	{lex_ordered=1};

        typedef std::pair<value_type, size_type>     point_type;
        typedef std::vector<point_type>              point_vec_type;
        typedef std::pair<size_type, point_vec_type> r2d_res_type;


    protected:

        mutable int_vector_type m_data;     // array keeps track of path offset in select-like methods

        void copy(const wt_fixedLength& wt) {
            m_data      = wt.m_data;
        }


    public:
	const size_type&       sigma = 1;


        //! Default constructor
        wt_fixedLength() {
        };

        //! Semi-external constructor
        /*! \param buf         File buffer of the int_vector for which the wt_int should be build.
         *  \param size        Size of the prefix of v, which should be indexed.
         *  \param max_level   Maximal level of the wavelet tree. If set to 0, determined automatically.
         *    \par Time complexity
         *        \f$ \Order{n\log|\Sigma|}\f$, where \f$n=size\f$
         *        I.e. we need \Order{n\log n} if rac is a permutation of 0..n-1.
         *    \par Space complexity
         *        \f$ n\log|\Sigma| + O(1)\f$ bits, where \f$n=size\f$.
         */
        template<uint8_t int_width>
        wt_fixedLength(int_vector_buffer<int_width>& buf, size_type size)  {

            if (0 == size)
                return;
            size_type n = buf.size();  // set n
            if (n < size) {
                throw std::logic_error("n="+util::to_string(n)+" < "+util::to_string(size)+"=size");
                return;
            }

            int_vector<0> rac(size, 0, buf.width());

			for (size_type i=0; i < size; ++i) {
				rac[i] = buf[i];
			}

			util::bit_compress(rac);
			m_data = int_vector_type(std::move(rac));

        }

        //! Copy constructor
        wt_fixedLength(const wt_fixedLength& wt) {
            copy(wt);
        }

        //! Copy constructor
        wt_fixedLength(wt_fixedLength&& wt) {
            *this = std::move(wt);
        }

        //! Assignment operator
        wt_fixedLength& operator=(const wt_fixedLength& wt) {
            if (this != &wt) {
                copy(wt);
            }
            return *this;
        }

        //! Assignment move operator
        wt_fixedLength& operator=(wt_fixedLength&& wt) {
            if (this != &wt) {
                m_data      = std::move(wt.m_data);
            }
            return *this;
        }

        //! Swap operator
        void swap(wt_fixedLength& wt) {
            if (this != &wt) {
            	m_data.swap(wt.m_data);
            }
        }

        //! Returns the size of the original vector.
        size_type size()const {
            return m_data.size();
        }

        //! Returns whether the wavelet tree contains no data.
        bool empty()const {
            return m_data.empty();
        }

        //! Recovers the i-th symbol of the original vector.
        /*! \param i The index of the symbol in the original vector.
         *  \returns The i-th symbol of the original vector.
         *  \par Precondition
         *       \f$ i < size() \f$
         */
        value_type operator[](size_type i)const {
            assert(i < size());
            return m_data[i];
        };

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


        //! For each symbol c in wt[i..j-1] get rank(i,c) and rank(j,c).
        /*!
         * \param i        The start index (inclusive) of the interval.
         * \param j        The end index (exclusive) of the interval.
         * \param k        Reference for number of different symbols in [i..j-1].
         * \param cs       Reference to a vector that will contain in
         *                 cs[0..k-1] all symbols that occur in [i..j-1] in
         *                 ascending order.
         * \param rank_c_i Reference to a vector which equals
         *                 rank_c_i[p] = rank(i,cs[p]), for \f$ 0 \leq p < k \f$.
         * \param rank_c_j Reference to a vector which equals
         *                 rank_c_j[p] = rank(j,cs[p]), for \f$ 0 \leq p < k \f$.
         * \par Time complexity
         *      \f$ \Order{\min{\sigma, k \log \sigma}} \f$
         *
         * \par Precondition
         *      \f$ i \leq j \leq size() \f$
         *      \f$ cs.size() \geq \sigma \f$
         *      \f$ rank_{c_i}.size() \geq \sigma \f$
         *      \f$ rank_{c_j}.size() \geq \sigma \f$
         */
        void interval_symbols(size_type i, size_type j, size_type& k,
                              std::vector<value_type>& cs,
                              std::vector<size_type>& rank_c_i,
                              std::vector<size_type>& rank_c_j) const {
        }

        //! How many symbols are lexicographic smaller/greater than c in [i..j-1].
        /*!
         * \param i       Start index (inclusive) of the interval.
         * \param j       End index (exclusive) of the interval.
         * \param c       Symbol c.
         * \return A triple containing:
         *         * rank(i,c)
         *         * #symbols smaller than c in [i..j-1]
         *         * #symbols greater than c in [i..j-1]
         *
         * \par Precondition
         *      \f$ i \leq j \leq size() \f$
         */
        template<class t_ret_type = std::tuple<size_type, size_type, size_type>>
        t_ret_type lex_count(size_type i, size_type j, value_type c)const {
            assert(i <= j and j <= size());
                        return t_ret_type {i, 0, 0};
        };

        //! How many symbols are lexicographic smaller than c in [0..i-1].
        /*!
         * \param i Exclusive right bound of the range.
         * \param c Symbol c.
         * \return A tuple containing:
         *         * rank(i,c)
         *         * #symbols smaller than c in [0..i-1]
         * \par Precondition
         *      \f$ i \leq size() \f$
         */
        template<class t_ret_type = std::tuple<size_type, size_type>>
        t_ret_type lex_smaller_count(size_type i, value_type c) const {

            return t_ret_type {i, 0};
        }

        //! range_search_2d searches points in the index interval [lb..rb] and value interval [vlb..vrb].
        /*! \param lb     Left bound of index interval (inclusive)
         *  \param rb     Right bound of index interval (inclusive)
         *  \param vlb    Left bound of value interval (inclusive)
         *  \param vrb    Right bound of value interval (inclusive)
         *  \param report Should the matching points be returned?
         *  \return Pair (#of found points, vector of points), the vector is empty when
         *          report = false.
         */
        std::pair<size_type, std::vector<std::pair<value_type, size_type>>>
        range_search_2d(size_type lb, size_type rb, value_type vlb, value_type vrb,
                        bool report=true) const {
                        return make_pair(0, std::vector<std::pair<value_type, size_type>>());
        }


        //! Returns a const_iterator to the first element.
        const_iterator begin()const {
            return const_iterator(m_data, 0);
        }

        //! Returns a const_iterator to the element after the last element.
        const_iterator end()const {
            return const_iterator(m_data, size());
        }


        //! Serializes the data structure into the given ostream
        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += m_data.serialize(out, child, "vector_data");

            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Loads the data structure from the given istream.
        void load(std::istream& in) {
            m_data.load(in);
        }

        //! return the path to the leaf for a given symbol
        std::pair<uint64_t,uint64_t> path(value_type c) const {
            return {0,c};
        }
};

}// end namespace sdsl
#endif


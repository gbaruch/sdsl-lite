/*! \file wt_cann.hpp
    \brief wt_cann.hpp contains a class for a Canonical Huffman shaped wavelet tree
                       over byte sequences.
    \author Simon Gog and Timo Beller
*/
#ifndef INCLUDED_SDSL_WT_CANN
#define INCLUDED_SDSL_WT_CANN

#include <iostream>
#include <vector>
#include <sdsl/wt_pc.hpp>
#include <tuple>
#include <bitset>
//! Namespace for the succinct data structure library.
namespace sdsl
{

// forward declaration
struct cann_shape;

//! A Huffman-shaped wavelet tree.
/*!
 * A wavelet tree is build for a vector of characters over the byte alphabet
 * \f$\Sigma\f$. If you need a wavelet tree for a integer alphabet you should
 * use `wt_int`.
 * The wavelet tree \f$wt\f$ consists of a tree of bitvectors and provides
 * three efficient methods:
 *   - The "[]"-operator: \f$wt[i]\f$ returns the i-th symbol of vector for
 *     which the wavelet tree was build for.
 *   - The rank method: \f$wt.rank(i,c)\f$ returns the number of occurrences
 *     of symbol \f$c\f$ in the prefix [0..i-1] in the vector for which the
 *     wavelet tree was build for.
 *   - The select method: \f$wt.select(j,c)\f$ returns the index
 *     \f$i\in [0..size()-1]\f$ of the j-th occurrence of symbol \f$c\f$.
 *
 * The idea of using a Huffman shaped wavelet was first mentioned on page 17
 * of the following technical report:
 * Veli Mäkinen and Gonzalo Navarro:
 * ,,Succinct Suffix Arrays based on Run-Length Encoding.''
 * Available under: http://swp.dcc.uchile.cl/TR/2005/TR_DCC-2005-004.pdf
 *
 * \tparam t_bitvector   Underlying bitvector structure.
 * \tparam t_rank        Rank support for pattern `1` on the bitvector.
 * \tparam t_select      Select support for pattern `1` on the bitvector.
 * \tparam t_select_zero Select support for pattern `0` on the bitvector.
 * \tparam t_dfs_shape   Layout of the tree structure in memory. Set 0
 *                       for BFS layout and 1 fro DFS layout.
 *
 * \par Space complexity
 *      \f$n H_0 + 2|\Sigma|\log n\f$ bits, where \f$n\f$ is the size
 *       of the vector the wavelet tree was build for.
 *
 *  @ingroup wt
 */
template<class t_bitvector   = bit_vector,
         class t_rank        = typename t_bitvector::rank_1_type,
         class t_select      = typename t_bitvector::select_1_type,
         class t_select_zero = typename t_bitvector::select_0_type,
         class t_tree_strat  = byte_tree<> >
using wt_cann = wt_pc<cann_shape,
      t_bitvector,
      t_rank,
      t_select,
      t_select_zero,
      t_tree_strat>;

// Canonical Huffman shape for wt_pc
template<class t_wt>
struct _cann_shape {
    typedef typename t_wt::size_type         size_type;
    typedef std::pair<size_type, size_type>  tPII;    // (freq, nodenr)-pair
    typedef std::priority_queue
    <tPII, std::vector<tPII>,
    std::greater<tPII>>                      tMPQPII; // min priority queue
    enum { lex_ordered = 0 };


    template<class t_rac>
    static void
    construct_tree(t_rac& C, std::vector<pc_node>& temp_nodes) {
        tMPQPII pq;
        size_type i = 0;
        // add leaves of Huffman tree
        std::for_each(std::begin(C), std::end(C), [&](decltype(*std::begin(C)) &freq) {
            if (freq > 0) {
                        pq.push(tPII(freq, temp_nodes.size()));// push (frequency, node pointer)
                        // initial bv_pos with number of occurrences and bv_pos_rank
                        // value with the code of the corresponding char, parent,
                        // child[0], and child[1] are set to undef
                        temp_nodes.emplace_back(pc_node(freq, i));
			}
			++i;
		});

		while (pq.size() > 1) {
			tPII v1, v2;
			v1 = pq.top(); pq.pop();
			v2 = pq.top(); pq.pop();
			temp_nodes[v1.second].parent = temp_nodes.size(); // parent is new node
			temp_nodes[v2.second].parent = temp_nodes.size(); // parent is new node
			size_type frq_sum = v1.first + v2.first;
			pq.push(tPII(frq_sum, temp_nodes.size()));
			temp_nodes.emplace_back(pc_node(frq_sum, 0, pc_node::undef,
											v1.second, v2.second));
		}


        std::vector<std::tuple<uint64_t, uint64_t, uint64_t>> depthSymbolAndFrequency;
    	for( size_type j = 0; j < temp_nodes.size(); j++)
    	{
    		if(temp_nodes[j].child[0] == pc_node::undef)
    		{
    			uint64_t depth = 0;
    			uint64_t current = j;
    			while(temp_nodes[current].parent != pc_node::undef)
    			{
    				depth++;
    				current = temp_nodes[current].parent;
    			}
    			depthSymbolAndFrequency.push_back(std::tuple<uint64_t, uint64_t,  uint64_t>(depth,temp_nodes[j].sym, temp_nodes[j].freq));
    		}
    	}
	/*	std::cout << depthSymbolAndFrequency.size() << " objects in depthSymbolAndFrequency." <<std::endl;;
    	for(auto it = depthSymbolAndFrequency.begin(); it != depthSymbolAndFrequency.end(); it++)
    	{
    		std::cout << "depth " <<std::get<0>(*it) << "freq " << std::get<2>(*it) <<std::endl;
    	}*/
    	temp_nodes.clear();
    	size_type currentNode = 0; //todo: handle single char tree
        sdsl::bit_vector bv;
        std::map<sdsl::bit_vector, uint64_t> chars_coding;
        std::sort(depthSymbolAndFrequency.begin(), depthSymbolAndFrequency.end());
        temp_nodes.emplace_back(pc_node(0, 0));
        size_type nodesToGoDown = 0;
    	for(std::vector<std::tuple<uint64_t, uint64_t, uint64_t>>::iterator it = depthSymbolAndFrequency.begin(); it != depthSymbolAndFrequency.end(); it++)
    	{
    		//todo: update freqs
    		size_type bv_length_diff =std::get<0>(*it) - bv.size();
        	if(bv_length_diff > 0)
        	{
        		nodesToGoDown += bv_length_diff;
        		bv.bit_resize(bv.size() + bv_length_diff);
        	}


        	if(nodesToGoDown > 0)
        	{
        		for(; nodesToGoDown > 0; nodesToGoDown--)
        		{
        			temp_nodes[currentNode].minimal_depth = nodesToGoDown;
					temp_nodes[currentNode].child[0] = currentNode + 1;
					temp_nodes.emplace_back(pc_node(0, 0, currentNode));
        			currentNode++;

        		}
        	}

        	//update frequencies of ancestors
        	size_type tempNode = currentNode;
        	while(tempNode != pc_node::undef)
        	{
        		temp_nodes[tempNode].freq += std::get<2>(*it);
        		tempNode = temp_nodes[tempNode].parent;
        	}
        	temp_nodes[currentNode].sym = std::get<1>(*it);
        	chars_coding[bv] = std::get<1>(*it);

        	//adding 1. going up if needed
        	uint64_t temp = bv.size() - 1;
        	bool finished = false;
			while(bv[temp] == 1)
			{
				nodesToGoDown++;
				currentNode = temp_nodes[currentNode].parent;
				temp_nodes[currentNode].maximal_depth = nodesToGoDown;
				bv[temp] = 0;
				if(temp == 0)
				{
					finished = true;
					break;
				}
				temp--;
			}
			if(!finished)
			{
				currentNode = temp_nodes[currentNode].parent;
				temp_nodes[currentNode].child[1] = temp_nodes.size();
				temp_nodes.emplace_back(pc_node(0, 0, currentNode));
				currentNode = temp_nodes.size() - 1;
				bv[temp] = 1;
			}
        }


    	pc_node root = temp_nodes[0];
    	//root.minimal_depth = temp_nodes[0].minimal_depth;
    	//root.maximal_depth = temp_nodes[0].maximal_depth;
    	for (uint64_t i=1; i < temp_nodes.size(); ++i) {
			temp_nodes[i-1] = temp_nodes[i];
			temp_nodes[i-1].parent = (temp_nodes[i-1].parent+temp_nodes.size()-1)%temp_nodes.size();
			temp_nodes[i-1].child[0] -= (temp_nodes[i-1].child[0] != pc_node::undef);
			temp_nodes[i-1].child[1] -= (temp_nodes[i-1].child[1] != pc_node::undef);
			//temp_nodes[i-1].minimal_depth = temp_nodes[i].minimal_depth;
			//temp_nodes[i-1].maximal_depth = temp_nodes[i].maximal_depth;
		}
		root.child[0] -= (root.child[0] != pc_node::undef);
		root.child[1] -= (root.child[1] != pc_node::undef);
		temp_nodes[temp_nodes.size()-1] = root;
		//temp_nodes.back().minimal_depth = root.minimal_depth;
		//temp_nodes.back().maximal_depth = root.maximal_depth;

		for(std::map<sdsl::bit_vector, uint64_t>::iterator it = chars_coding.begin(); it != chars_coding.end();
				it++)
		{
			std::cout << "code of length "<< it->first.size()<< " for " << std::bitset<12>(it->first.get_int(0, it->first.size())).to_string() << " is " <<it->second << std::endl;
		}
    	/*//printing the tree
    	for(size_type i = 0; i < temp_nodes.size(); i++)
    	{
    		char symbol;
    		if(temp_nodes[i].sym)
    			symbol = (char)temp_nodes[i].sym;
    		else
    			symbol = '*';
    		std::cout << i << ": symbol: " << symbol << ", parent is " << temp_nodes[i].parent
    				<< ", left child is " << temp_nodes[i].child[0]  << ", right child is " << temp_nodes[i].child[1]
    				<< ", freq is " << temp_nodes[i].freq
    				<< ", minimal depth is " << temp_nodes[i].minimal_depth
    				<< ", maximal depth is " << temp_nodes[i].maximal_depth<< std::endl;
    	}*/
    }
};

struct cann_shape {
    template<class t_wt>
    using type = _cann_shape<t_wt>;
};


}// end namespace sdsl
#endif

/*
 * SkeletonWaveletTree.h
 *
 *  Created on: Apr 12, 2015
 *      Author: gilad
 */

#ifndef NEW_SKELETONWAVELETTREE_H_
#define NEW_SKELETONWAVELETTREE_H_

#include <iostream>
#include <sdsl/wt_helper.hpp>
#include <sdsl/wt_pc.hpp>
//#include <sdsl/wt_cann.hpp>
#include <cmath>

namespace sdsl
{
using namespace std::chrono;
using timer = std::chrono::high_resolution_clock;
struct pc_node_ex {
    uint64_t  freq;     // frequency of symbol sym
    uint64_t  sym;      // symbol
    uint64_t  parent;   // pointer to the parent
    uint64_t  child[2]; // pointer to the children
    uint64_t  minimal_depth;
    uint64_t  maximal_depth;

    enum :uint64_t {undef = 0xFFFFFFFFFFFFFFFFULL}; // max uint64_t value

    pc_node_ex(uint64_t freq_=0, uint64_t sym_=0, uint64_t parent_=undef,
            uint64_t child_left=undef, uint64_t child_right=undef)
    {
    	freq = freq_;
		sym = sym_;
		parent = parent_;
		child[0] = child_left;
		child[1] = child_right;
		minimal_depth = undef;
		maximal_depth = undef;
    }

    pc_node_ex& operator=(const pc_node_ex& v)
    {
    	if (this != &v) {
    	freq = v.freq;
    	sym = v.sym;
    	parent = v.parent;
    	child[0] = v.child[0];
    	child[1] = v.child[1];
    	minimal_depth = v.minimal_depth;
    	maximal_depth = v.maximal_depth;
    	}
    	return *this;
    }

    operator pc_node() const {
    	pc_node v;
    	v.freq = freq;
    	v.sym = sym;
    	v.parent = parent;
    	v.child[0] = child[0];
    	v.child[1] = child[1];
    	return v;
    }
};

typedef std::tuple<int_vector<>::size_type, uint64_t,int_vector<>::size_type>  tPII;    // (freq, sym ,nodenr)-tuple
typedef std::priority_queue<tPII, std::vector<tPII>,std::greater<tPII>>                      tMPQPII;
template<class t_shape,
         class t_bitvector   = bit_vector,
         class t_rank        = typename t_bitvector::rank_1_type,
         class t_select      = typename t_bitvector::select_1_type,
         class t_select_zero = typename t_bitvector::select_0_type,
         class t_tree_strat  = int_tree<>
         >
	class wt_new_sk: public sdsl::wt_pc<t_shape, t_bitvector, t_rank, t_select, t_select_zero, t_tree_strat>
	{
	public:
		typedef typename
		t_tree_strat::template type<wt_pc<t_shape, t_bitvector, t_rank, t_select, t_select_zero, t_tree_strat>>            tree_strat_type;
        typedef typename
        t_tree_strat::template type<wt_pc<t_shape, t_bitvector, t_rank, t_select, t_select_zero, t_tree_strat>>            tree_strat_type_pc;
		typedef int_vector<>::size_type               size_type;
		typedef typename
		tree_strat_type::value_type                   value_type;
		typedef typename t_bitvector::difference_type difference_type;
		typedef random_access_const_iterator<wt_pc<t_shape, t_bitvector, t_rank, t_select, t_select_zero, t_tree_strat>>   const_iterator;
		typedef const_iterator                        iterator;
		typedef t_bitvector                           bit_vector_type;
		typedef t_rank                                rank_1_type;
		typedef t_select                              select_1_type;
		typedef t_select_zero                         select_0_type;
		typedef wt_tag                                index_category;
		typedef typename
		tree_strat_type::alphabet_category            alphabet_category;
		typedef typename
		t_shape::template type<wt_pc<t_shape, t_bitvector, t_rank, t_select, t_select_zero, t_tree_strat>>                 shape_type;
		enum { lex_ordered=shape_type::lex_ordered };
		using node_type = typename tree_strat_type::node_type;

		//std::map<uint64_t ,sdsl::bit_vector> 	m_chars_coding;
		std::map<uint64_t, uint64_t>			m_path_to_char;
		bit_vector		 						m_leaves_bv;

		virtual ~wt_new_sk(){}

		//creates a canonical huffman tree, represented in temp_nodes in DFS order
		void construct_cann_tree(const std::vector<size_type>& C, std::vector<pc_node_ex>& temp_nodes, std::map<uint64_t ,sdsl::bit_vector>& chars_coding,
				std::vector<uint64_t>& chars_pruned_length) {
			tMPQPII pq;
			size_type i = 0;
			// add leaves of Huffman tree
			std::for_each(std::begin(C), std::end(C), [&](decltype(*std::begin(C)) &freq) {
				if (freq > 0) {
							pq.push(tPII(freq, C.size() - i, temp_nodes.size()));// push (frequency, node pointer)
							// initial bv_pos with number of occurrences and bv_pos_rank
							// value with the code of the corresponding char, parent,
							// child[0], and child[1] are set to undef
							temp_nodes.emplace_back( pc_node_ex(freq, i));
				}
				++i;
			});
			while (pq.size() > 1) {
				tPII v1, v2;
				v1 = pq.top(); pq.pop();
				v2 = pq.top(); pq.pop();
				temp_nodes[std::get<2>(v1)].parent = temp_nodes.size(); // parent is new node
				temp_nodes[std::get<2>(v2)].parent = temp_nodes.size(); // parent is new node
				size_type frq_sum = std::get<0>(v1) + std::get<0>(v2);
				pq.push(tPII(frq_sum, -1, temp_nodes.size()));
				temp_nodes.emplace_back(pc_node_ex(frq_sum, 0, pc_node::undef,
						std::get<2>(v1), std::get<2>(v2)));
			}

			std::vector<std::tuple<uint64_t, uint64_t, uint64_t>> depthFrequencyAndSymbol;
			for( size_type j = 0; j < temp_nodes.size(); j++)
			{
				if(temp_nodes[j].child[0] == pc_node::undef)
				{
					size_type depth = 0;
					size_type current = j;
					while(temp_nodes[current].parent != pc_node::undef)
					{
						depth++;
						current = temp_nodes[current].parent;
					}
					depthFrequencyAndSymbol.push_back(std::tuple<uint64_t, uint64_t,  uint64_t>(depth, ((uint64_t)-1) - temp_nodes[j].freq, temp_nodes[j].sym));
				}
			}

			temp_nodes.clear();
			size_type currentNode = 0;
			bit_vector bv;

			std::sort(depthFrequencyAndSymbol.begin(), depthFrequencyAndSymbol.end());
			/*for(std::vector<std::tuple<uint64_t, uint64_t, uint64_t>>::iterator it = depthFrequencyAndSymbol.begin(); it != depthFrequencyAndSymbol.end(); it++)
						{
				std::cout << "depth " << std::get<0>(*it) << " freq is " << ((uint64_t)-1) - std::get<1>(*it) << " sym is " <<(char)std::get<2>(*it) <<std::endl;
						}*/
			temp_nodes.emplace_back(pc_node_ex(0, 0));
			size_type nodesToGoDown = 0;
			for(std::vector<std::tuple<uint64_t, uint64_t, uint64_t>>::iterator it = depthFrequencyAndSymbol.begin(); it != depthFrequencyAndSymbol.end(); it++)
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
						temp_nodes.emplace_back(pc_node_ex(0, 0, currentNode));
						currentNode++;

					}
				}

				//update frequencies of ancestors
				size_type tempNode = currentNode;
				while(tempNode != pc_node::undef)
				{
					temp_nodes[tempNode].freq += ((uint64_t)-1) - std::get<1>(*it);
					tempNode = temp_nodes[tempNode].parent;
				}
				temp_nodes[currentNode].sym = std::get<2>(*it);

				chars_coding[std::get<2>(*it)] = bv;
				if(chars_pruned_length.size() <= std::get<2>(*it))
				{
					chars_pruned_length.resize(std::get<2>(*it) + 1);
				}
				chars_pruned_length[std::get<2>(*it)] = 0;
				temp_nodes[currentNode].minimal_depth = temp_nodes[currentNode].maximal_depth = 0;

if(bv.size() > 0)
{
				//adding 1. going up if needed
				size_type temp = bv.size() - 1;
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
					temp_nodes.emplace_back(pc_node_ex(0, 0, currentNode));
					currentNode = temp_nodes.size() - 1;
					bv[temp] = 1;
				}
}
			}


			pc_node_ex root = temp_nodes[0];
			root.minimal_depth = temp_nodes[0].minimal_depth;
			root.maximal_depth = temp_nodes[0].maximal_depth;
			for (size_type i = 1; i < temp_nodes.size(); ++i) {
				temp_nodes[i-1] = temp_nodes[i];
				temp_nodes[i-1].parent = (temp_nodes[i-1].parent+temp_nodes.size()-1)%temp_nodes.size();
				temp_nodes[i-1].child[0] -= (temp_nodes[i-1].child[0] != pc_node::undef);
				temp_nodes[i-1].child[1] -= (temp_nodes[i-1].child[1] != pc_node::undef);
				temp_nodes[i-1].minimal_depth = temp_nodes[i].minimal_depth;
				temp_nodes[i-1].maximal_depth = temp_nodes[i].maximal_depth;
			}
			root.child[0] -= (root.child[0] != pc_node::undef);
			root.child[1] -= (root.child[1] != pc_node::undef);
			temp_nodes[temp_nodes.size()-1] = root;
			temp_nodes[temp_nodes.size()-1].minimal_depth = root.minimal_depth;
			temp_nodes[temp_nodes.size()-1].maximal_depth = root.maximal_depth;
	/*
			//printing the tree
			for(int i = 0; i < temp_nodes.size(); i++)
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

		//prunes full sub trees f the canonical tree
		void prune(std::vector<pc_node_ex>& temp_nodes, std::vector<uint64_t>& chars_pruned_length)
		{
			//m_h.resize(temp_nodes.size(), 0);
			std::vector<pc_node_ex> pruned_nodes(temp_nodes.begin(), temp_nodes.end());
			size_type delta = 0;
			//std::cout << temp_nodes.size() << " temp nodes size" <<std::endl;
			//std::cout.flush();

			pc_node_ex root = temp_nodes[temp_nodes.size() - 1] ;

			//checking if is a full tree
			if(root.minimal_depth == root.maximal_depth && root.minimal_depth > 0)
			{
				pruned_nodes.erase(pruned_nodes.begin(),  pruned_nodes.begin() + temp_nodes.size() - 1);
				//update pruned length of leaves under current node
				for(size_type i = 0; i < temp_nodes.size() - 1; i++)
				{
					if(temp_nodes[i].child[0] == pc_node::undef)
					{
						//chars_pruned_length.insert(std::pair<uint64_t, uint64_t>(temp_nodes[i].sym, root.minimal_depth));
						if(chars_pruned_length.size() < temp_nodes[i].sym + 1)
						{
							chars_pruned_length.resize(temp_nodes[i].sym + 1);
						}
						chars_pruned_length[temp_nodes[i].sym] = root.minimal_depth;
					}
				}
				pruned_nodes[0].child[0] = pc_node::undef;
				pruned_nodes[0].child[1] = root.minimal_depth;
				pruned_nodes[0].sym = 0;
			}
			else
			{
				for(size_type currentNode = 0; currentNode < temp_nodes.size() - 1; currentNode++)
				{
					if(temp_nodes[currentNode].minimal_depth == temp_nodes[currentNode].maximal_depth
							&& temp_nodes[currentNode].minimal_depth != 0)
					{
					  /*	std::cout << "Node " << currentNode << " is a root of a full subtree of depth " <<
								temp_nodes[currentNode].minimal_depth;
						std::cout << std::endl;
					  */
						//pruning the siblings of current node
						size_type current_delta = pow(2., (double)temp_nodes[currentNode].minimal_depth + 1) - 2;
						//std::cout << "current delta " <<current_delta <<std::endl;

						pruned_nodes[currentNode - delta].child[1] = temp_nodes[currentNode].minimal_depth;

						//std::cout << "pruning from " << currentNode + 1 - delta << " to "<<currentNode + 1 + current_delta - delta <<std::endl;
						pruned_nodes.erase(pruned_nodes.begin() + currentNode + 1 - delta,  pruned_nodes.begin() + currentNode + 1 + current_delta - delta);

						//update pruned length of leaves under current node
						for(size_type i = currentNode + 1; i < currentNode + 1 + current_delta; i++)
						{
							if(temp_nodes[i].child[0] == pc_node::undef)
							{
								if(chars_pruned_length.size() < temp_nodes[i].sym + 1)
								{
									chars_pruned_length.resize(temp_nodes[i].sym + 1);
								}
								//chars_pruned_length.insert(std::pair<uint64_t, uint64_t>(temp_nodes[i].sym, temp_nodes[currentNode].minimal_depth));
								chars_pruned_length[temp_nodes[i].sym] = temp_nodes[currentNode].minimal_depth;
							}
						}

						//uint64_t s = pruned_nodes.size();
						//std::cout << s <<std::endl;
						for(size_type i =0 ; i< pruned_nodes.size() ; i++)
						{
							if(pruned_nodes[i].child[0] >= currentNode + 1 - delta && pruned_nodes[i].child[0] != pc_node::undef)
							{
								pruned_nodes[i].child[0] -= current_delta;
							}
							if(pruned_nodes[i].child[1] >= currentNode + 1 - delta && pruned_nodes[i].child[0] != pc_node::undef)//intentionally checking child[0] because child[1] might represent the depth
							{
								pruned_nodes[i].child[1] -= current_delta;
							}
							if(pruned_nodes[i].parent >= currentNode + 1 - delta && pruned_nodes[i].parent != pc_node::undef)
							{
								pruned_nodes[i].parent -= current_delta;
							}
						}
						pruned_nodes[currentNode - delta].child[0] = pc_node::undef;
						pruned_nodes[currentNode - delta].child[1] = pruned_nodes[currentNode - delta].maximal_depth;
						//pruned_nodes[currentNode - delta].minimal_depth = pruned_nodes[current_child - delta].maximal_depth = 0;
						pruned_nodes[currentNode - delta].sym = 0;

						delta += current_delta;

						currentNode += current_delta; // +1 will come from the loop iteration
					}
				}
			}
			std::swap(pruned_nodes, temp_nodes);
			/*//printing the tree
			std::cout<< "full" <<std::endl;
			for(int i = 0; i < pruned_nodes.size(); i++)
			{
				char symbol;
				if(pruned_nodes[i].sym)
					symbol = (char)pruned_nodes[i].sym;
				else
					symbol = '*';
				std::cout << i << ": symbol: " << symbol << ", parent is " << pruned_nodes[i].parent
						<< ", left child is " << pruned_nodes[i].child[0]  << ", right child is " << pruned_nodes[i].child[1]
						<< ", freq is " << pruned_nodes[i].freq
						<< ", minimal depth is " << pruned_nodes[i].minimal_depth
						<< ", maximal depth is " << pruned_nodes[i].maximal_depth<< std::endl;
			}
			std::cout<< "pruned" <<std::endl;
			for(int i = 0; i < temp_nodes.size(); i++)
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

		//! Swap operator
		void swap(wt_new_sk& wt) {
			if (this != &wt) {
				std::swap(this->m_size, wt.m_size);
				std::swap(this->m_sigma,  wt.m_sigma);
				this->m_bv.swap(wt.m_bv);
				util::swap_support(this->m_bv_rank, wt.m_bv_rank,
								   &this->m_bv, &(wt.m_bv));

				util::swap_support(this->m_bv_select1, wt.m_bv_select1,
								   &this->m_bv, &(wt.m_bv));
				util::swap_support(this->m_bv_select0, wt.m_bv_select0,
								   &this->m_bv, &(wt.m_bv));
				this->m_tree.swap(wt.m_tree);

				//std::swap(m_chars_coding, wt.m_chars_coding);
				std::swap(m_path_to_char, wt.m_path_to_char);
				std::swap(m_leaves_bv, wt.m_leaves_bv);
			}
		}

		// calculates the tree shape returns the size of the WT bit vector
		virtual size_type construct_tree_shape(const std::vector<size_type>& C, std::map<uint64_t ,sdsl::bit_vector>& chars_coding, std::vector<uint64_t>& chars_pruned_length,
				std::vector<pc_node_ex>& temp_nodes) {
			// vector  for node of the tree
			construct_cann_tree(C, temp_nodes, chars_coding, chars_pruned_length);
			prune(temp_nodes, chars_pruned_length);
			// Convert code tree into BFS order in memory and
			// calculate bv_pos values
			size_type bv_size = 0;
			std::vector<pc_node> temp(temp_nodes.begin(), temp_nodes.end());
			tree_strat_type temp_tree(temp, bv_size, this);
			this->m_tree.swap(temp_tree);

			return bv_size;
		}

		//creates the tree and the bit_vector, and initializing rand & select
		wt_new_sk(int_vector_buffer<tree_strat_type::int_width>& input_buf,
				size_type size) {
			this->m_size = size;
		    if (0 == this->m_size)
		        return;

		    // O(n + |\Sigma|\log|\Sigma|) algorithm for calculating node sizes
		    // TODO: C should also depend on the tree_strategy. C is just a mapping
		    // from a symbol to its frequency. So a map<uint64_t,uint64_t> could be
		    // used for integer alphabets...
		    std::vector<size_type> C;
		    // 1. Count occurrences of characters
		    calculate_character_occurences(input_buf, this->m_size, C);
		    // 2. Calculate effective alphabet size
		    calculate_effective_alphabet_size(C, this->m_sigma);
		    // 3. Generate tree shape
		    std::map<uint64_t ,bit_vector> chars_coding;
		    std::vector<uint64_t> chars_pruned_length;
		    std::vector<pc_node_ex> temp_nodes;
		    size_type tree_size = this->construct_tree_shape(C, chars_coding, chars_pruned_length, temp_nodes);
		    // 4. Generate wavelet tree bit sequence m_bv
		    std::vector<uint64_t> bv_node_pos(this->m_tree.size(), 0);


		    //create BFS\DFS mapping [temp_nodes and m_tree]
		    node_type node_in_tree = this->m_tree.child(this->m_tree.root(),0);
		    size_type node_in_temp = 0;
		    std::map<size_type, node_type> temp_to_tree;
		    temp_to_tree[temp_nodes.size() - 1] = 0;
		    while(node_in_temp != temp_nodes.size() - 1) // check we didn't reach the root.
		    {
		    	temp_to_tree[node_in_temp] = node_in_tree;

		    	if(temp_nodes[node_in_temp].minimal_depth > 0 && temp_nodes[node_in_temp].child[0] != pc_node::undef)
		    	{
		    		node_in_temp = temp_nodes[node_in_temp].child[0];
					node_in_tree = this->m_tree.child(node_in_tree, 0);
		    	}
		    	else
		    	{
		    		size_type parent;
					bool is_left_child;
					do
					{
						parent = temp_nodes[node_in_temp].parent;
						is_left_child = parent != pc_node::undef && temp_nodes[parent].child[0] == node_in_temp;
						if(is_left_child)
						{
							node_in_temp = temp_nodes[parent].child[1];
							node_in_tree = this->m_tree.child(this->m_tree.parent(node_in_tree), 1);
						}
						else if (parent != pc_node::undef)
						{
							node_in_temp = parent;
							node_in_tree = this->m_tree.parent(node_in_tree);
						}
					}
					while(!is_left_child && parent != pc_node::undef);
		    	}
		    }

		    size_type leaves_bv_size = 0;
		    std::vector<size_type> h(this->m_tree.size());
		    for(size_type node_in_temp = 0; node_in_temp < temp_nodes.size(); node_in_temp++)
		    {
		    	if(temp_nodes[node_in_temp].minimal_depth == temp_nodes[node_in_temp].maximal_depth && temp_nodes[node_in_temp].minimal_depth > 0)
				{
					//this node is a root of a full subtree
		    		size_type current_extra = (temp_nodes[node_in_temp].minimal_depth) * temp_nodes[node_in_temp].freq;
					 //uint64_t current_saving = temp_nodes[node_in_temp].freq;//not needed because was not generated to leaves
					 //tree_size -= current_saving; //not needed because was not generated to leaves
					 leaves_bv_size += current_extra;
					 for(int node_in_tree = temp_to_tree[node_in_temp] + 1; node_in_tree < this->m_tree.size() ; node_in_tree++)
					 {
						 if(this->m_tree.is_leaf(node_in_tree))
						 {
							 bv_node_pos[node_in_tree] += current_extra;
						 }
						 else
						 {
							 //bv_node_pos[node_in_tree] -= current_saving; //not needed because was not generated to leaves
						 }
					 }
				}
		    }

		    //std::swap(h, m_h);

		    bit_vector temp_bv(tree_size, 0);
		    bit_vector temp_leaves_bv(leaves_bv_size, 0);

		    // Initializing starting position of wavelet tree nodes

		    for (size_type v=0; v < this->m_tree.size(); ++v) {
		    	if(!this->m_tree.is_leaf(v))
		    	{
		    		bv_node_pos[v] += this->m_tree.bv_pos(v);
		    	}
		        this->m_tree.set_bv_pos(v, bv_node_pos[v]);
		        //std::cout << v << " bv_node_pos is " << bv_node_pos[v] <<std::endl;
		    }

		    if (input_buf.size() < size) {
		        throw std::logic_error("Stream size is smaller than size!");
		        return;
		    }

		    for(auto it = chars_coding.begin(); it != chars_coding.end(); it++)
			{
		    	auto size = it->second.size();
				auto path = it->second.get_int(0, size);

				this->m_tree.set_bit_path(it->first, path, size);
				path = bits::rev(path);
				path = path >> (64-size); // remove the length
				path = path | (size << 56);
				m_path_to_char[path] = it->first;
			}

		    value_type old_chr = input_buf[0];
		    size_type times = 0;
		    for (size_type i=0; i < this->m_size; ++i) {
		        value_type chr = input_buf[i];
		        if (chr != old_chr) {
		        	this->insert_char(old_chr, bv_node_pos, times, temp_bv, temp_leaves_bv, chars_pruned_length);
		            times = 1;
		            old_chr = chr;
		        } else { // chr == old_chr
		            ++times;
		            if (times == 64) {
		            	this->insert_char(old_chr, bv_node_pos, times, temp_bv, temp_leaves_bv, chars_pruned_length);
		                times = 0;
		            }
		        }
		    }
		    if (times > 0) {
		    	this->insert_char(old_chr, bv_node_pos, times, temp_bv, temp_leaves_bv, chars_pruned_length);
		    }
		    this->m_bv = bit_vector_type(std::move(temp_bv));
		    std::swap(this->m_leaves_bv, temp_leaves_bv);

		    /*int j = 0;
		    for(int i =0 ; i < this->m_bv.size(); i++)
		    {
		    	while(this->m_tree.is_leaf(j))
		    	{
		    		j++;
		    	}
		    	while(i == bv_node_pos[j])
		    	{
		    		std::cout << std::endl;

		    		j++;
		    	}
		    	std::cout << this->m_bv[i];
		    }
		    std::cout << std::endl;

		    j = 0;
			for(int i =0 ; i < this->m_leaves_bv.size(); i++)
			{
				while(!this->m_tree.is_leaf(j))
				{
					j++;
				}
				while(i == bv_node_pos[j])
				{
					std::cout << std::endl;

					j++;
				}
				std::cout << this->m_leaves_bv[i];
			}
			std::cout << std::endl;*/
		    // 5. Initialize rank and select data structures for m_bv
		    this->construct_init_rank_select();
		    // 6. Finish inner nodes by precalculating the bv_pos_rank values
		    this->m_tree.init_node_ranks(this->m_bv_rank);

		    //std::swap(chars_coding, m_chars_coding);
		    //m_chars_coding = chars_coding;

		}
		wt_new_sk(){}


		// insert a character into the wavelet tree
		virtual void insert_char(value_type old_chr, std::vector<size_type>& bv_node_pos,
						 size_type times, bit_vector& bv, bit_vector& leaves_bv,
						 std::vector<size_type>& chars_pruned_length) {
			//bit_vector char_code = chars_coding[old_chr];
			size_type pruned_length = chars_pruned_length[old_chr];
			//auto code = chars_coding.find(old_chr);
			node_type v = this->m_tree.root();

			//auto path_pair = this->path(old_chr);
			//uint64_t path_len =  path_pair.first;
			//uint64_t p = bits::rev(path_pair.second) >> (64 - path_len);
			uint64_t p = this->m_tree.bit_path(old_chr);
			uint64_t path_len = p >> 56;

			for (uint32_t l=0; l< path_len - pruned_length ; ++l, p >>= 1) {
				if (p&1) {
					bv.set_int(bv_node_pos[v], 0xFFFFFFFFFFFFFFFFULL,times);
				}
				bv_node_pos[v] += times;
				v = this->m_tree.child(v, p&1);
			}


			if(pruned_length > 0)
			{
				uint64_t p_rev = bits::rev(p) >> (64 - pruned_length);

				for(int j = 0; j < times; j++)
				{
					leaves_bv.set_int(bv_node_pos[v], p_rev, pruned_length);
					bv_node_pos[v] += pruned_length;
				}
			}
		}

        //! Recovers the i-th symbol of the original vector.
        /*!
         * \param i Index in the original vector.
         * \return The i-th symbol of the original vector.
         * \par Time complexity
         *      \f$ \Order{H_0} \f$ on average, where \f$ H_0 \f$ is the
         *      zero order entropy of the sequence
         *
         * \par Precondition
         *      \f$ i < size() \f$
         */
        value_type operator[](size_type i)const {
            assert(i < this->size());
            // which stores how many of the next symbols are equal
            // with the current char
            uint64_t path =0;
            uint64_t path_len = 0;
            node_type v = this->m_tree.root(); // start at root node
            while (!this->m_tree.is_leaf(v)) {   // while  not a leaf
            	//current_bv.resize(current_bv.size() + 1);
            	path = path << 1;
            	path_len++;
                if (this->m_bv[ this->m_tree.bv_pos(v) + i]) {  // goto right child
                    i = this->m_bv_rank(this->m_tree.bv_pos(v) + i)
                        - this->m_tree.bv_pos_rank(v);
                    v = this->m_tree.child(v,1);
                    path = path | 1;
                } else { // goto the left child
                    i -= (this->m_bv_rank(this->m_tree.bv_pos(v) + i)
                          - this->m_tree.bv_pos_rank(v));
                    v = this->m_tree.child(v,0);
                }
            }

            //add the pruned bits
            if(this->m_tree.child(v,1) != (node_type)-1)
            {
            	size_type h = this->m_tree.child(v,1);
            	path_len += h;

            	uint64_t pruned = m_leaves_bv.get_int(this->m_tree.bv_pos(v) + h * i, h);
            	path = path << h | pruned;
            }
            path |= (path_len << 56);
            auto f = m_path_to_char.find(path);

            if(f != m_path_to_char.end())
            	return f->second;
            std::cout << "failed to find needed path:" << std::bitset<64>(path).to_string() << std::endl;
            for(auto it = m_path_to_char.begin(); it != m_path_to_char.end(); it++)
            {
            	std::cout<< "path: " <<std::bitset<64>(it->first).to_string() << ".\tsymbol: " <<it->second <<".\tsymbol bits: " <<std::bitset<64>(it->second).to_string() <<std::endl;
            }
            std::cout.flush();
            return 0;
			/*//todo: improve
			for (std::map<uint64_t ,sdsl::bit_vector>::const_iterator it = m_chars_coding.begin(); it != m_chars_coding.end(); ++it )
				if (it->second == current_bv)
					return it->first;

			// if v is an original leaf bv_pos_rank returns symbol itself
			return this->m_tree.bv_pos_rank(v);*/

        }

        //! Serializes the data structure into the given ostream
		virtual size_type serialize(std::ostream& out, structure_tree_node* v=nullptr,
							std::string name="") const {

			structure_tree_node* child = structure_tree::add_child(
											 v, name, util::class_name(*this));
			size_type written_bytes = 0;
			written_bytes += write_member(this->m_size,out,child, "size");
			written_bytes += write_member(this->m_sigma,out,child, "sigma");
			written_bytes += this->m_bv.serialize(out,child,"bv");
			written_bytes += this->m_bv_rank.serialize(out,child,"bv_rank");
			written_bytes += this->m_bv_select1.serialize(out,child,"bv_select_1");
			written_bytes += this->m_bv_select0.serialize(out,child,"bv_select_0");
			written_bytes += this->m_tree.serialize(out,child,"tree");
			written_bytes += this->m_leaves_bv.serialize(out,child,"leaves_bv");

			//written_bytes += write_member(m_path_to_char, out, child, "path_to_char");
			//todo: add it to the count in the graphics also
			//std::map<uint64_t ,sdsl::bit_vector>::iterator
			for(auto it = m_path_to_char.begin(); it != m_path_to_char.end(); it++)
			{
//				out.write((char*)&it->first, sizeof(it->first));
				written_bytes += write_member(it->first,out,child, "path");
				written_bytes += write_member(it->second,out,child, "char");
			}

			structure_tree::add_size(child, written_bytes);
			return written_bytes;
		}

		//! Loads the data structure from the given istream.
		virtual void load(std::istream& in) {
			read_member(this->m_size, in);
			read_member(this->m_sigma, in);
			this->m_bv.load(in);
			this->m_bv_rank.load(in, &this->m_bv);
			this->m_bv_select1.load(in, &this->m_bv);
			this->m_bv_select0.load(in, &this->m_bv);
			this->m_tree.load(in);
			this->m_leaves_bv.load(in);
			//m_path_to_char.load(in);
/*			uint64_t m_path_to_char_size = 0;
			read_member(m_path_to_char_size, in);
			m_path_to_char = std::vector<uint64_t>(m_path_to_char_size);
			load_vector(m_path_to_char, in);*/
			for(int i = 0; i< this->m_sigma; i++)
			{
				size_type symbol, path;

				read_member(path, in);
				read_member(symbol, in);
				m_path_to_char[path] = symbol;
			}
		}
		;
	};
}
#endif /* NEW_SKELETONWAVELETTREE_H_ */







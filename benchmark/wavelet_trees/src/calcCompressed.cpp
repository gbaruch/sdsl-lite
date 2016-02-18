/*
 * calcCompressed.cpp
 *
 *  Created on: Jun 17, 2015
 *      Author: gilad
 */

#include<sdsl/wavelet_trees.hpp>
#include<sdsl/wt_helper.hpp>

#include <iostream>
using namespace sdsl;
using namespace std;

int main(int argc, char* argv[])
{
    if (argc < 4) {
        std::cout << "Usage: file num_bytes testcase_name wt_id" << endl;
        return 1;
    }
    uint8_t type = argv[2][0]=='d' ? 'd' : argv[2][0]-'0';

    wt_huff<> huff;
    construct(huff, argv[1], type);

    cout << "# COMPRESSED_SIZE = " << huff.m_bv.size() / 8 << endl;
    return 0;
}

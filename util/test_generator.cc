
// Copyright (c) Simon Fraser University & The Chinese University of Hong Kong. All rights reserved.
// Licensed under the MIT license.
/*
* Used to test the key generator (different key distribution)
*/
#include <random>
#include <cstdint>
#include "key_generator.hpp"
#include <map>
#include <ctime>
#include <iostream>

int main(){
	uint64_t N = 100;
	uint64_t key_num = 100;
	key_generator_t *kg1;
	key_generator_t *kg2;
	kg1 = new zipfian_key_generator_t(100, 0.8);
	kg2 = new uniform_key_generator_t();

	std::map<uint64_t,int> counting_;
	std::map<uint64_t,int> counting2_;

	for (int i = 0; i < key_num; ++i){
		uint64_t random_num = kg1->next_uint64();
		if(counting_.find(random_num) != counting_.end()){
			counting_[random_num] += 1;
		}else{
			counting_.insert(std::pair<uint64_t, int>(random_num, 1));
		}
		
		random_num = kg2->next_uint64();
		if(counting2_.find(random_num) != counting_.end()){
			counting2_[random_num] += 1;
		}else{
			counting2_.insert(std::pair<uint64_t, int>(random_num, 1));
		}

	}
	uint64_t count = 0;
	uint64_t single = 0;
	for(auto iter = counting_.begin(); iter != counting_.end(); iter++){
		single++;
		std::cout << iter->first <<  " : " << iter->second << std::endl;
		count += iter->second;
	}

	std::cout << "SUM = " << count << std::endl;
	std::cout << "Single = " << single << std::endl;
	count = 0;
	for(auto iter = counting2_.begin(); iter != counting2_.end(); iter++){
		std::cout << iter->first <<  " : " << iter->second << std::endl;
		count += iter->second;
	}
	std::cout << "SUM = " << count << std::endl;	
}

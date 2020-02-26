#pragma once

#include "zipfian_int_distribution.hpp"
#include "random.h"

class key_generator_t{
public:
 key_generator_t(){

 }

 virtual uint64_t next_uint64() = 0;
};

class uniform_key_generator_t : public key_generator_t{
public:
    uniform_key_generator_t(){
        unsigned long long init[4]={0x12345ULL, 0x23456ULL, 0x34567ULL, 0x45678ULL}, length=4;
        init_by_array64(init, length);
    }

    uint64_t next_uint64(){
        return genrand64_int64();
    }
};

class range_key_generator_t : public key_generator_t{
public:
    range_key_generator_t(size_t start){
        start_number = start;
    }

    uint64_t next_uint64(){
        return start_number++;
    }
private:
    size_t start_number = 0;
};

class zipfian_key_generator_t : public key_generator_t{
public:
    zipfian_key_generator_t(size_t start, size_t end, float skew) : dist_(start, end, skew), generator_(time(0)){

    }

    uint64_t next_uint64(){
        return dist_(generator_);
    }

private:
    std::default_random_engine generator_;
    zipfian_int_distribution<uint64_t> dist_;
};
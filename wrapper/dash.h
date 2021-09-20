#ifndef __DASH_WRAPPER_HPP__
#define __DASH_WRAPPER_HPP__

#include "tree_api.hpp"
#include <thread>
#include "../src/ex_finger.h"
#include "../src/allocator.h"

static const char *wrapper_path = "/mnt/pmem0/dash.data";

class dash_wrapper : public tree_api
{
public:
    dash_wrapper() { 
        // Create an instance of Dash
        bool file_exist = false;
        if (FileExists(wrapper_path)) {
            file_exist = true;
            std::cout << "Dash pool already exists, please remove it for testing" << std::endl;
        }

        Allocator::Initialize(wrapper_path, 1024ul * 1024ul * 1024ul * 10ul);
        dash_instance_ = reinterpret_cast<extendible::Finger_EH<uint64_t> *>(
            Allocator::GetRoot(sizeof(extendible::Finger_EH<uint64_t>)));

        if (!file_exist) {
            // During initialization phase, allocate 64 segments for Dash-EH
            size_t segment_number = 64;
            new (dash_instance_) extendible::Finger_EH<uint64_t>(
                segment_number, Allocator::Get()->pm_pool_);
        }else{
            new (dash_instance_) extendible::Finger_EH<uint64_t>();
        }
    }
    
    virtual ~dash_wrapper() { }
    
    virtual bool find(const char* key, size_t key_sz, char* value_out) override
    {
        // For now only test 8 bytes key and value (uint64_t)
        thread_local int insert_counter(0);
        if(insert_counter == 0){
            // Enroll into the epoch
            Allocator::Protect();
        }
        insert_counter = (insert_counter + 1) & ((1 << 10) - 1);

        auto value = dash_instance_->Get(*reinterpret_cast<uint64_t*>(const_cast<char*>(key)));

        if(insert_counter == 1023){
            // Exit the epoch
            Allocator::Unprotect();
        }
        if (value == nullptr)
            return false;
        memcpy(value_out, &value, sizeof(value));
        return true;
    }

    virtual bool insert(const char* key, size_t key_sz, const char* value, size_t value_sz) override
    {
        thread_local int insert_counter(0);
        if(insert_counter == 0){
            // Enroll into the epoch
            Allocator::Protect();
        }
        insert_counter = (insert_counter + 1) & ((1 << 10) - 1);

        auto K = *reinterpret_cast<uint64_t*>(const_cast<char*>(key));
        auto ret = dash_instance_->Insert(K, value);
        
        if(insert_counter == 1023){
            // Exit the epoch
            Allocator::Unprotect();
        }
        return ((ret == 0) ? true : false);
    }

    virtual bool update(const char* key, size_t key_sz, const char* value, size_t value_sz) override
    {
        // Dash does not provides this interface, although easy to implement        
        return true;
    }

    virtual bool remove(const char* key, size_t key_sz) override
    {
        thread_local int insert_counter(0);
        if(insert_counter == 0){
            // Enroll into the epoch
            Allocator::Protect();
        }
        insert_counter = (insert_counter + 1) & ((1 << 10) - 1);

        auto ret = dash_instance_->Delete(*reinterpret_cast<uint64_t*>(const_cast<char*>(key)));
        
        if(insert_counter == 1023){
            // Exit the epoch
            Allocator::Unprotect();
        }
        return ret;
    }

    virtual int scan(const char* key, size_t key_sz, int scan_sz, char*& values_out) override
    {   
        // Dash does not support scan
        return 0;
    }
private:
    extendible::Finger_EH<uint64_t>* dash_instance_;
};

#endif
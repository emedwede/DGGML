#ifndef DGGML_HISTOBUCKET_HPP
#define DGGML_HISTOBUCKET_HPP

#include <iostream> 
#include "Utlities/MemoryManager.hpp"

namespace DGGML
{

template <typename KeyType>
struct Histobucket 
{
    using count_type = std::size_t;
    using key_type = KeyType;

    Histobucket() : num_bins(0) {
        counts = DGGML::MemoryManager::allocate_std<count_type>(0);
        offsets = DGGML::MemoryManager::allocate_std<count_type>(0);
        capacities = DGGML::MemoryManager::allocate_std<count_type>(0);
        bucket_data = DGGML::MemoryManager::allocate_std<key_type>(0);
    }
    
    Histobucket(std::size_t n) : num_bins(n) {
        counts = DGGML::MemoryManager::allocate_std<count_type>(n);
        offsets = DGGML::MemoryManager::allocate_std<count_type>(n);
        capacities = DGGML::MemoryManager::allocate_std<count_type>(n);
        bucket_data = DGGML::MemoryManager::allocate_std<key_type>(n);

        zero_init();
    }
    
    void reset_and_realloc(std::size_t n)
    {
        DGGML::MemoryManager::deallocate_std(counts);
        DGGML::MemoryManager::deallocate_std(offsets);
        DGGML::MemoryManager::deallocate_std(capacities);

        num_bins = n;
        counts = DGGML::MemoryManager::allocate_std<count_type>(n);
        offsets = DGGML::MemoryManager::allocate_std<count_type>(n);
        capacities = DGGML::MemoryManager::allocate_std<count_type>(n);
         
        zero_init();

    }

    ~Histobucket() 
    {
        DGGML::MemoryManager::deallocate_std(counts);
        DGGML::MemoryManager::deallocate_std(offsets);
        DGGML::MemoryManager::deallocate_std(capacities);
        DGGML::MemoryManager::deallocate_std(bucket_data);
    }
   
    count_type numBin() const 
    {
        return num_bins;
    }

    count_type binSize(const count_type bin_id) const 
    {
        return counts[bin_id];
    }
    
    count_type maxSize()
    {
        std::size_t max_size = binSize(0);                                                                                                                                        
        for(auto i = 1; i < numBin(); i++)                                                                                                                                        
        {                                                                                                                                                                                     
            if(binSize(i) > max_size) max_size = binSize(i);                                                                            
        } 
        return max_size;
    }
    
    count_type totalSize()
    {
        std::size_t total = 0;
        for(auto i = 0; i < numBin(); i++)
        {
            total += binSize(i);
        }
        return total;
    }
    void incrementBin(const count_type bin_id) const 
    {
        counts[bin_id]++;
    }
    
    void reset_counts() const
    {
        for (auto i = 0; i < num_bins; i++)
        {
            counts[i] = 0;
        }
    }

    std::size_t totalBinCounts() const 
    {
        std::size_t sum = 0;
        for(auto i = 0; i < num_bins; i++)
        {
            sum += binSize(i);
        }
        return sum;
    }

    std::size_t totalBinCapacity() const 
    {
        std::size_t c = 0;
        for(auto i = 0; i < num_bins; i++)
        {
            c += binCapacity(i);
        }
        return c;
    }
    
    count_type binBegin(const count_type bin_id) const
    {
        return offsets[bin_id];
    }
    
    count_type binEnd(const count_type bin_id)
    {
        return offsets[bin_id] + counts[bin_id];
    }

    count_type binCapacity(const count_type bin_id) const
    {
        return capacities[bin_id];
    }
    
    void setBinCapacity(count_type c, const count_type bin_id)
    {
        capacities[bin_id] = c;
    }

    void setBinOffset(count_type o, const count_type bin_id)
    {
        offsets[bin_id] = o;
    }

    void bin_insert(key_type key, const count_type bin_id) 
    {
        if(binSize(bin_id) < binCapacity(bin_id)) //then there is room to insert
        {
            bucket_data[binBegin(bin_id)+binSize(bin_id)] = key;
            incrementBin(bin_id);
        }
    }
    
    void build_buckets()
    {
        DGGML::MemoryManager::deallocate_std<key_type>(bucket_data);
        DGGML::MemoryManager::allocate_std<key_type>(totalBinCapacity());
    }
    friend std::ostream& operator<<(std::ostream& os, Histobucket<KeyType>& bins)
    {
        os << "Number of bins: " << bins.num_bins << "\n";

        return os;
    }

    void zero_init()
    {
        for(auto i = 0; i < num_bins; i++)
        {
            counts[i] = 0;
            offsets[i] = 0;
            capacities[i] = 0;
        }
    }


    count_type num_bins; // how many bins are there
    count_type* counts;  // how many items are in a bucket
    count_type* offsets; // offsets for binning
    count_type* capacities; // extra space for the bucket
    key_type* bucket_data; //basically the edge relationships 
};

} // end namespace DGGML

#endif

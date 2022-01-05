#ifndef CAJETE_HISTOBUCKET_HPP
#define CAJETE_HISTOBUCKET_HPP 

#include <iostream> 
#include "MemoryManager.hpp"

namespace Cajete 
{

template <typename KeyType>
class Histobucket 
{
    public:
        using count_type = std::size_t;
        using key_type = KeyType;

        Histobucket() : num_bins(0) {
            counts = Cajete::MemoryManager::allocate_std<count_type>(0);
            offsets = Cajete::MemoryManager::allocate_std<count_type>(0);
            capacities = Cajete::MemoryManager::allocate_std<count_type>(0);
        }
        
        Histobucket(std::size_t n) : num_bins(n) {
            counts = Cajete::MemoryManager::allocate_std<count_type>(n);
            offsets = Cajete::MemoryManager::allocate_std<count_type>(n);
            capacities = Cajete::MemoryManager::allocate_std<count_type>(n);
            
            zero_init();
        }
        
        void reset_and_realloc(std::size_t n)
        {
            Cajete::MemoryManager::deallocate_std(counts);
            Cajete::MemoryManager::deallocate_std(offsets);
            Cajete::MemoryManager::deallocate_std(capacities);
            
            num_bins = n;
            counts = Cajete::MemoryManager::allocate_std<count_type>(n);
            offsets = Cajete::MemoryManager::allocate_std<count_type>(n);
            capacities = Cajete::MemoryManager::allocate_std<count_type>(n);
             
            zero_init();

        }

        ~Histobucket() 
        {
            Cajete::MemoryManager::deallocate_std(counts);
            Cajete::MemoryManager::deallocate_std(offsets);
            Cajete::MemoryManager::deallocate_std(capacities);
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

        void add_count(const count_type bin_id) const 
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
                add_count(bin_id);
            }
        }

        friend std::ostream& operator<<(std::ostream& os, Histobucket<KeyType>& bins)
        {
            os << "Number of bins: " << bins.num_bins << "\n";

            return os;
        }
    private:

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

} // end namespace Cajete 

#endif

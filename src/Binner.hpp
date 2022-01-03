#ifndef CAJETE_BINNER_HPP
#define CAJETE_BINNER_HPP 

#include <iostream> 
#include "MemoryManager.hpp"

namespace Cajete 
{

template <typename KeyType>
class Binner 
{
    public:
        using count_type = std::size_t;
        using index_type = KeyType;

        Binner() : num_bins(0) {
            counts = Cajete::MemoryManager::allocate_std<count_type>(0);
            offsets = Cajete::MemoryManager::allocate_std<count_type>(0);
            capacities = Cajete::MemoryManager::allocate_std<count_type>(0);
        }
        
        Binner(std::size_t n) : num_bins(n) {
            counts = Cajete::MemoryManager::allocate_std<count_type>(n);
            offsets = Cajete::MemoryManager::allocate_std<count_type>(n);
            capacities = Cajete::MemoryManager::allocate_std<count_type>(n);

            for(auto i = 0; i < n; i++)
            {
                counts[i] = 0;
                offsets[i] = 0;
                capacities[i] = 0;
            }
        }

        ~Binner() 
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

        count_type binOffset(const count_type bin_id) const 
        {
            return offsets[bin_id];
        }

        count_type binCapacity(const count_type bin_id) const
        {
            return capacities[bin_id];
        }


        friend std::ostream& operator<<(std::ostream& os, Binner<KeyType>& bins)
        {
            os << "Number of bins: " << bins.num_bins << "\n";

            return os;
        }
    private:
        count_type num_bins; // how many bins are there
        count_type* counts;  // how many items are in a bucket
        count_type* offsets; // offsets for binning
        count_type* capacities; // extra space for the bucket
};

} // end namespace Cajete 

#endif

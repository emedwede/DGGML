#ifndef CAJETE_MEMORYMANAGER_HPP
#define CAJETE_MEMORYMANAGER_HPP

#include <iostream> // TODO: replace with fmt library?

#include <memory> // For smart pointers 

// The memory Manager is part of the Library namespace
namespace Cajete 
{
// The memory manager has it's own namespace just to make
// it more distinct 
namespace MemoryManager
{
	// Standard host side non-smart style serial allocator
	template <typename T>
	T* allocate_std(std::size_t size) 
	{
		
		T* ptr;
	
	// These compiler statements enable flexibility
	// to add new allocators for different device spaces 
	#if defined CAJETE_ENABLE_HOST_SERIAL
			ptr = new T[size];
			#if defined(CAJETE_DEBUG)
				std::cout << "Successfully allocated " << sizeof(ptr)*size 
					<< " bytes at " << ptr << "\n";
			#endif
	//TODO: add new allocators here
	// #else / #elif
	#endif
		return ptr;
	}

	// Standard host side smart style serial allocator
	template <typename T>
	std::shared_ptr<T> allocate_smart(std::size_t size) {
			//Pre: C++ 17: 
			std::shared_ptr<T> ptr(new T[size], std::default_delete<T[]>());
			
			//C++ 17 Way:
			//std::shared_ptr<T> ptr(new T[size]);

			#if defined(CAJETE_DEBUG)
				std::cout << "Successfully allocated " << sizeof(ptr)*size 
					<< " bytes at " << ptr << "\n";
			#endif

			return ptr;
	}

	// Standard deallocator
	template <typename T>
	void deallocate_std(T *&ptr)
	{
		// if not null
		if(ptr) {
		// These compiler statements enable flexibility
		// to add new deallocators for different device spaces
		#if defined(CAJETE_ENABLE_HOST_SERIAL)
			delete [] ptr;
			#if defined(CAJETE_DEBUG) 
				std::cout << "Successfully deallocated bytes at " << ptr << "\n";
			#endif
		// TODO: add new deallocators here
		// #else / #elif
		#endif
			ptr = nullptr;
		}
	}
	
	//Smart pointers don't need deallocators

} //End namespace MemoryManager
} //End namespace Cajete 

#endif

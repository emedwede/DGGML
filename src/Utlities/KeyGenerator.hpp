#ifndef DGGML_KEYGENERATOR_HPP
#define DGGML_KEYGENERATOR_HPP

namespace DGGML
{

template<typename KeyType>
struct KeyGenerator {
    using key_type = KeyType;

    key_type current_key;

    KeyGenerator(KeyType key) : current_key(key) {}

    KeyGenerator() : current_key(0) {}

    key_type get_key() { return current_key++; }
};

}
#endif //DGGML_KEYGENERATOR_HPP

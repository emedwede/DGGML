#ifndef DGGML_RANDOMFUNCTIONS_HPP
#define DGGML_RANDOMFUNCTIONS_HPP

#include <random>

auto RandomRealsBetween(double low, double high)
{
    auto randomFunc =
            [distribution_ = std::uniform_real_distribution<double>(low, high),
                    random_engine_ = std::mt19937{std::random_device{}() }]() mutable
            {
                return distribution_(random_engine_);
            };
    return randomFunc;
};

auto RandomIntsBetween(int low, int high)
{
    auto randomFunc =
            [distribution_ = std::uniform_int_distribution<int>(low, high),
                    random_engine_ = std::mt19937{std::random_device{}() }]() mutable
            {
                return distribution_(random_engine_);
            };
    return randomFunc;
};

#endif //DGGML_RANDOMFUNCTIONS_HPP

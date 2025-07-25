#pragma once

#include <boost/unordered/unordered_flat_map.hpp>

inline namespace euler
{
/// Type for universal recursive memoizer.
template <typename Fun, typename... Args> struct memoize_t
{
    static_assert(std::is_invocable_v<Fun, std::reference_wrapper<memoize_t>, Args...>);

    using key_type = std::tuple<Args...>;
    using value_type = std::invoke_result_t<Fun, std::reference_wrapper<memoize_t>, Args...>;

    Fun f;
    boost::unordered_flat_map<key_type, value_type> cache;

    value_type operator()(Args... args)
    {
        key_type key{args...};
        if (auto const it = cache.find(key); it != cache.end())
            return it->second;
        auto res = f(std::ref(*this), std::move(args)...);
        cache.emplace(std::move(key), res);
        return res;
    }
};

/// Universal recursive memoizer.
/// Usage: `memoize<Args...>([](auto &&f, ...) -> T { ... })`.
template <typename... Args, typename Fun> memoize_t<Fun, Args...> memoize(Fun f) { return {std::move(f), {}}; }
} // namespace euler

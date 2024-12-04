#pragma once

#include <boost/unordered/unordered_flat_map.hpp>

inline namespace euler
{
/// Universal recursive memoizer.
/// Usage: `memoize([](auto &&f, ...) -> T { ... })`.
template <template <typename...> typename Map = boost::unordered_flat_map, typename Fun> constexpr auto memoize(Fun f)
{
    return [f = std::move(f)]<typename Self, typename... Args>(this Self &&self, Args &&...args) {
        static Map<std::tuple<std::decay_t<Args>...>, std::invoke_result_t<Fun, Self, Args...>> cache;
        if (auto it = cache.find({args...}); it != cache.end())
            return it->second;
        return cache[{args...}] = f(std::forward<Self>(self), std::forward<Args>(args)...);
    };
}
} // namespace euler

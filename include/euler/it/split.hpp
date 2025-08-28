#pragma once

#include "base.hpp"

inline namespace euler
{
namespace it
{
/// Splits a string at a character. TODO: make this work for wstrings too.
class split : public it_base
{
    std::string _s;
    char _d;

  public:
    using value_type = std::string;

    split(std::string s, char d) : _s(std::move(s)), _d(d) {}

    template <std::invocable<value_type> Fun> constexpr it::result_t operator()(Fun f) const
    {
        std::istringstream iss(_s);
        for (std::string token; std::getline(iss, token, _d);)
            if (!it::callbackResult(f, std::move(token)))
                return it::result_break;
        return it::result_continue;
    }
};
} // namespace it
} // namespace euler

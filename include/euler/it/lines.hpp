#pragma once

#include "base.hpp"
#include <fstream>

inline namespace euler
{
namespace it
{
/// Gets lines from a file.
class lines : public it_base
{
  public:
    using value_type = std::string;

    lines(std::string fileName) : _fileName(std::move(fileName)) {}

    template <std::invocable<value_type> Fun> constexpr it::result_t operator()(Fun f) const
    {
        std::ifstream fin(_fileName);
        for (std::string line; std::getline(fin, line);)
            if (!it::callbackResult(f, std::move(line)))
                return it::result_break;
        return it::result_continue;
    }

  private:
    std::string _fileName;
};
} // namespace it
} // namespace euler

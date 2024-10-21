#pragma once

#include "euler/matrix.hpp"

/// Convert from OKLCH to RGB.
/// @param l Lightness, between 0 and 1.
/// @param c Chroma, usually between 0 and 0.2 (but can go up to 0.37).
/// @param h Hue, between 0 and 1.
inline std::array<uint8_t, 3> oklchToRgb(float l, float c, float h)
{
    using std::pow;

    constexpr Matrix<float, 3> m1inv{{{{1.22701385110352, -0.557799980651822, 0.281256148966468},
                                       {-0.0405801784232806, 1.11225686961683, -0.0716766786656012},
                                       {-0.0763812845057069, -0.421481978418013, 1.58616322044079}}}};
    constexpr Matrix<float, 3> m2inv{{{{1, 0.396337792173768, 0.215803758060759},
                                       {1, -0.105561342323656, -0.0638541747717059},
                                       {1, -0.0894841820949657, -1.29148553786409}}}};
    constexpr Matrix<float, 3> xyzToRgb{
        {{{3.2406, -1.5372, -0.4986}, {-0.9689, 1.8758, 0.0415}, {0.0557, -0.2040, 1.0570}}}};
    auto hr = 2 * std::numbers::pi * h;
    auto a = c * cosf(hr);
    auto b = c * sinf(hr);
    Vector lab{l, a, b};
    auto lms = m2inv * lab;
    for (auto &&x : lms)
        x = pow(x, 3);
    auto xyz = m1inv * lms;
    auto rgb = xyzToRgb * xyz;
    for (auto &&x : rgb)
    {
        if (x <= 0.0031308F)
            x *= 12.92F;
        else
            x = 1.055F * pow(x, 1.0F / 2.4F) - 0.055F;
        x = std::round(255 * std::clamp(x, 0.F, 1.F));
    }
    return {(uint8_t)rgb[0], (uint8_t)rgb[1], (uint8_t)rgb[2]};
}

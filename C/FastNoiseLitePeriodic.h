// MIT License
//
// Copyright(c) 2023 Jordan Peck (jordan.me2@gmail.com)
// Copyright(c) 2023 Contributors
// Copyright(c) 2026 tritao (periodic subset)
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files(the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and / or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions :
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//

/* ProceduralKit periodic subset v1. Based on FastNoiseLite's value hash,
 * Hermite interpolation and axis order. Separate from upstream fnl_state.
 * Coordinates/periods must satisfy the fieldkit contract before calling.
 */
#ifndef FNL_PERIODIC_VALUE_H
#define FNL_PERIODIC_VALUE_H
#include <stdint.h>
#include <math.h>

static inline uint32_t fnlp_wrap(int32_t cell, uint32_t period) {
    uint32_t r = (uint32_t)(cell < 0 ? -cell : cell) % period;
    return cell < 0 && r != 0u ? period - r : r;
}
static inline uint32_t fnlp_hash(uint32_t seed, uint32_t x, uint32_t y, uint32_t z) {
    uint32_t h = seed ^ (x * UINT32_C(501125321))
                      ^ (y * UINT32_C(1136930381))
                      ^ (z * UINT32_C(1720413743));
    h *= UINT32_C(0x27d4eb2d);
    return h;
}
static inline float fnlp_value(uint32_t seed, uint32_t x, uint32_t y, uint32_t z) {
    uint32_t h = fnlp_hash(seed, x, y, z);
    h *= h;
    h ^= h << 19;
    /* Interpret the signed two's-complement value without signed overflow
     * or an implementation-defined uint32 -> int32 conversion. */
    return h < UINT32_C(0x80000000)
        ? (float)h * 0x1p-31f
        : -(float)(~h + UINT32_C(1)) * 0x1p-31f;
}
static inline float fnlp_hermite(float t) { return t * t * (3.0f - 2.0f * t); }
static inline float fnlp_lerp(float a, float b, float t) { return a + t * (b - a); }
static inline float fnlPeriodicValue(uint32_t seed, const float *p,
                                     const uint32_t *period, uint32_t dimensions) {
    uint32_t c0[3] = {0, 0, 0}, c1[3] = {0, 0, 0};
    float t[3] = {0, 0, 0}, v[8];
    for (uint32_t a = 0; a < dimensions; a++) {
        float cell = floorf(p[a]);
        c0[a] = fnlp_wrap((int32_t)cell, period[a]);
        c1[a] = (c0[a] + 1u) % period[a];
        t[a] = fnlp_hermite(p[a] - cell);
    }
    for (uint32_t i = 0; i < (1u << dimensions); i++)
        v[i] = fnlp_value(seed, (i & 1u) ? c1[0] : c0[0],
                               (i & 2u) ? c1[1] : c0[1],
                               (i & 4u) ? c1[2] : c0[2]);
    float xy0 = fnlp_lerp(fnlp_lerp(v[0], v[1], t[0]), fnlp_lerp(v[2], v[3], t[0]), t[1]);
    if (dimensions == 2u) return xy0;
    float xy1 = fnlp_lerp(fnlp_lerp(v[4], v[5], t[0]), fnlp_lerp(v[6], v[7], t[0]), t[1]);
    return fnlp_lerp(xy0, xy1, t[2]);
}
#endif

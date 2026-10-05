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

// ProceduralKit periodic subset v1. Same formulas/order as the C mirror.
// Callers validate finite |coordinates| <= 2^20, dimensions 2/3,
// and active periods 1..65536 before invoking this pure kernel.
#ifndef FNL_PERIODIC_VALUE_GLSL
#define FNL_PERIODIC_VALUE_GLSL
uint fnlp_wrap(int cell, uint period) {
    // GLSL signed remainder with negative operands is not portable.
    uint r = uint(cell < 0 ? -cell : cell) % period;
    return cell < 0 && r != 0u ? period - r : r;
}
float fnlp_value(uint seed, uint x, uint y, uint z) {
    uint h = seed ^ (x * 501125321u) ^ (y * 1136930381u) ^ (z * 1720413743u);
    h *= 0x27d4eb2du;
    h *= h;
    h ^= h << 19;
    return h < 0x80000000u ? float(h) * (1.0 / 2147483648.0)
                         : -float(~h + 1u) * (1.0 / 2147483648.0);
}
float fnlp_hermite(float t) { return t * t * (3.0 - 2.0 * t); }
float fnlp_lerp(float a, float b, float t) { return a + t * (b - a); }
float fnlPeriodicValue(uint seed, vec3 p, uvec3 period, uint dimensions) {
    uvec3 c0 = uvec3(0), c1 = uvec3(0);
    vec3 t = vec3(0);
    float v[8];
    for (uint a = 0u; a < dimensions; a++) {
        float cell = floor(p[a]);
        c0[a] = fnlp_wrap(int(cell), period[a]);
        c1[a] = (c0[a] + 1u) % period[a];
        t[a] = fnlp_hermite(p[a] - cell);
    }
    for (uint i = 0u; i < (1u << dimensions); i++)
        v[i] = fnlp_value(seed, (i & 1u) != 0u ? c1[0] : c0[0],
                               (i & 2u) != 0u ? c1[1] : c0[1],
                               (i & 4u) != 0u ? c1[2] : c0[2]);
    float xy0 = fnlp_lerp(fnlp_lerp(v[0], v[1], t[0]), fnlp_lerp(v[2], v[3], t[0]), t[1]);
    if (dimensions == 2u) return xy0;
    float xy1 = fnlp_lerp(fnlp_lerp(v[4], v[5], t[0]), fnlp_lerp(v[6], v[7], t[0]), t[1]);
    return fnlp_lerp(xy0, xy1, t[2]);
}
#endif

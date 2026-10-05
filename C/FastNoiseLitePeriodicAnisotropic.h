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
#ifndef FASTNOISELITE_PERIODIC_ANISOTROPIC_H
#define FASTNOISELITE_PERIODIC_ANISOTROPIC_H
#include "FastNoiseLitePeriodicPerlin.h"
#include "FastNoiseLitePeriodicWhite.h"
/* Validated integer transform, source periods and transformed coordinates. */
static inline float fnlPeriodicAnisotropic(uint32_t seed,const float *p,const int32_t *matrix,
        const uint32_t *source_period,uint32_t dimensions,uint32_t basis) {
    float q[3]={0,0,0};
    for(uint32_t r=0;r<dimensions;r++) for(uint32_t c=0;c<dimensions;c++)
        q[r]+=(float)matrix[r*3+c]*p[c];
    if(basis==0) return fnlPeriodicValue(seed,q,source_period,dimensions);
    if(basis==1) return fnlPeriodicPerlin(seed,q,source_period,dimensions);
    return fnlPeriodicWhite(seed,q,source_period,dimensions);
}
#endif

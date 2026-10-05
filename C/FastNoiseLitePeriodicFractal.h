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

#ifndef FNL_PERIODIC_FRACTAL_H
#define FNL_PERIODIC_FRACTAL_H
#include "FastNoiseLitePeriodicPerlin.h"
#include "FastNoiseLitePeriodicCellular.h"
#include "FastNoiseLitePeriodicSimplexLattice.h"
/* Validated inputs only. basis: value/perlin/cellular=0/1/2.
 * fold: fbm/ridged/turbulence=0/1/2. channel: F1/F2/F2-F1=0/1/2. */
static inline float fnlPeriodicFractal(uint32_t seed,const float *p,const uint32_t *period,
    uint32_t dimensions,uint32_t basis,uint32_t fold,uint32_t octaves,
    uint32_t lacunarity,float gain,uint32_t metric,uint32_t channel) {
    float q[3]={0,0,0};uint32_t lattice[3]={1,1,1};
    for(uint32_t a=0;a<dimensions;a++){q[a]=p[a];lattice[a]=period[a];}
    float sum=0,weight=0,amplitude=1;
    for(uint32_t i=0;i<octaves;i++) {
        float n;
        if(basis==0u) n=fnlPeriodicValue(seed+i,q,lattice,dimensions);
        else if(basis==1u) n=fnlPeriodicPerlin(seed+i,q,lattice,dimensions);
        else if(basis==3u || basis==4u) n=fnlPeriodicSimplexLattice(seed+i,q,lattice,dimensions);
        else {
            fnlp_cellular_result c=fnlPeriodicCellular(seed+i,q,lattice,dimensions,metric);
            n=channel==0u?c.f1:channel==1u?c.f2:c.f2-c.f1;
        }
        if(fold==1u) n=1.0f-fabsf(n);
        else if(fold==2u) n=fabsf(n);
        sum+=n*amplitude;weight+=amplitude;amplitude*=gain;
        if(i+1u<octaves) for(uint32_t a=0;a<dimensions;a++) {
            q[a]*=(float)lacunarity;lattice[a]*=lacunarity;
        }
    }
    return sum/weight;
}
#endif

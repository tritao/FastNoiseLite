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
/* Validated simplex lattice coordinates; explicit triangular/tetrahedral basis. */
#ifndef FASTNOISELITE_PERIODIC_SIMPLEX_LATTICE_H
#define FASTNOISELITE_PERIODIC_SIMPLEX_LATTICE_H
#include "FastNoiseLitePeriodicPerlin.h"
static inline float fnlp_simplex_corner(uint32_t seed,const int32_t *cell,const float *delta,
        const uint32_t *period,uint32_t d) {
    float radius=.5f;
    for(uint32_t a=0;a<d;a++) radius-=delta[a]*delta[a];
    if(radius<=0) return 0;
    uint32_t identity[3]={0};
    for(uint32_t a=0;a<d;a++) identity[a]=fnlp_wrap(cell[a],period[a]);
    float gradient=fnlp_gradient(seed,identity[0],identity[1],identity[2],delta,d);
    radius*=radius;
    return radius*radius*gradient;
}
static inline float fnlPeriodicSimplexLattice(uint32_t seed,const float *p,const uint32_t *period,uint32_t d) {
    int32_t base[3]={0};float fraction[3]={0};uint32_t order[3]={0,1,2};
    for(uint32_t a=0;a<d;a++) {float b=floorf(p[a]);base[a]=(int32_t)b;fraction[a]=p[a]-b;}
    /* Descending fractional coordinates select the simplex in a lattice cube.
     * Stable axis order owns exact rank ties. */
    for(uint32_t a=1;a<d;a++) for(uint32_t b=a;b>0 && fraction[order[b]]>fraction[order[b-1]];b--) {
        uint32_t t=order[b];order[b]=order[b-1];order[b-1]=t;
    }
    float g=d==2?.21132486540518711775f:1.f/6.f;
    int32_t offsets[3]={0};float value=0;
    for(uint32_t corner=0;corner<=d;corner++) {
        float r[3]={0},sum=0;int32_t cell[3]={0};
        for(uint32_t a=0;a<d;a++) {r[a]=fraction[a]-(float)offsets[a];sum+=r[a];cell[a]=base[a]+offsets[a];}
        for(uint32_t a=0;a<d;a++) r[a]-=g*sum;
        value+=fnlp_simplex_corner(seed,cell,r,period,d);
        if(corner<d) offsets[order[corner]]++;
    }
    /* Conservative absolute bound from independent per-corner maxima. */
    return value*(d==2?32.f:16.f);
}
#endif

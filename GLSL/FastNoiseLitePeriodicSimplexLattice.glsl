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
#ifndef FASTNOISELITE_PERIODIC_SIMPLEX_LATTICE_GLSL
#define FASTNOISELITE_PERIODIC_SIMPLEX_LATTICE_GLSL
float fnlp_simplex_corner(uint seed,ivec3 cell,vec3 delta,uvec3 period,uint d) {
    float radius=.5;
    for(uint a=0u;a<d;a++) radius-=delta[a]*delta[a];
    if(radius<=0.) return 0.;
    uvec3 identity=uvec3(0);
    for(uint a=0u;a<d;a++) identity[a]=fnlp_wrap(cell[a],period[a]);
    float gradient=fnlp_gradient(seed,identity.x,identity.y,identity.z,delta,d);
    radius*=radius;
    return radius*radius*gradient;
}
float fnlPeriodicSimplexLattice(uint seed,vec3 p,uvec3 period,uint d) {
    ivec3 base=ivec3(0);vec3 fraction=vec3(0);uvec3 order=uvec3(0,1,2);
    for(uint a=0u;a<d;a++) {float b=floor(p[a]);base[a]=int(b);fraction[a]=p[a]-b;}
    // Fixed pairwise ranks avoid a lane-dependent insertion-sort loop.
    for(uint a=0u;a<d;a++) {
        uint rank=0u;
        for(uint b=0u;b<d;b++) if(fraction[b]>fraction[a] || (fraction[b]==fraction[a] && b<a)) rank++;
        order[rank]=a;
    }
    float g=d==2u?.21132486540518711775:1./6.;
    ivec3 offsets=ivec3(0);float value=0.;
    for(uint corner=0u;corner<=d;corner++) {
        vec3 r=vec3(0);float sum=0.;ivec3 cell=ivec3(0);
        for(uint a=0u;a<d;a++) {r[a]=fraction[a]-float(offsets[a]);sum+=r[a];cell[a]=base[a]+offsets[a];}
        for(uint a=0u;a<d;a++) r[a]-=g*sum;
        value+=fnlp_simplex_corner(seed,cell,r,period,d);
        if(corner<d) offsets[order[corner]]++;
    }
    return value*(d==2u?32.:16.);
}
#endif

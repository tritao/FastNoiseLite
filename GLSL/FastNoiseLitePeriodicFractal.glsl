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

#ifndef FNL_PERIODIC_FRACTAL_GLSL
#define FNL_PERIODIC_FRACTAL_GLSL
// Include the periodic value, Perlin and cellular mirrors first.
float fnlPeriodicFractal(uint seed,vec3 p,uvec3 period,uint dimensions,
    uint basis,uint fold,uint octaves,uint lacunarity,float gain,uint metric,uint channel) {
    vec3 q=p;uvec3 lattice=period;
    float sum=0.0,weight=0.0,amplitude=1.0;
    for(uint i=0u;i<octaves;i++) {
        float n;
        if(basis==0u) n=fnlPeriodicValue(seed+i,q,lattice,dimensions);
        else if(basis==1u) n=fnlPeriodicPerlin(seed+i,q,lattice,dimensions);
        else if(basis==3u || basis==4u) n=fnlPeriodicSimplexLattice(seed+i,q,lattice,dimensions);
        else {
            fnlp_cellular_result c=fnlPeriodicCellular(seed+i,q,lattice,dimensions,metric);
            n=channel==0u?c.distances[0]:channel==1u?c.distances[1]:c.distances[1]-c.distances[0];
        }
        if(fold==1u) n=1.0-abs(n);
        else if(fold==2u) n=abs(n);
        sum+=n*amplitude;weight+=amplitude;amplitude*=gain;
        if(i+1u<octaves) for(uint a=0u;a<dimensions;a++) {
            q[a]*=float(lacunarity);lattice[a]*=lacunarity;
        }
    }
    return sum/weight;
}
#endif

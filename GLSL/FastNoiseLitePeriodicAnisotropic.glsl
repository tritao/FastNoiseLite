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
#ifndef FASTNOISELITE_PERIODIC_ANISOTROPIC_GLSL
#define FASTNOISELITE_PERIODIC_ANISOTROPIC_GLSL
/* Integer-only mantissa multiplication prevents fused multiply-add from
 * changing transformed coordinates. Coefficients are bounded to 65536. */
float fnlp_integer_product(float p,int m) {
    uint bits=floatBitsToUint(p),mag=uint(m<0?-m:m);
    uint sign=(bits&0x80000000u)^(m<0?0x80000000u:0u);
    uint exponent=(bits>>23)&255u,mantissa=bits&0x7fffffu;
    if(exponent!=0u) mantissa|=0x800000u;
    if(mag==0u || mantissa==0u) return uintBitsToFloat(sign);
    uint low=(mantissa&65535u)*mag;
    uint upper=(mantissa>>16)*mag+(low>>16);
    uint high=upper>>16;low=(upper<<16)|(low&65535u);
    int top=0;
    for(int i=0;i<=40;i++) {
        uint bit=i<32?(low>>uint(i))&1u:(high>>uint(i-32))&1u;
        if(bit!=0u) top=i;
    }
    int power=exponent==0u?-149:int(exponent)-150;
    int out_exp=top+power+127;
    if(out_exp<=0) return uintBitsToFloat(sign|low);
    int shift=top-23;uint rounded;
    if(shift<=0) rounded=low<<uint(-shift);
    else {
        uint n=uint(shift);rounded=(low>>n)|(high<<(32u-n));
        uint remainder=low&((1u<<n)-1u),halfway=1u<<(n-1u);
        if(remainder>halfway || (remainder==halfway && (rounded&1u)!=0u)) rounded++;
        if(rounded==0x1000000u) {rounded>>=1;out_exp++;}
    }
    return uintBitsToFloat(sign|(uint(out_exp)<<23)|(rounded&0x7fffffu));
}
/* Explicit rows avoid GLSL column-major matrix ambiguity. */
float fnlPeriodicAnisotropic(uint seed,vec3 p,ivec3 row0,ivec3 row1,ivec3 row2,
        uvec3 source_period,uint dimensions,uint basis) {
    ivec3 rows[3]=ivec3[3](row0,row1,row2);vec3 q=vec3(0);
    for(uint r=0u;r<dimensions;r++) for(uint c=0u;c<dimensions;c++)
        q[r]+=fnlp_integer_product(p[c],rows[r][c]);
    if(basis==0u) return fnlPeriodicValue(seed,q,source_period,dimensions);
    if(basis==1u) return fnlPeriodicPerlin(seed,q,source_period,dimensions);
    return fnlPeriodicWhite(seed,q,source_period,dimensions);
}
#endif

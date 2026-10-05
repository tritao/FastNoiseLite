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
#ifndef FASTNOISELITE_PERIODIC_CELLULAR_DETAIL_H
#define FASTNOISELITE_PERIODIC_CELLULAR_DETAIL_H
#include "FastNoiseLitePeriodicCellular.h"
static float fnlp_detail_clamp(float x,float lo,float hi) {return x<lo?lo:x>hi?hi:x;}
/* Projection onto each piece of the competing Manhattan site's region.
 * Low saturated coordinates either stay fixed or enter the affine region.
 * Each selected piece is a box intersected with one halfspace. */
static float fnlp_detail_manhattan_pair(const float *a,const float *b,unsigned d) {
    float best=0x1.fffffep127f;
    for(unsigned mask=0;mask<(1u<<d);mask++) {
        float u[3],lo[3],hi[3],knots[7]={0},fixed=0;unsigned n=0,nk=1;
        for(unsigned j=0;j<d;j++) {
            float h=fabsf(b[j]-a[j])*.5f;
            if(h==0) continue;
            float sign=b[j]>=a[j]?1:-1;
            float v=-sign*((a[j]+b[j])*.5f);
            if(v>=h) {fixed+=h;continue;}
            if(v< -h && !(mask&(1u<<j))) {fixed-=h;continue;}
            u[n]=v;lo[n]=v< -h?-h:v;hi[n]=h;
            knots[nk++]=lo[n]-v;knots[nk++]=h-v;n++;
        }
        float target=-fixed,upper=0;
        for(unsigned j=0;j<n;j++) upper+=hi[j];
        if(upper<target) continue;
        for(unsigned j=1;j<nk;j++) for(unsigned k=j;k>0 && knots[k]<knots[k-1];k--) {
            float t=knots[k];knots[k]=knots[k-1];knots[k-1]=t;
        }
        float lambda=0;
        for(unsigned k=0;k<nk;k++) {
            float left=knots[k],sum=0;unsigned slope=0;
            for(unsigned j=0;j<n;j++) {
                float v=u[j]+left;sum+=fnlp_detail_clamp(v,lo[j],hi[j]);
                if(left>=lo[j]-u[j] && left<hi[j]-u[j]) slope++;
            }
            lambda=left;
            if(sum>=target) break;
            if(slope && k+1<nk) {
                float candidate=left+(target-sum)/(float)slope;
                if(candidate<=knots[k+1]) {lambda=candidate;break;}
            }
        }
        float cost=0;
        for(unsigned j=0;j<n;j++) {float v=fnlp_detail_clamp(u[j]+lambda,lo[j],hi[j])-u[j];cost+=v*v;}
        if(cost<best) best=cost;
    }
    return sqrtf(best);
}
/* A Chebyshev competitor region is a union of 2*d cones. Projecting
 * onto one cone reduces to a bounded one-dimensional quadratic. */
static float fnlp_detail_chebyshev_pair(const float *a,const float *b,unsigned d) {
    float best=0x1.fffffep127f;
    for(unsigned j=0;j<d;j++) for(int sign=-1;sign<=1;sign+=2) {
        float separation=(float)sign*(b[j]-a[j]);
        if(separation<0) continue;
        float threshold=separation*.5f,tp=-(float)sign*a[j],dist[2];unsigned n=0;
        for(unsigned i=0;i<d;i++) if(i!=j) dist[n++]=fabsf(b[i]);
        if(n==2 && dist[0]<dist[1]) {float t=dist[0];dist[0]=dist[1];dist[1]=t;}
        float sum=0;
        for(unsigned k=0;k<=n;k++) {
            float lower=fmaxf(threshold,k<n?dist[k]:0),upper=k?dist[k-1]:0x1.fffffep127f;
            if(lower<=upper) {
                float t=fnlp_detail_clamp((tp+sum)/(float)(k+1),lower,upper),cost=(t-tp)*(t-tp);
                for(unsigned i=0;i<n;i++) {float v=fmaxf(dist[i]-t,0);cost+=v*v;}
                if(cost<best) best=cost;
            }
            if(k<n) sum+=dist[k];
        }
    }
    return sqrtf(best);
}


typedef struct fnlp_cellular_detail_result {
    fnlp_cellular_result cell;
    float offset[3],border;
} fnlp_cellular_detail_result;
static inline fnlp_cellular_detail_result fnlPeriodicCellularDetail(uint32_t seed,const float *query,const uint32_t *period,uint32_t dimensions,uint32_t metric) {
    unsigned d=dimensions;
    float jitter=d==2?.43701595f:.39614353f;
    int radius=metric==0?3:metric==1?(d==2?4:5):(d==2?3:4),zr=d==3?radius:0;

        fnlp_cellular_result a=fnlPeriodicCellular(seed,query,period,d,metric);
        fnlp_cellular_detail_result result={0};
        result.cell.f1=a.f1;result.cell.f2=a.f2;result.cell.hash=a.hash;
        float aa=0,fraction[3]={0};int32_t center[3]={0};
        for(unsigned j=0;j<3;j++) {result.cell.cell[j]=a.cell[j];result.offset[j]=-a.delta[j];aa+=a.delta[j]*a.delta[j];}
        for(unsigned j=0;j<d;j++) {float base=floorf(query[j]+.5f);center[j]=(int32_t)base;fraction[j]=query[j]-base;}
        float best=0x1.fffffep127f;
        for(int x=-radius;x<=radius;x++) for(int y=-radius;y<=radius;y++) for(int z=-zr;z<=zr;z++) {
            uint32_t cell[3]={fnlp_wrap(center[0]+x,period[0]),fnlp_wrap(center[1]+y,period[1]),0};
            if(d==3) cell[2]=fnlp_wrap(center[2]+z,period[2]);
            uint32_t hash=fnlp_hash(seed,cell[0],cell[1],cell[2]),idx=hash&(d==2?510u:1020u);
            const float *vectors=d==2?FNLP_RAND_VECS_2D:FNLP_RAND_VECS_3D;
            float b[3]={((float)x-fraction[0])+vectors[idx]*jitter,((float)y-fraction[1])+vectors[idx+1]*jitter,0};
            if(d==3) b[2]=((float)z-fraction[2])+vectors[idx+2]*jitter;
            float bb=0,ab=0;for(unsigned j=0;j<d;j++) {float v=b[j]-a.delta[j];ab+=v*v;bb+=b[j]*b[j];}
            if(ab==0) continue;
            float border;
            if(metric==0) border=fmaxf(0,(bb-aa)/(2*sqrtf(ab)));
            else if(metric==1) border=fnlp_detail_manhattan_pair(a.delta,b,d);
            else border=fnlp_detail_chebyshev_pair(a.delta,b,d);
            if(border<best) best=border;
        }
        result.border=best;return result;
}
#endif

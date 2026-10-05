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
#ifndef FASTNOISELITE_PERIODIC_CELLULAR_DETAIL_GLSL
#define FASTNOISELITE_PERIODIC_CELLULAR_DETAIL_GLSL
float fnlp_detail_clamp(float x,float lo,float hi) {return x<lo?lo:x>hi?hi:x;}
/* Projection onto each piece of the competing Manhattan site's region.
 * Low saturated coordinates either stay fixed_sum or enter the affine region.
 * Each selected piece is a box intersected with one halfspace. */
float fnlp_detail_manhattan_pair(vec3 a,vec3 b,uint d) {
    float best=3.402823466e38;
    for(uint mask=0u;mask<(1u<<d);mask++) {
        float u[3],lo[3],hi[3],knots[7]=float[7](0.,0.,0.,0.,0.,0.,0.),fixed_sum=0.;uint n=0u,nk=1u;
        for(uint j=0u;j<d;j++) {
            float h=abs(b[j]-a[j])*.5f;
            if(h==0.) continue;
            float sign=b[j]>=a[j]?1.:-1.;
            float v=-sign*((a[j]+b[j])*.5f);
            if(v>=h) {fixed_sum+=h;continue;}
            if(v< -h && ((mask&(1u<<j))==0u)) {fixed_sum-=h;continue;}
            u[n]=v;lo[n]=v< -h?-h:v;hi[n]=h;
            knots[nk++]=lo[n]-v;knots[nk++]=h-v;n++;
        }
        float target=-fixed_sum,upper=0.;
        for(uint j=0u;j<n;j++) upper+=hi[j];
        if(upper<target) continue;
        for(uint j=1u;j<nk;j++) for(uint k=j;k>0u && knots[k]<knots[k-1u];k--) {
            float t=knots[k];knots[k]=knots[k-1u];knots[k-1u]=t;
        }
        float lambda=0.;
        for(uint k=0u;k<nk;k++) {
            float left=knots[k],sum=0.;uint slope=0u;
            for(uint j=0u;j<n;j++) {
                float v=u[j]+left;sum+=fnlp_detail_clamp(v,lo[j],hi[j]);
                if(left>=lo[j]-u[j] && left<hi[j]-u[j]) slope++;
            }
            lambda=left;
            if(sum>=target) break;
            if(slope>0u && k+1u<nk) {
                float candidate=left+(target-sum)/float(slope);
                if(candidate<=knots[k+1u]) {lambda=candidate;break;}
            }
        }
        float cost=0.;
        for(uint j=0u;j<n;j++) {float v=fnlp_detail_clamp(u[j]+lambda,lo[j],hi[j])-u[j];cost+=v*v;}
        if(cost<best) best=cost;
    }
    return sqrt(best);
}
/* A Chebyshev competitor region is a union of 2*d cones. Projecting
 * onto one cone reduces to a bounded one-dimensional quadratic. */
float fnlp_detail_chebyshev_pair(vec3 a,vec3 b,uint d) {
    float best=3.402823466e38;
    for(uint j=0u;j<d;j++) for(int sign=-1;sign<=1;sign+=2) {
        float separation=float(sign)*(b[j]-a[j]);
        if(separation<0.) continue;
        float threshold=separation*.5f,tp=-float(sign)*a[j],dist[2];uint n=0u;
        for(uint i=0u;i<d;i++) if(i!=j) dist[n++]=abs(b[i]);
        if(n==2u && dist[0]<dist[1]) {float t=dist[0];dist[0]=dist[1];dist[1]=t;}
        float sum=0.;
        for(uint k=0u;k<=n;k++) {
            float lower=max(threshold,k<n?dist[k]:0.),upper=k>0u?dist[k-1u]:3.402823466e38;
            if(lower<=upper) {
                float t=fnlp_detail_clamp((tp+sum)/float(k+1u),lower,upper),cost=(t-tp)*(t-tp);
                for(uint i=0u;i<n;i++) {float v=max(dist[i]-t,0.);cost+=v*v;}
                if(cost<best) best=cost;
            }
            if(k<n) sum+=dist[k];
        }
    }
    return sqrt(best);
}


struct fnlp_cellular_detail_result {fnlp_cellular_result cell;vec3 offset;float border;};
fnlp_cellular_detail_result fnlPeriodicCellularDetail(uint seed,vec3 query,uvec3 period,uint dimensions,uint metric) {
    uint d=dimensions;
    float jitter=d==2u?.43701595:.39614353;
    int radius=metric==0u?3:metric==1u?(d==2u?4:5):(d==2u?3:4),zr=d==3u?radius:0;
    fnlp_cellular_result a=fnlPeriodicCellular(seed,query,period,d,metric);
    fnlp_cellular_detail_result result;result.cell=a;result.offset=-a.delta;
    float aa=0.;vec3 fraction=vec3(0);ivec3 center=ivec3(0);
    for(uint j=0u;j<3u;j++) aa+=a.delta[j]*a.delta[j];
    for(uint j=0u;j<d;j++) {float base=floor(query[j]+.5);center[j]=int(base);fraction[j]=query[j]-base;}
    float best=3.402823466e38;
    for(int x=-radius;x<=radius;x++) for(int y=-radius;y<=radius;y++) for(int z=-zr;z<=zr;z++) {
        uvec3 cell=uvec3(fnlp_wrap(center[0]+x,period[0]),fnlp_wrap(center[1]+y,period[1]),0);
        if(d==3u) cell[2]=fnlp_wrap(center[2]+z,period[2]);
        uint hash=fnlp_hash(seed,cell[0],cell[1],cell[2]),idx=hash&(d==2u?510u:1020u);
        float vx=d==2u?FNLP_RAND_VECS_2D[idx]:FNLP_RAND_VECS_3D[idx];
        float vy=d==2u?FNLP_RAND_VECS_2D[idx+1u]:FNLP_RAND_VECS_3D[idx+1u];
        vec3 b=vec3((float(x)-fraction[0])+vx*jitter,(float(y)-fraction[1])+vy*jitter,0);
        if(d==3u) b[2]=(float(z)-fraction[2])+FNLP_RAND_VECS_3D[idx+2u]*jitter;
        float bb=0.,ab=0.;for(uint j=0u;j<d;j++) {float v=b[j]-a.delta[j];ab+=v*v;bb+=b[j]*b[j];}
        if(ab==0.) continue;
        float border;
        if(metric==0u) border=max(0.,(bb-aa)/(2.*sqrt(ab)));
        else if(metric==1u) border=fnlp_detail_manhattan_pair(a.delta,b,d);
        else border=fnlp_detail_chebyshev_pair(a.delta,b,d);
        if(border<best) best=border;
    }
    result.border=best;return result;
}
#endif

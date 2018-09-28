/* Electronic Strings by NR4^QM/Team210 - 64k Intro at Deadline 2018
 * Copyright (C) 2018 Alexander Kraus <nr4@z10.info>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#version 130
uniform float iTime;
uniform vec2 iResolution;
uniform sampler2D iFont;
uniform float iFontWidth;

float iScale = 0.; //TODO: uniform this

const vec3 c = vec3(1.,0.,-1.), 
    tc = vec3(234.,82.,35.)/255.;
const float pi = acos(-1.);

const float t_logo = 30.;
const float t_tunnel = 7000.;

float t = 0.;

// Extract specific byte from font texture 
float fdata(float byte)
{
    float byte_n = mod(byte, 4.), 
        rgba_n = byte - byte_n,
        x = mod(rgba_n, iFontWidth), 
        y = (rgba_n-x)/iFontWidth;
    return byte_n;
}

float dletter(vec2 x, int which, float size)
{
    // first byte is number of contained glyphs
    
    
//     +texture2D(iFont, x)
    return 0.;
}
float rand(vec2 a0)
{
    return fract(sin(dot(a0.xy, vec2(12.9898,78.233)))*43758.5453);
}

//distance to quadratic bezier spline with parameter t
float dist(vec2 p0,vec2 p1,vec2 p2,vec2 x,float t)
{
    t = clamp(t, 0., 1.);
    return length(x-pow(1.-t,2.)*p0-2.*(1.-t)*t*p1-t*t*p2);
}

//minimum distance to quadratic bezier spline
float dsp(vec2 p0, vec2 p1, vec2 p2, vec2 x)
{
    //coefficients for 0 = t^3 + a * t^2 + b * t + c
    vec2 E = x-p0, F = p2-2.*p1+p0, G = p1-p0;
    vec3 ai = vec3(3.*dot(G,F), 2.*dot(G,G)-dot(E,F), -dot(E,G))/dot(F,F);

    //discriminant and helpers
    float tau = ai.x/3., p = ai.y-tau*ai.x, q = - tau*(tau*tau+p)+ai.z, dis = q*q/4.+p*p*p/27.;
    
    //triple real root
    if(dis > 0.) 
    {
        vec2 ki = -.5*q*c.xx+sqrt(dis)*c.xz, ui = sign(ki)*pow(abs(ki), c.xx/3.);
        return dist(p0,p1,p2,x,ui.x+ui.y-tau);
    }
    
    //three distinct real roots
    float fac = sqrt(-4./3.*p), arg = acos(-.5*q*sqrt(-27./p/p/p))/3.;
    vec3 t = c.zxz*fac*cos(arg*c.xxx+c*pi/3.)-tau;
    return min(
        dist(p0,p1,p2,x, t.x),
        min(
            dist(p0,p1,p2,x,t.y),
            dist(p0,p1,p2,x,t.z)
        )
    );
}

//minimum distance to linear bezier spline
float dsg(vec2 p0, vec2 p1, vec2 x)
{
    vec2 d = p1-p0;
    float t = clamp(dot(x-p0,d)/dot(d,d),0.,1.);
    return length(x-mix(p0,p1,t));
}

mat3 rot(vec3 theta)
{
    vec3 c = cos(theta), s = sin(theta);
    return mat3(c.x,0.,s.x,0.,1.,0.,-s.x,0.,c.x)
        *mat3(1., 0., 0., 0., c.y, -s.y, 0., s.y, c.y)
        *mat3(c.z,-s.z,0., s.z,c.z,0., 0.,0.,1.);
}

float zextrude(float z, float d2d, float h)
{
    vec2 d = abs(vec2(min(d2d, 0.),z))-h*c.yx;
    return min(max(d.x,d.y),0.)+length(max(d,0.));
}

float circle(vec2 x, float r, float  w)
{
    return w-abs(length(x)-r);
}

float circleseg(vec2 x, float r, float w, float p1, float p2)
{
    float p = clamp(sign(x.y)*acos(x.x/length(x)), p1, p2);
    vec2 y = r*vec2(cos(p), sin(p));
    return circle(x-y,0.,w);
}

float line(vec2 x, vec2 p1, vec2 p2, float w)
{
    vec2 d = p2-p1;
    return length(x-mix(p1,p2,clamp(dot(x-p1,d)/dot(d,d),0.,1.)))-w;
}

float line(vec3 x, vec3 p1, vec3 p2, float w)
{
    vec3 d = p2-p1;
    return length(x-mix(p1,p2,clamp(dot(x-p1,d)/dot(d,d),0.,1.)))-w;
}

float pict210(vec3 x, float h, float r, float w)
{
    float ret = min(
        zextrude(x.z,circle(x.xy-r*c.xy,r,w),h),
        zextrude(x.z, -line(x.xy,-r*c.yx,r*c.yx,w), h)
    );
    return min(ret,
        zextrude(x.z,circleseg(x.xy+r*c.xy,r,w,-pi/2.,pi/2.),h)
    );
}

vec2 logo(vec3 x)
{
    //coord transforms
    vec3 x0 = x;
    x = rot(c.yxy*pi/4.)*x;
    
    //parameters
    float h = .05;//+.15*iScale;

    //big 210 logo
    vec2 sdb, sdf = vec2(pict210(x-.11*c.yyx,h,.2,.05),2.);

    // add grid with 210s in the background on a plane
    vec3 xa = x;
    x += c.yxy*t*.1+10.*c.xxy;
    sdb = vec2(x.z+0.21,2.);
    sdf = mix(sdb, sdf, step(sdf.x,sdb.x));
    vec3 y = vec3(mod(x.xy, .1)-.05, x.z),
        ind = x-y;
    float r = rand(ind.xy);
    h = 4.*r*h*abs(8.*sin(.05*2.*pi*(ind.x-.12)));
    if(r > .5) y = rot(.5*pi*c.yyx)*y;
    sdb = vec2(pict210(y+0.21*c.yyx,.5*h,.022,.005),r);
    sdf = mix(sdb, sdf, step(sdf.x, sdb.x));
    
    float guard = -length(max(abs(y)-.05*c.xxy-(.1+.01*abs(ind.x))*c.yyx,0.));
    guard = abs(guard)+.1*.1;
    sdf.x = min(sdf.x, guard);

    x = xa;

    //bullshit flying through the eye
    x -= .1*(t-12.)*c.yyx;    

    float dw = .04, b = step(-.15,x0.z)*.12;//*(.12+.08*sin(2.*pi*x0.z));

    y = mod(x,dw)-.5*dw;
    ind = x-y;
    float cone = -b+length(ind.xy-.2*c.xy);

    if(cone < 2.*dw)
    {
        for(float i=-1.; i<=1.; i+=1.)
        {
            for(float j=-1.; j<=1.; j+=1.)
            {
                for(float k=-1.;k<=1.; k+=1.)
                {
                    vec3 indd = ind + vec3(i,j,k);
                    float ka = -b+length(indd.xy-.2*c.xy);
                    if(ka > .5*dw) continue;
                    r = rand(10.*indd.xy+21.*indd.yz+35.*indd.zx);
                    vec3 z = rot(1.e1*(t+12.)*vec3(.1,.2,.3)*r)*(y+dw*vec3(i,j,k));
                    
                    float w = .2*dw+dw*.15*r*(.5+.5*sin(1.e1*t+1.33e1*(r+length(x.xy))));
                    
                    sdb = vec2(length(max(abs(z)-w*c.xxx,0.)), r);
                    sdf = mix(sdb, sdf, step(sdf.x,sdb.x));
                }
            }
        }
        float guard = -length(max(abs(y)-.5*dw*c.xxx,0.));
        guard = abs(guard)+dw*.1;
        sdf.x = min(sdf.x, guard);
    }
    else 
    {
        sdb = vec2(cone, 1.);
        sdf = mix(sdf, sdb, step(sdb.x, sdf.x));
    }
    return sdf;
}

vec2 tunnel(vec3 x)
{
    //build the tunnel with floor
    float R = .5;

    vec2 sdf = vec2(zextrude(x.y,-abs(R-length(x.xz))+.01, 2000.), 3.),
        sda = vec2(x.z+.8*R, 3.);
    sdf = mix(sda,sdf,step(sdf.x,sda.x));

    //add profile to the floor; TODO: better shit.
    float dr = .04;
    vec3 y = vec3(mod(x.xy, dr)-.5*dr, x.z);
    sda = vec2(zextrude(y.z+.4, .4*dr-length(y.xy), .005), 3.);
    sdf = mix(sda, sdf, step(sdf.x, sda.x));

    //add pipes to the wall
    float pipe_radius = .02,
        phi = 0.,
        rout = R-4.*pipe_radius;
    vec3 ka = vec3(cos(phi), 0., sin(phi));
    sda = vec2(line(x, R*ka, rout*ka, pipe_radius), 3.);
    sdf = mix(sda, sdf, step(sdf.x, sda.x));
    sda = vec2(line(x, rout*ka, rout*ka+.5*c.yxy, pipe_radius), 3.);
    sdf = mix(sda, sdf, step(sdf.x, sda.x));

    //vec3 px = vec3(R, acos(max(x.x,x.z)/R), x.y);
    //float dr = circle(

    return sdf;
}

vec2 scene(vec3 x)
{
    //TODO: add time-dependent scenes here
    if(t < t_logo) return logo(x);
    else if(t < t_tunnel) return tunnel(x);
}

const float dx = 1.e-4;
vec3 normal(vec3 x)
{
    float s = scene(x).x;
    return normalize(vec3(scene(x+dx*c.xyy).x-s, scene(x+dx*c.yxy).x-s, scene(x+dx*c.yyx).x-s));
}

vec3 background(vec2 uv)
{
    return c.yyy;
}

vec4 add2(vec4 sdf, vec4 sda)
{
    return vec4(
        min(sdf.x, sda.x), 
        mix(sda.gba, sdf.gba, smoothstep(-1.5/iResolution.y, 1.5/iResolution.y, sda.x))
    );
}

float blend(float tin, float tout)
{
    return smoothstep(tin-.5, tin+.5, t)*(1.-smoothstep(tout-.5, tout+.5, t));
}

vec4 textbox(vec2 x, vec2 p0, vec2 p1, float w)
{
    float dw = .004;

    vec4 sdf = c.xyyy, sda; 

    // Text background
    sda = vec4(line(x, p1, p1+w*c.xy, 10.*dw), c.xyy);
    sdf = add2(sdf, sda);
    sdf = add2(sdf, vec4(abs(sda.x-2.*dw) - dw, c.xxy));

    // Connection line and dots
    // dot at p1
    sda = vec4(length(x-p1)-2.*dw, tc);
    sdf = add2(sdf, sda);
    sdf = add2(sdf, vec4(abs(sda.x-2.*dw) - dw, c.xxy));

    // dot at p0
    sda = vec4(length(x-p0)-2.*dw, tc);
    sdf = add2(sdf, sda);
    sdf = add2(sdf, vec4(abs(sda.x-2.*dw) - dw, c.xxy));
    
    //line
    sdf = add2(sdf, vec4(line(x,p0,p1,dw), tc));

    return sdf;
}

vec4 textlayer(vec2 x)
{
    float alpha = 0.;
    vec4 sdt =c.xyyy, sda;
    float d = 10.;
    vec2 uv = x;
    {
const vec2 lin[26] = vec2[26](vec2(6.63e-01,3.75e-01),vec2(6.63e-01,4.28e-01),vec2(6.50e-01,4.28e-01),vec2(6.77e-01,4.28e-01),vec2(7.04e-01,3.75e-01),vec2(7.17e-01,3.75e-01),vec2(6.91e-01,3.88e-01),vec2(7.17e-01,3.88e-01),vec2(7.58e-01,3.75e-01),vec2(7.58e-01,4.02e-01),vec2(7.72e-01,3.75e-01),vec2(7.72e-01,4.02e-01),vec2(7.99e-01,3.75e-01),vec2(7.99e-01,3.88e-01),vec2(8.25e-01,3.75e-01),vec2(8.25e-01,3.88e-01),vec2(8.39e-01,3.75e-01),vec2(8.66e-01,3.75e-01),vec2(8.80e-01,4.15e-01),vec2(8.93e-01,4.28e-01),vec2(8.93e-01,4.28e-01),vec2(8.93e-01,3.75e-01),vec2(9.07e-01,3.88e-01),vec2(9.07e-01,4.15e-01),vec2(9.34e-01,3.88e-01),vec2(9.34e-01,4.15e-01)),
quad[48] = vec2[48](vec2(7.17e-01,3.88e-01),vec2(7.17e-01,4.02e-01),vec2(7.04e-01,4.02e-01),vec2(7.04e-01,4.02e-01),vec2(6.91e-01,4.02e-01),vec2(6.91e-01,3.88e-01),vec2(6.91e-01,3.88e-01),vec2(6.91e-01,3.75e-01),vec2(7.04e-01,3.75e-01),vec2(7.45e-01,4.02e-01),vec2(7.31e-01,4.02e-01),vec2(7.31e-01,3.88e-01),vec2(7.31e-01,3.88e-01),vec2(7.31e-01,3.75e-01),vec2(7.45e-01,3.75e-01),vec2(7.45e-01,3.75e-01),vec2(7.58e-01,3.75e-01),vec2(7.58e-01,3.88e-01),vec2(7.45e-01,4.02e-01),vec2(7.58e-01,4.02e-01),vec2(7.58e-01,3.88e-01),vec2(7.72e-01,3.88e-01),vec2(7.72e-01,4.02e-01),vec2(7.85e-01,4.02e-01),vec2(7.85e-01,4.02e-01),vec2(7.99e-01,4.02e-01),vec2(7.99e-01,3.88e-01),vec2(7.99e-01,3.88e-01),vec2(7.99e-01,4.02e-01),vec2(8.12e-01,4.02e-01),vec2(8.12e-01,4.02e-01),vec2(8.25e-01,4.02e-01),vec2(8.25e-01,3.88e-01),vec2(8.39e-01,4.28e-01),vec2(8.79e-01,4.42e-01),vec2(8.39e-01,3.75e-01),vec2(9.07e-01,3.88e-01),vec2(9.07e-01,3.75e-01),vec2(9.21e-01,3.75e-01),vec2(9.21e-01,3.75e-01),vec2(9.34e-01,3.75e-01),vec2(9.34e-01,3.88e-01),vec2(9.34e-01,4.15e-01),vec2(9.34e-01,4.28e-01),vec2(9.21e-01,4.28e-01),vec2(9.21e-01,4.28e-01),vec2(9.07e-01,4.28e-01),vec2(9.07e-01,4.15e-01));
for(int i=0; i<13;++i) d=min(d,dsg(lin[2*i], lin[2*i+1], uv));
for(int i=0; i<16; ++i) d=min(d,dsp(quad[3*i], quad[3*i+1], quad[3*i+2], uv));
    }
    sda = textbox(x, vec2(.44,.155), vec2(.6,.4), .5);
    vec4 sdb = vec4(d-.004, tc);
    sda = add2(sda, sdb); 
    sda = add2(sda, vec4(abs(sdb.x)-.002, c.xxy));
    sdt = add2(sdt, sda);
    alpha = max(alpha, step(sda.x, 1.5/iResolution.y)*blend(11., 25.));
    
    d = 10.;
    {
const vec2 lin[48] = vec2[48](vec2(6.50e-01,1.88e-01),vec2(6.50e-01,2.15e-01),vec2(6.77e-01,1.88e-01),vec2(6.77e-01,2.15e-01),vec2(6.91e-01,1.75e-01),vec2(6.91e-01,2.28e-01),vec2(7.44e-01,1.75e-01),vec2(7.44e-01,2.28e-01),vec2(6.91e-01,2.28e-01),vec2(7.17e-01,2.02e-01),vec2(7.17e-01,2.02e-01),vec2(7.44e-01,2.28e-01),vec2(7.58e-01,1.88e-01),vec2(7.58e-01,1.88e-01),vec2(7.58e-01,2.15e-01),vec2(7.58e-01,2.15e-01),vec2(7.99e-01,1.75e-01),vec2(7.85e-01,1.75e-01),vec2(7.99e-01,2.28e-01),vec2(7.85e-01,2.28e-01),vec2(7.72e-01,1.88e-01),vec2(7.72e-01,2.15e-01),vec2(8.67e-01,2.02e-01),vec2(8.80e-01,2.02e-01),vec2(8.67e-01,1.75e-01),vec2(8.80e-01,1.75e-01),vec2(8.80e-01,1.75e-01),vec2(8.80e-01,2.28e-01),vec2(9.07e-01,1.75e-01),vec2(9.21e-01,1.75e-01),vec2(8.94e-01,1.88e-01),vec2(9.21e-01,1.88e-01),vec2(9.35e-01,1.75e-01),vec2(9.61e-01,2.28e-01),vec2(9.75e-01,1.75e-01),vec2(9.89e-01,1.75e-01),vec2(1.00e+00,2.28e-01),vec2(9.89e-01,2.28e-01),vec2(1.02e+00,1.75e-01),vec2(1.02e+00,2.28e-01),vec2(1.02e+00,2.28e-01),vec2(1.04e+00,2.28e-01),vec2(1.02e+00,2.02e-01),vec2(1.04e+00,2.02e-01),vec2(1.06e+00,1.75e-01),vec2(1.08e+00,2.28e-01),vec2(1.06e+00,2.28e-01),vec2(1.08e+00,1.75e-01)),
quad[60] = vec2[60](vec2(6.50e-01,1.88e-01),vec2(6.50e-01,1.75e-01),vec2(6.63e-01,1.75e-01),vec2(6.63e-01,1.75e-01),vec2(6.77e-01,1.75e-01),vec2(6.77e-01,1.88e-01),vec2(6.77e-01,2.15e-01),vec2(6.77e-01,2.28e-01),vec2(6.63e-01,2.28e-01),vec2(6.63e-01,2.28e-01),vec2(6.50e-01,2.28e-01),vec2(6.50e-01,2.15e-01),vec2(6.63e-01,1.88e-01),vec2(6.63e-01,1.75e-01),vec2(6.77e-01,1.75e-01),vec2(7.85e-01,1.75e-01),vec2(7.72e-01,1.75e-01),vec2(7.72e-01,1.88e-01),vec2(7.72e-01,2.15e-01),vec2(7.72e-01,2.28e-01),vec2(7.85e-01,2.28e-01),vec2(8.26e-01,1.75e-01),vec2(8.13e-01,1.75e-01),vec2(8.13e-01,1.88e-01),vec2(8.13e-01,1.88e-01),vec2(8.13e-01,2.02e-01),vec2(8.26e-01,2.02e-01),vec2(8.26e-01,2.02e-01),vec2(8.39e-01,2.02e-01),vec2(8.39e-01,1.88e-01),vec2(8.39e-01,1.88e-01),vec2(8.39e-01,1.75e-01),vec2(8.26e-01,1.75e-01),vec2(8.67e-01,1.75e-01),vec2(8.53e-01,1.75e-01),vec2(8.53e-01,1.88e-01),vec2(8.53e-01,1.88e-01),vec2(8.53e-01,2.02e-01),vec2(8.67e-01,2.02e-01),vec2(9.21e-01,1.88e-01),vec2(9.21e-01,2.02e-01),vec2(9.07e-01,2.02e-01),vec2(9.07e-01,2.02e-01),vec2(8.94e-01,2.02e-01),vec2(8.94e-01,1.88e-01),vec2(8.94e-01,1.88e-01),vec2(8.94e-01,1.75e-01),vec2(9.07e-01,1.75e-01),vec2(9.89e-01,1.75e-01),vec2(1.00e+00,1.75e-01),vec2(1.00e+00,1.88e-01),vec2(1.00e+00,1.88e-01),vec2(1.00e+00,2.02e-01),vec2(9.89e-01,2.02e-01),vec2(9.89e-01,2.02e-01),vec2(9.75e-01,2.02e-01),vec2(9.75e-01,2.15e-01),vec2(9.75e-01,2.15e-01),vec2(9.75e-01,2.28e-01),vec2(9.89e-01,2.28e-01));
for(int i=0; i<24;++i) d=min(d,dsg(lin[2*i], lin[2*i+1], uv));
for(int i=0; i<20; ++i) d=min(d,dsp(quad[3*i], quad[3*i+1], quad[3*i+2], uv));
    }
    sda = textbox(x, vec2(.485, .1), vec2(.6, .2), .5);
    sdb = vec4(d-.004, tc);
    sda = add2(sda, sdb); 
    sda = add2(sda, vec4(abs(sdb.x)-.002, c.xxy));
    sdt = add2(sdt, sda);
    alpha = max(alpha, step(sda.x, 1.5/iResolution.y)*blend(12., 25.));
    
    d = 10.;
    {
const vec2 lin[60] = vec2[60](vec2(6.50e-01,-2.50e-02),vec2(6.50e-01,2.83e-02),vec2(6.50e-01,2.83e-02),vec2(6.77e-01,-2.50e-02),vec2(6.77e-01,-2.50e-02),vec2(6.77e-01,2.83e-02),vec2(6.91e-01,-2.50e-02),vec2(6.91e-01,2.83e-02),vec2(6.91e-01,2.83e-02),vec2(7.04e-01,2.83e-02),vec2(6.91e-01,1.67e-03),vec2(7.04e-01,1.67e-03),vec2(7.17e-01,-1.17e-02),vec2(7.17e-01,-2.50e-02),vec2(7.58e-01,2.83e-02),vec2(7.31e-01,-1.17e-02),vec2(7.31e-01,-1.17e-02),vec2(7.58e-01,-1.17e-02),vec2(7.58e-01,1.67e-03),vec2(7.58e-01,-2.50e-02),vec2(7.72e-01,-1.17e-02),vec2(7.72e-01,-1.17e-02),vec2(7.72e-01,1.50e-02),vec2(7.72e-01,1.50e-02),vec2(8.13e-01,-2.50e-02),vec2(7.99e-01,-2.50e-02),vec2(8.13e-01,2.83e-02),vec2(7.99e-01,2.83e-02),vec2(7.86e-01,-1.17e-02),vec2(7.86e-01,1.50e-02),vec2(8.81e-01,1.67e-03),vec2(8.94e-01,1.67e-03),vec2(8.81e-01,-2.50e-02),vec2(8.94e-01,-2.50e-02),vec2(8.94e-01,-2.50e-02),vec2(8.94e-01,2.83e-02),vec2(9.21e-01,-2.50e-02),vec2(9.35e-01,-2.50e-02),vec2(9.08e-01,-1.17e-02),vec2(9.35e-01,-1.17e-02),vec2(9.49e-01,-2.50e-02),vec2(9.75e-01,2.83e-02),vec2(1.00e+00,1.67e-03),vec2(1.02e+00,1.67e-03),vec2(1.02e+00,1.67e-03),vec2(1.02e+00,-1.17e-02),vec2(9.89e-01,-1.17e-02),vec2(9.89e-01,1.50e-02),vec2(1.00e+00,2.83e-02),vec2(1.02e+00,2.83e-02),vec2(1.03e+00,-2.50e-02),vec2(1.03e+00,2.83e-02),vec2(1.03e+00,2.83e-02),vec2(1.06e+00,2.83e-02),vec2(1.03e+00,1.67e-03),vec2(1.06e+00,1.67e-03),vec2(1.07e+00,-2.50e-02),vec2(1.10e+00,2.83e-02),vec2(1.07e+00,2.83e-02),vec2(1.10e+00,-2.50e-02)),
quad[51] = vec2[51](vec2(7.04e-01,2.83e-02),vec2(7.17e-01,2.83e-02),vec2(7.17e-01,1.50e-02),vec2(7.17e-01,1.50e-02),vec2(7.17e-01,1.67e-03),vec2(7.04e-01,1.67e-03),vec2(7.04e-01,1.67e-03),vec2(7.17e-01,1.67e-03),vec2(7.17e-01,-1.17e-02),vec2(7.99e-01,-2.50e-02),vec2(7.86e-01,-2.50e-02),vec2(7.86e-01,-1.17e-02),vec2(7.86e-01,1.50e-02),vec2(7.86e-01,2.83e-02),vec2(7.99e-01,2.83e-02),vec2(8.40e-01,-2.50e-02),vec2(8.27e-01,-2.50e-02),vec2(8.27e-01,-1.17e-02),vec2(8.27e-01,-1.17e-02),vec2(8.27e-01,1.67e-03),vec2(8.40e-01,1.67e-03),vec2(8.40e-01,1.67e-03),vec2(8.53e-01,1.67e-03),vec2(8.53e-01,-1.17e-02),vec2(8.53e-01,-1.17e-02),vec2(8.53e-01,-2.50e-02),vec2(8.40e-01,-2.50e-02),vec2(8.81e-01,-2.50e-02),vec2(8.67e-01,-2.50e-02),vec2(8.67e-01,-1.17e-02),vec2(8.67e-01,-1.17e-02),vec2(8.67e-01,1.67e-03),vec2(8.81e-01,1.67e-03),vec2(9.35e-01,-1.17e-02),vec2(9.35e-01,1.67e-03),vec2(9.21e-01,1.67e-03),vec2(9.21e-01,1.67e-03),vec2(9.08e-01,1.67e-03),vec2(9.08e-01,-1.17e-02),vec2(9.08e-01,-1.17e-02),vec2(9.08e-01,-2.50e-02),vec2(9.21e-01,-2.50e-02),vec2(9.89e-01,1.50e-02),vec2(9.89e-01,2.83e-02),vec2(1.00e+00,2.83e-02),vec2(1.00e+00,-2.50e-02),vec2(9.89e-01,-2.50e-02),vec2(9.89e-01,-1.17e-02),vec2(1.00e+00,-2.50e-02),vec2(1.02e+00,-2.50e-02),vec2(1.02e+00,-1.17e-02));
for(int i=0; i<30;++i) d=min(d,dsg(lin[2*i], lin[2*i+1], uv));
for(int i=0; i<17; ++i) d=min(d,dsp(quad[3*i], quad[3*i+1], quad[3*i+2], uv));
    }
    sda = textbox(x, vec2(.485, .03), vec2(.6,.0), .5);
    sdb = vec4(d-.004, tc);
    sda = add2(sda, sdb); 
    sda = add2(sda, vec4(abs(sdb.x)-.002, c.xxy));
    sdt = add2(sdt, sda);
    alpha = max(alpha, step(sda.x, 1.5/iResolution.y)*blend(13., 25.));
    
    d = 10.;
    {
const vec2 lin[50] = vec2[50](vec2(6.50e-01,-2.25e-01),vec2(6.50e-01,-1.72e-01),vec2(7.03e-01,-2.25e-01),vec2(7.03e-01,-1.72e-01),vec2(6.50e-01,-1.72e-01),vec2(6.77e-01,-1.98e-01),vec2(6.77e-01,-1.98e-01),vec2(7.03e-01,-1.72e-01),vec2(7.31e-01,-2.25e-01),vec2(7.31e-01,-1.72e-01),vec2(7.17e-01,-2.25e-01),vec2(7.44e-01,-2.25e-01),vec2(7.17e-01,-1.72e-01),vec2(7.44e-01,-1.72e-01),vec2(7.85e-01,-2.25e-01),vec2(7.71e-01,-2.25e-01),vec2(7.85e-01,-1.72e-01),vec2(7.71e-01,-1.72e-01),vec2(7.58e-01,-2.12e-01),vec2(7.58e-01,-1.85e-01),vec2(7.99e-01,-2.12e-01),vec2(7.99e-01,-2.12e-01),vec2(7.99e-01,-1.85e-01),vec2(7.99e-01,-1.85e-01),vec2(8.13e-01,-2.25e-01),vec2(8.26e-01,-2.25e-01),vec2(8.39e-01,-1.72e-01),vec2(8.26e-01,-1.72e-01),vec2(8.53e-01,-1.98e-01),vec2(8.53e-01,-2.12e-01),vec2(8.80e-01,-1.98e-01),vec2(8.80e-01,-2.25e-01),vec2(8.94e-01,-2.52e-01),vec2(8.94e-01,-1.98e-01),vec2(8.94e-01,-1.98e-01),vec2(9.07e-01,-1.98e-01),vec2(8.94e-01,-2.25e-01),vec2(9.07e-01,-2.25e-01),vec2(9.35e-01,-2.52e-01),vec2(9.35e-01,-1.98e-01),vec2(9.35e-01,-1.98e-01),vec2(9.48e-01,-1.98e-01),vec2(9.35e-01,-2.25e-01),vec2(9.48e-01,-2.25e-01),vec2(1.02e+00,-2.25e-01),vec2(1.02e+00,-1.98e-01),vec2(1.06e+00,-2.12e-01),vec2(1.06e+00,-1.72e-01),vec2(1.04e+00,-1.98e-01),vec2(1.07e+00,-1.98e-01)),
quad[54] = vec2[54](vec2(7.71e-01,-2.25e-01),vec2(7.58e-01,-2.25e-01),vec2(7.58e-01,-2.12e-01),vec2(7.58e-01,-1.85e-01),vec2(7.58e-01,-1.72e-01),vec2(7.71e-01,-1.72e-01),vec2(8.26e-01,-2.25e-01),vec2(8.39e-01,-2.25e-01),vec2(8.39e-01,-2.12e-01),vec2(8.39e-01,-2.12e-01),vec2(8.39e-01,-1.98e-01),vec2(8.26e-01,-1.98e-01),vec2(8.26e-01,-1.98e-01),vec2(8.13e-01,-1.98e-01),vec2(8.13e-01,-1.85e-01),vec2(8.13e-01,-1.85e-01),vec2(8.13e-01,-1.72e-01),vec2(8.26e-01,-1.72e-01),vec2(8.53e-01,-2.12e-01),vec2(8.53e-01,-2.25e-01),vec2(8.67e-01,-2.25e-01),vec2(8.67e-01,-2.25e-01),vec2(8.80e-01,-2.25e-01),vec2(8.80e-01,-2.12e-01),vec2(9.07e-01,-2.25e-01),vec2(9.21e-01,-2.25e-01),vec2(9.21e-01,-2.12e-01),vec2(9.21e-01,-2.12e-01),vec2(9.21e-01,-1.98e-01),vec2(9.07e-01,-1.98e-01),vec2(9.48e-01,-2.25e-01),vec2(9.61e-01,-2.25e-01),vec2(9.61e-01,-2.12e-01),vec2(9.61e-01,-2.12e-01),vec2(9.61e-01,-1.98e-01),vec2(9.48e-01,-1.98e-01),vec2(9.89e-01,-2.25e-01),vec2(9.75e-01,-2.25e-01),vec2(9.75e-01,-2.12e-01),vec2(9.75e-01,-2.12e-01),vec2(9.75e-01,-1.98e-01),vec2(9.89e-01,-1.98e-01),vec2(9.89e-01,-1.98e-01),vec2(1.00e+00,-1.98e-01),vec2(1.00e+00,-2.12e-01),vec2(1.00e+00,-2.12e-01),vec2(1.00e+00,-2.25e-01),vec2(9.89e-01,-2.25e-01),vec2(1.02e+00,-2.12e-01),vec2(1.02e+00,-1.98e-01),vec2(1.03e+00,-1.98e-01),vec2(1.06e+00,-2.12e-01),vec2(1.06e+00,-2.25e-01),vec2(1.07e+00,-2.25e-01));
for(int i=0; i<25;++i) d=min(d,dsg(lin[2*i], lin[2*i+1], uv));
for(int i=0; i<18; ++i) d=min(d,dsp(quad[3*i], quad[3*i+1], quad[3*i+2], uv));
    }
    sda = textbox(x, vec2(.435, -.02), vec2(.6,-.2), .5);
    sdb = vec4(d-.004, tc);
    sda = add2(sda, sdb); 
    sda = add2(sda, vec4(abs(sdb.x)-.002, c.xxy));
    sdt = add2(sdt, sda);
    alpha = max(alpha, step(sda.x, 1.5/iResolution.y)*blend(14., 25.));
    
    d = 10.;
    {
const vec2 lin[48] = vec2[48](vec2(6.50e-01,-3.98e-01),vec2(6.50e-01,-4.12e-01),vec2(6.77e-01,-3.98e-01),vec2(6.77e-01,-4.12e-01),vec2(7.03e-01,-3.98e-01),vec2(7.03e-01,-4.12e-01),vec2(7.17e-01,-3.98e-01),vec2(7.17e-01,-4.12e-01),vec2(7.44e-01,-3.98e-01),vec2(7.44e-01,-4.12e-01),vec2(7.71e-01,-3.98e-01),vec2(7.71e-01,-4.12e-01),vec2(7.85e-01,-3.98e-01),vec2(7.85e-01,-4.12e-01),vec2(8.11e-01,-3.98e-01),vec2(8.11e-01,-4.12e-01),vec2(8.38e-01,-3.98e-01),vec2(8.38e-01,-4.12e-01),vec2(8.52e-01,-4.25e-01),vec2(8.52e-01,-4.25e-01),vec2(8.93e-01,-3.98e-01),vec2(8.66e-01,-3.98e-01),vec2(8.93e-01,-3.98e-01),vec2(8.66e-01,-4.25e-01),vec2(8.66e-01,-4.25e-01),vec2(8.93e-01,-4.25e-01),vec2(9.07e-01,-3.85e-01),vec2(9.20e-01,-3.72e-01),vec2(9.20e-01,-3.72e-01),vec2(9.20e-01,-4.25e-01),vec2(9.34e-01,-4.12e-01),vec2(9.34e-01,-3.85e-01),vec2(9.61e-01,-4.12e-01),vec2(9.61e-01,-3.85e-01),vec2(9.75e-01,-4.25e-01),vec2(9.75e-01,-4.25e-01),vec2(9.89e-01,-4.12e-01),vec2(9.89e-01,-3.98e-01),vec2(9.89e-01,-3.85e-01),vec2(9.89e-01,-3.85e-01),vec2(1.02e+00,-4.25e-01),vec2(1.02e+00,-3.98e-01),vec2(1.04e+00,-4.12e-01),vec2(1.04e+00,-4.25e-01),vec2(1.07e+00,-4.12e-01),vec2(1.07e+00,-3.85e-01),vec2(1.06e+00,-3.98e-01),vec2(1.08e+00,-3.98e-01)),
quad[75] = vec2[75](vec2(6.50e-01,-4.12e-01),vec2(6.50e-01,-4.25e-01),vec2(6.63e-01,-4.25e-01),vec2(6.63e-01,-4.25e-01),vec2(6.77e-01,-4.25e-01),vec2(6.77e-01,-4.12e-01),vec2(6.77e-01,-4.12e-01),vec2(6.77e-01,-4.25e-01),vec2(6.90e-01,-4.25e-01),vec2(6.90e-01,-4.25e-01),vec2(7.03e-01,-4.25e-01),vec2(7.03e-01,-4.12e-01),vec2(7.17e-01,-4.12e-01),vec2(7.17e-01,-4.25e-01),vec2(7.31e-01,-4.25e-01),vec2(7.31e-01,-4.25e-01),vec2(7.44e-01,-4.25e-01),vec2(7.44e-01,-4.12e-01),vec2(7.44e-01,-4.12e-01),vec2(7.44e-01,-4.25e-01),vec2(7.57e-01,-4.25e-01),vec2(7.57e-01,-4.25e-01),vec2(7.71e-01,-4.25e-01),vec2(7.71e-01,-4.12e-01),vec2(7.85e-01,-4.12e-01),vec2(7.85e-01,-4.25e-01),vec2(7.98e-01,-4.25e-01),vec2(7.98e-01,-4.25e-01),vec2(8.11e-01,-4.25e-01),vec2(8.11e-01,-4.12e-01),vec2(8.11e-01,-4.12e-01),vec2(8.11e-01,-4.25e-01),vec2(8.25e-01,-4.25e-01),vec2(8.25e-01,-4.25e-01),vec2(8.38e-01,-4.25e-01),vec2(8.38e-01,-4.12e-01),vec2(9.34e-01,-4.12e-01),vec2(9.34e-01,-4.25e-01),vec2(9.47e-01,-4.25e-01),vec2(9.47e-01,-4.25e-01),vec2(9.61e-01,-4.25e-01),vec2(9.61e-01,-4.12e-01),vec2(9.61e-01,-3.85e-01),vec2(9.61e-01,-3.72e-01),vec2(9.47e-01,-3.72e-01),vec2(9.47e-01,-3.72e-01),vec2(9.34e-01,-3.72e-01),vec2(9.34e-01,-3.85e-01),vec2(9.89e-01,-4.12e-01),vec2(9.89e-01,-4.25e-01),vec2(1.00e+00,-4.25e-01),vec2(1.03e+00,-3.98e-01),vec2(1.04e+00,-3.98e-01),vec2(1.04e+00,-4.12e-01),vec2(1.02e+00,-4.12e-01),vec2(1.02e+00,-3.98e-01),vec2(1.03e+00,-3.98e-01),vec2(1.06e+00,-4.25e-01),vec2(1.07e+00,-4.25e-01),vec2(1.07e+00,-4.12e-01),vec2(1.07e+00,-3.85e-01),vec2(1.07e+00,-3.72e-01),vec2(1.08e+00,-3.72e-01),vec2(1.11e+00,-4.25e-01),vec2(1.10e+00,-4.25e-01),vec2(1.10e+00,-4.12e-01),vec2(1.10e+00,-4.12e-01),vec2(1.10e+00,-3.98e-01),vec2(1.11e+00,-3.98e-01),vec2(1.11e+00,-3.98e-01),vec2(1.12e+00,-3.98e-01),vec2(1.12e+00,-4.12e-01),vec2(1.12e+00,-4.12e-01),vec2(1.12e+00,-4.25e-01),vec2(1.11e+00,-4.25e-01));
for(int i=0; i<24;++i) d=min(d,dsg(lin[2*i], lin[2*i+1], uv));
for(int i=0; i<25; ++i) d=min(d,dsp(quad[3*i], quad[3*i+1], quad[3*i+2], uv));
    }
    sda = textbox(x, vec2(.365, -.045), vec2(.6,-.4), .5);
    sdb = vec4(d-.004, tc);
    sda = add2(sda, sdb); 
    sda = add2(sda, vec4(abs(sdb.x)-.002, c.xxy));
    sdt = add2(sdt, sda);
    alpha = max(alpha, step(sda.x, .01)*blend(15., 25.));

    return vec4(sdt.gba, alpha)*smoothstep(1.5/iResolution.y, -1.5/iResolution.y, sdt.x);
}

vec3 raymarch(vec2 uv, float time)
{
    t = time;
    
    int nmax = 200;
    float d=0., dmax=50., prec=5.e-4;
    vec2 s; 
    vec3 ro = c.yzy, ta = c.yyy, r = c.xyy, u = c.yyx, rt = ta + uv.x*r + uv.y*u, rd = normalize(rt-ro), x, col;

    //TODO: add scene dependent raymarching setup here 
    if(t<t_logo) //TODO: remove if it does not affect precision
    {
        dmax = 20.;
        nmax = 100;
    }

    int quit = 0;
    for(int i=0; i<nmax; ++i)
    {
        x = ro + d * rd;
        s = scene(x);
        if(s.x < prec) break;
        if((d > dmax) || (i == nmax-1)) 
        {
            quit = 1;
            break;
        }
        d += s.x;
    }

    if(quit == 1)
        col = background(uv);
    else
    {
        vec3 n = normal(x); 
        col=c.yyy;

        //TODO: add scene dependent lighting setup and fog here
        if(t<t_tunnel)
        {
            if(s.y <= 1.)
            {
                vec3 l = c.yyx, re = normalize(reflect(-l,n)), v = normalize(x-ro),
                    amb = abs(rot(5.e-2*vec3(10.,21.,35.)*s.y+t+x.z)*c.xxy),
                    dif = abs(rot(5.e-2*vec3(21.,35.,10.)*s.y+t+x.z)*c.yxx),
                    spc = abs(rot(5.e-2*vec3(35.,10.,21.)*s.y+t+x.z)*c.xyx);
                col = .1*amb + dot(l,n)*.3*dif + pow(abs(dot(re,v)),2.)*.9*spc;

                col = mix(col, .01*c.xxx, tanh(2.8*x.y));
            }
            else if(s.y <= 2.)
            {
                vec3 l = c.yyx, re = normalize(reflect(-l,n)), v = normalize(x-ro),
                    amb = .1*c.xxx*length(abs(rot(5.e-2*vec3(10.,21.,35.)*s.y+t+x.z)*c.xxy)),
                    dif = .1*c.xxx*length(abs(rot(5.e-2*vec3(21.,35.,10.)*s.y+t+x.z)*c.yxx)),
                    spc = .2*c.xxx*length(abs(rot(5.e-2*vec3(35.,10.,21.)*s.y+t+x.z)*c.xyx));
                col = .1*amb + dot(l,n)*.3*dif + pow(abs(dot(re,v)),2.)*.9*spc;

                col = mix(col, .01*c.xxx, tanh(2.8*x.y));
            }
            else if(s.y <= 3.)
            {
                vec3 l = normalize(c.yxx), re = normalize(reflect(-l,n)), v = normalize(x-ro),
                    amb = abs(rot(5.e-2*vec3(10.,21.,35.)*s.y+t+x.z)*c.xxy),
                    dif = abs(rot(5.e-2*vec3(21.,35.,10.)*s.y+t+x.z)*c.yxx),
                    spc = abs(rot(5.e-2*vec3(35.,10.,21.)*s.y+t+x.z)*c.xyx);
                col = .1*amb + dot(l,n)*.2*dif + pow(abs(dot(re,v)),2.)*.9*spc;

                col = mix(col, .01*c.xxx, tanh(.8*x.y));
            }
        }
    }

    //TODO: add text here
    vec4 text = textlayer(uv);
    col = mix(col, text.rgb, 1.*text.a);
    
    return col;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    iScale = 2.*mod(iTime, .5);
    vec2 uv = fragCoord/iResolution.yy-.5;
    const float dt = 1.e-2;

    vec3 col = c.yyy;
    // 4x FSAA
    if(false)
    {
        float dx = .75/iResolution.y;
        for(float i = -dx; i < dx; i += dx)
            for(float j = -dx; j < dx; j += dx)
                 col += raymarch(uv+vec2(i,j), iTime);
        col *= .25;
    }
    // 2x TAA
    else if(false)
    {
        float dt = 1.e-2;
        col = .5*raymarch(uv, iTime)+.5*raymarch(uv, iTime+dt);
    }
    else col = raymarch(uv, iTime);
    
    col += texture(iFont, uv+.5).rgb;
    
    fragColor = vec4(col,1.0);
}




void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}

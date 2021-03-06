 float iScale;
 float iNBeats;
 
const float pi = acos(-1.);
const vec3 c = vec3(1.,0.,-1.);

// Hash function
float rand(vec2 x)
{
    return fract(sin(dot(x-1. ,vec2(12.9898,78.233)))*43758.5453);
}

float rand(vec3 x)
{
    return fract(sin(dot(x-1. ,vec3(12.9898,78.233,33.1818)))*43758.5453);
}

vec3 rand3(vec3 x)
{
    return vec3(rand(x.x*c.xx),rand(x.y*c.xx),rand(x.z*c.xx));
}

/* compute voronoi distance and closest point.
 * x: coordinate
 * return value: vec3(distance, coordinate of control point)
 */
vec3 vor(vec2 x)
{
    vec2 y = floor(x);
   	float ret = 1.;
    
    //find closest control point. ("In which cell am I?")
    vec2 pf=c.yy, p;
    float df=10., d;
    
    for(int i=-1; i<=1; i+=1)
        for(int j=-1; j<=1; j+=1)
        {
            p = y + vec2(float(i), float(j));
            p += rand(p);
            
            d = length(x-p);
            
            if(d < df)
            {
                df = d;
                pf = p;
            }
        }
    
    //compute voronoi distance: minimum distance to any edge
    for(int i=-1; i<=1; i+=1)
        for(int j=-1; j<=1; j+=1)
        {
            p = y + vec2(float(i), float(j));
            p += rand(p);
            
            vec2 o = p - pf;
            d = length(.5*o-dot(x-pf, o)/dot(o,o)*o);
            ret = min(ret, d);
        }
    
    return vec3(ret, pf);
}

float valuenoise(vec2 x)
{
    vec2 y = floor(x);
    x = fract(x);
    float r00 = -1.+2.*rand(y),
        r10 = -1.+2.*rand(y+c.xy),
        r01 = -1.+2.*rand(y+c.yx),
        r11 = -1.+2.*rand(y+c.xx);
    return mix(
        mix(r00, r10, x.x),
        mix(r01, r11, x.x),
        x.y
    );
}

float valuenoise(vec3 x)
{
    vec3 y = floor(x);
    x = fract(x);
    float r000 = -1.+2.*rand(y),
        r100 = -1.+2.*rand(y+c.xyy),
        r010 = -1.+2.*rand(y+c.yxy),
        r001 = -1.+2.*rand(y+c.yyx),
        r110 = -1.+2.*rand(y+c.xxy),
        r011 = -1.+2.*rand(y+c.yxx),
        r101 = -1.+2.*rand(y+c.xyx),
        r111 = -1.+2.*rand(y+c.xxx);
    return 	mix(
        		mix(
            		mix(r000, r100, x.x),
                    mix(r010, r110, x.x),
                    x.y
                ),
        		mix(
                    mix(r001, r101, x.x),
                    mix(r011, r111, x.x),
                    x.y
                ),
        		x.z);
        
}

vec2 add(vec2 sda, vec2 sdb)
{
    return mix(sda, sdb, step(sdb.x, sda.x));
}

vec2 sub(vec2 sda, vec2 sdb)
{
    return mix(-sda, sdb, step(sda.x, sdb.x));
}

// Distance to line segment
float linesegment(vec3 x, vec3 p0, vec3 p1)
{
    vec3 d = p1-p0;
    float t = clamp(dot(x-p0,d)/dot(d,d),0.,1.);
    return length(x-mix(p0,p1,t));
}

// Stroke
float stroke(float sdf, float w)
{
    return abs(sdf)-w;
}

float zextrude(float z, float d2d, float h)
{
    vec2 d = abs(vec2(min(d2d, 0.),z))-h*c.yx;
    return min(max(d.x,d.y),0.)+length(max(d,0.));
}

/* compute voronoi distance and closest point.
 * x: coordinate
 * return value: vec3(distance, coordinate of control point)
 */
vec4 vor(vec3 x)
{
    vec3 y = floor(x);
   	float ret = 10.;
    
    //find closest control point. ("In which cell am I?")
    vec3 pf=c.yyy, p;
    float df=100., d;
    
    for(int i=-1; i<=1; i+=1)
        for(int j=-1; j<=1; j+=1)
            for(int k=-1; k<=1; k+=1)
            {
                p = y + vec3(float(i), float(j), float(k));
                p += rand3(p);

                d = length(x-p);
				
                if(d < df)
                {
                    df = d;
                    pf = p;
                }
            }
    
    //compute voronoi distance: minimum distance to any edge
    for(int i=-1; i<=1; i+=1)
        for(int j=-1; j<=1; j+=1)
            for(int k=-1; k<=1; k+=1)
            {
                p = y + vec3(float(i), float(j), float(k));
                p += rand3(p);

                vec3 o = p - pf;
                d = abs(.5-dot(x-pf, o)/length(o));
                ret = min(ret, d);
            }
    return vec4(ret, pf);
}

mat3 rot(vec3 p)
{
    return mat3(c.xyyy, cos(p.x), sin(p.x), 0., -sin(p.x), cos(p.x))
        *mat3(cos(p.y), 0., -sin(p.y), c.yxy, sin(p.y), 0., cos(p.y))
        *mat3(cos(p.z), -sin(p.z), 0., sin(p.z), cos(p.z), c.yyyx);
}

//BUILD A WORLD
//THEN BREAK IT
vec3 ind;
vec2 scene(vec3 x)
{
    //x = rot(.05*vec3(1.,2.,3.)*iTime+iNBeats)*x;
    //x += c.yxy*(-30.+.1*iTime);
    vec2 dis = 8.*vec2((.1+.05*iScale)*valuenoise(x-2.-1.*iTime),(.1+.15*iScale)*valuenoise(x.xy-5.-1.*iTime));
    
    
    vec3 v = vor(4.*x.xy-iNBeats-dis);
    float d = stroke(zextrude(x.z+.1-.2*x.y, stroke(v.x,.15), .1+.5*rand(v.yz+7.+iNBeats)+.1*valuenoise(1.*v.yz-.4*iTime)),.05);
    //d = max(-d2, d);
    //artificial guards for artifacts
    float dr = .04;
    vec3 y = mod(x, dr)-.5*dr;
    float guard = -length(max(abs(y)-vec3(.5*dr*c.xx, .6),0.));
    guard = abs(guard)+dr*.1;
    d = min(d, guard);
    /*
    vec3 z = (x-c.yxx);
    z = vec3(mod(z.xy, .2)-.1, z.z);
    d = min(d, length(z)-.1);
    */
    //floor
    //d = min(d, x.z-.1); 
    
    //add spikes
    //vec4 v2 = abs(x.z-valuenoise(4.*x.xy)-.5)*c.xxxx, 
    //vec4 w  = vor(4.*x-3.*rand(iNBeats*c.xx)-(.2+.1*iScale)*.1*valuenoise(x.xy-2.-1.*iTime));
    ind = vec3(v.yz, 0.);
    //d = min(stroke(.4*w.x,5.e-3+1.e-3*iScale), d);
    //d = max(-stroke(v2.x,1.e-1), d);
    //d = min(length(x)-1., d);
    //d = max(-stroke(.4*w.x,/*5.e-3+*/5.e-3*iScale), d);
    return add(vec2(d, 1.), vec2(x.z, 2.));
}

vec3 stdcolor(vec2 x)
{
	return 0.5 + 0.5*cos(iTime+x.xyx+vec3(0,2,4));
}

// Distance to circle
float circle(vec2 x, float r)
{
    return length(x)-r;
}

// Distance to circle segment
float circlesegment(vec2 x, float r, float p0, float p1)
{
    float p = atan(x.y, x.x);
    p = clamp(p, p0, p1);
    return length(x-r*vec2(cos(p), sin(p)));
}

// Distance to line segment
float linesegment(vec2 x, vec2 p0, vec2 p1)
{
    vec2 d = p1-p0;
    float t = clamp(dot(x-p0,d)/dot(d,d),0.,1.);
    return length(x-mix(p0,p1,t));
}

// Distance to 210 logo
float logo(vec2 x, float r)
{
    return min(
        min(circle(x+r*c.zy, r), linesegment(x,r*c.yz, r*c.yx)),
        circlesegment(x+r*c.xy, r, -.5*pi, .5*pi)
    );
}

//performs raymarching
//scene: name of the scene function
//xc: 	 name of the coordinate variable
//ro:	 name of the ray origin variable
//d:	 name of the distance variable
//dir:	 name of the direction variable
//s:	 name of the scenestruct variable
//N:	 number of iterations used
//eps:	 exit criterion
//flag:  name of the flag to set if raymarching succeeded
#define raymarch(scene, xc, ro, d, dir, s, N, eps, flag) \
	flag = false;\
	for(int i=0; i<N; ++i)\
    {\
        xc = ro + d*dir;\
        s = scene(xc);\
        if(s.x < eps)\
        {\
            flag = true;\
            break;\
        }\
        d += s.x;\
    }

//computes normal with finite differences
//scene: name of the scene function
//n:	 name of the normal variable
//eps:	 precision of the computation
//xc:	 location of normal evaluation
#define calcnormal(scene, n, eps, xc) \
	{\
        float ss = scene(xc).x;\
        n = normalize(vec3(scene(xc+eps*c.xyy).xc-ss,\
                           scene(xc+eps*c.yxy).xc-ss,\
                           scene(xc+eps*c.yyx).xc-ss));\
    }

//camera setup
//camera: camera function with camera(out vec3 ro, out vec3 r, out vec3 u, out vec3 t)
//ro:	  name of the ray origin variable
//r:	  name of the right variable
//u:	  name of the up variable
//t:	  name of the target variable
//uv:	  fragment coordinate
//dir:	  name of the dir variable
#define camerasetup(camera, ro, r, u, t, uv, dir) \
	{\
        camera(ro, r, u, t);\
        t += uv.x*r+uv.y*u;\
        dir = normalize(t-ro);\
    }

//camera for scene 1
void camera1(out vec3 ro, out vec3 r, out vec3 u, out vec3 t)
{
    ro = c.yyx;
    r = c.xyy;
    u = c.yxx;
    t = c.yxy;
}

vec3 synthcol(float scale, float phase)
{
    vec3 c2 = vec3(207.,30.,102.)/255.,
        c3 = vec3(245., 194., 87.)/255.;
    mat3 r1 = rot((5.e-1*phase)*vec3(1.1,1.3,1.5));
    return 
        (
            1.1*mix
            (
                -(cross(c2, r1*c2)),
                -(r1*c2), 
                scale
            )
        );
}

vec3 color(float rev, float ln, float index, vec2 uv, vec3 x)
{
    vec3 col = c.yyy;
    if(index == 1.)
    {
   		vec3 c1 = stdcolor(x.xy+.5*rand(ind.xy+17.)+iNBeats), 
        	c2 = stdcolor(x.xy+x.yz+x.zx+.5*rand(ind.xy+12.)+iNBeats+11.+uv), 
            c3 = stdcolor(x.xy+x.yz+x.zx+.5*rand(ind.xy+15.)+iNBeats+23.+uv);
		col = .1*c1*vec3(1.,1.,1.) + .2*c1*vec3(1.,1.,1.)*ln + vec3(1.,1.,1.)*pow(rev,2.*(2.-1.5*clamp(iScale,0.,1.))) + 2.*c1*pow(rev, 8.)+3.*c1*pow(rev, 16.);
		col = clamp(.33*col, 0., 1.);
        //col = abs(col);
	}
    else if(index == 2.)
    {
        return .1*c.xxx;
    }
    return col;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 uv = fragCoord/iResolution.yy-.5;
	
    iScale = mod(iTime, .5)/.5;
    iNBeats = floor(iTime);
    
    vec3 col = c.yyy;
    
    if(iTime < 1000.) //scene 1
    {
        //use raymarching in this scene
    	vec3 ro, r, u, t, x, dir;
    	camerasetup(camera1, ro, r, u, t, uv, dir);
    	
        float d = -(ro.z-.7)/dir.z;
    
        bool hit;
        vec2 s;
        raymarch(scene, x, ro, d, dir, s, 300, 1.e-4, hit);
        if(hit == false)
        {
            fragColor = c.yyyx;
            return;
        }
    
        vec3 n;
        calcnormal(scene, n, 1.e-3, x);
    
        vec3 l = x+2.*c.yyx, re = normalize(reflect(-l,n)), v = normalize(x-ro);
        float rev = abs(dot(re,v)), ln = abs(dot(l,n));
        
        col = color(rev, ln, s.y, uv, x);
        
        //reflections
        dir = normalize(reflect(dir, n));
        //dir = normalize(refract(dir, n, .9));
        d = 5.e-2;
        ro = x;
        raymarch(scene, x, ro, d, dir, s, 150, 5.e-4, hit);
        
        if(hit == false)
        {
            fragColor = c.yyyx;
            //return;
        }
        else
        {
            calcnormal(scene, n, 1.e-3, x);
            l = x+2.*c.yyx;
            re = normalize(reflect(-l,n)); 
            v = normalize(x-ro);
            rev = abs(dot(re,v));
            ln = abs(dot(l,n));

            col = mix(col, color(rev, ln, s.y, uv, x), .7);

            //second reflection
            dir = normalize(reflect(dir, n));
	        d = 5.e-2;
            ro = x;
            raymarch(scene, x, ro, d, dir, s, 150, 5.e-4, hit);

            if(hit == false)
            {
                fragColor = c.yyyx;
                //return;
            }
            else
            {
                calcnormal(scene, n, 1.e-3, x);
                l = x+2.*c.yyx;
                re = normalize(reflect(-l,n)); 
                v = normalize(x-ro);
                rev = abs(dot(re,v));
                ln = abs(dot(l,n));

                col = mix(col, color(rev, ln, s.y, uv, x), .5);
                //third reflection
                
                dir = normalize(reflect(dir, n));
                d = 5.e-2;
                ro = x;
                raymarch(scene, x, ro, d, dir, s, 150, 5.e-4, hit);

                if(hit == false)
                {
                    fragColor = c.yyyx;
                    //return;
                }
                else
                {
                    calcnormal(scene, n, 1.e-3, x);
                    l = x+2.*c.yyx;
                    re = normalize(reflect(-l,n)); 
                    v = normalize(x-ro);
                    rev = abs(dot(re,v));
                    ln = abs(dot(l,n));

                    col = mix(col, color(rev, ln, s.y, uv, x), .3);
                }
				
            }
        }
        
        
        
        //Portability
        //col = clamp((tanh(.7*col)), 0., 1.);
        //fog
        col = mix(col, c.yyy, tanh(2.e-1*(abs(x.y+x.z))));
        
    }
    
    //210 logo
    col = mix(clamp(col,c.yyy,c.xxx), c.xxx, smoothstep(1.5/iResolution.y, -1.5/iResolution.y, stroke(logo(uv-vec2(-.45,.45),.02),.005)));
    //trendy display lines
    col += vec3(0., 0.05, 0.1)*sin(uv.y*1050.+ 5.*iTime);
    
    fragColor = vec4(col,1.0);
}

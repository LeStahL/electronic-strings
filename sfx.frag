#version 130

uniform float iBlockOffset;
uniform float iSampleRate;
uniform float iVolume;
#define PI radians(180.)
float _sin(float a) { return sin(2. * PI * mod(a,1.)); }
float _unisin(float a,float b) { return (.5*_sin(a) + .5*_sin((1.+b)*a)); }
float _sq(float a) { return sign(2.*fract(a) - 1.); }
float _squ(float a,float pwm) { return sign(2.*fract(a) - 1. + pwm); }
float _tri(float a) { return (4.*abs(fract(a)-.5) - 1.); }
float _saw(float a) { return (2.*fract(a) - 1.); }
float quant(float a,float div,float invdiv) { return floor(div*a+.5)*invdiv; }
float quanti(float a,float div) { return floor(div*a+.5)/div; }
float clip(float a) { return clamp(a,-1.,1.); }
float theta(float x) { return smoothstep(0., 0.01, x); }
float freqC1(float note){ return 32.7 * pow(2.,note/12.); }
float minus1hochN(int n) { return (1. - 2.*float(n % 2)); }

const float BPM = 80.; 
const float BPS = BPM/60.;
const float SPB = 60./BPM;

float doubleslope(float t, float a, float d, float s)
{
    return smoothstep(-.00001,a,t) - (1.-s) * smoothstep(0.,d,t-a);
}

float s_atan(float a) { return 2./PI * atan(a); }
float s_crzy(float amp) { return clamp( s_atan(amp) - 0.1*cos(0.9*amp*exp(amp)), -1., 1.); }
float squarey(float a, float edge) { return abs(a) < edge ? a : floor(4.*a+.5)*.25; } 

float TRISQ(float t, float f, int MAXN, float MIX, float INR, float NDECAY, float RES, float RES_Q)
{
    float ret = 0.;
    
    int Ninc = 8; // try this: leaving out harmonics...
    
    for(int N=0; N<=MAXN; N+=Ninc)
    {
        float mode     = 2.*float(N) + 1.;
        float inv_mode = 1./mode; 		// avoid division? save table of Nmax <= 20 in some array or whatever
        float comp_TRI = (N % 2 == 1 ? -1. : 1.) * inv_mode*inv_mode;
        float comp_SQU = inv_mode;
        float filter_N = pow(1. + pow(float(N) * INR,2.*NDECAY),-.5) + RES * exp(-pow(float(N)*INR*RES_Q,2.));

        ret += (MIX * comp_TRI + (1.-MIX) * comp_SQU) * filter_N * _sin(mode * f * t);
    }
    
    return ret;
}

float QTRISQ(float t, float f, float QUANT, int MAXN, float MIX, float INR, float NDECAY, float RES, float RES_Q)
{
    return TRISQ(quant(t,QUANT,1./QUANT), f, MAXN, MIX, INR, NDECAY, RES, RES_Q);
}

float env_ADSR(float x, float L, float A, float D, float S, float R)
{
    float att = pow(x/A,8.);
    float dec = S + (1.-S) * exp(-(x-A)/D);
    float rel = (x < L-R) ? 1. : pow((L-x)/R,4.);

    return (x < A ? att : dec) * rel;
    
}

float macesaw(float t, float f, float CO, float Q, float det1, float det2, float res, float resQ)
{
    float s = 0.;
    float inv_CO = 1./CO;
    float inv_resQ = 1./resQ;
    float p = f*t;
        for(int N=1; N<=200; N++)
        {
            // saw
            float sawcomp = 2./PI * (1. - 2.*float(N % 2)) * 1./float(N);
            float filterN  = pow(1. + pow(float(N)*f*inv_CO,Q),-.5)
                     + res * exp(-pow((float(N)*f-CO)*inv_resQ,2.));
            
            if(abs(filterN*sawcomp) < 1e-6) break;
        		
            if(det1 > 0. || det2 > 0.)
            {
                s += 0.33 * (_sin(float(N)*p) + _sin(float(N)*p*(1.+det1)) + _sin(float(N)*p*(1.+det2)));
            }
            else
            {
                s += filterN * sawcomp * _sin(float(N)*p);
            }
        }
    return s;
}

float snare(float t, float t_on, float vel)
{
    // #define _tri(a) (4.*abs(fract(a)-.5) - 1.)
    t = t - min(t, t_on);
    float f1 = 6000.;
    float f2 = 800.;
    float f3 = 350.;
    float dec12 = 0.01;
    float dec23 = 0.01;
    float rel = 0.1;
    float snr = _tri(t * (f3 + (f1-f2)*smoothstep(-dec12,0.,-t)
                             + (f2-f3)*smoothstep(-dec12-dec23,-dec12,-t))) * smoothstep(-rel,-dec12-dec23,-t);
        
    //noise part
    float noise = 2. * fract(sin(t * 90.) * 45000.) * doubleslope(t,0.05,0.3,0.3);
   
    return vel * clamp(1.7 * (2. * snr + noise), -1.5, 1.5) * doubleslope(t,0.0,0.25,0.3);
}

float hut(float t, float t_on, float vel)
{
    t = t - min(t, t_on);
    float noise = fract(sin(t * 90.) * 45000.);
    noise = 1./(1.+noise);
    return vel * 2. * noise * doubleslope(t,0.,0.12,0.0);
    
    // might think of this one! - maybe tune length / pitch
    //float kick_blubb = (1.-exp(-1000.*t))*exp(-30.*t) * _sin((400.-200.*t)*t * _saw(4.*f*t));
}

float shake(float t, float t_on, float vel) // shaker is just some mod of hihat (hut)
{
    t = t - min(t, t_on);
    return vel * 0.5 * fract(sin(t * 90.) * 45000.) * doubleslope(t,0.03,0.15,0.15);
}

float hoskins_noise(float t) // thanks to https://www.shadertoy.com/view/4sjSW1 !
{
    float p = floor(t * (1500.0 * exp(-t*.100)));
	vec2 p2 = fract(vec2(p * 5.3983, p * 5.4427));
    p2 += dot(p2.yx, p2.xy + vec2(21.5351, 14.3137));
	return fract(p2.x * p2.y * 3.4337) * .5 * smoothstep(-.3,0.,-t);    
}

float facekick(float t, float t_on, float vel)
{
    t = t - min(t, t_on); // reset time to Bon event
    
    float f   = 50. + 150. * smoothstep(-0.12, 0., -t);
    float env = smoothstep(0.,0.015,t) * smoothstep(-0.08, 0., 0.16 - t);
    
    float kick_body = env * TRISQ(t, f, 3, 1., 0.8, 8., 4., 1.); // more heavy bass drum: increase reso parameters?
    
    float kick_click = 0.4 * step(t,0.03) * _sin(t*1100. * _saw(t*800.));
    
    float kick_blobb = (1.-exp(-1000.*t))*exp(-40.*t) * _sin((400.-200.*t)*t * _sin(1.*f*t));
    
	return vel * (kick_body + kick_blobb + 0.1*kick_click);
}

float hardkick(float t, float t_on, float vel)
{
    t = t - min(t, t_on); // reset time to Bon event
    
    float f   = 60. + 150. * smoothstep(-0.3, 0., -t);
    float env = smoothstep(0.,0.01,t) * smoothstep(-0.1, 0.2, 0.3 - t);
    
    float kick_body = env * .1*TRISQ(t, f, 100, 1., 1., .1, 16., 10.); // more heavy bass drum: increase reso parameters?
   
    kick_body += .7 * (smoothstep(0.,0.01,t) * smoothstep(-0.2, 0.2, 0.3 - t)) * _sin(t*f*.5);

    float kick_click = 1.5*step(t,0.05) * _sin(t*5000. * _saw(t*1000.));
    
    kick_click = s_atan(40.*(1.-exp(-1000.*t))*exp(-80.*t) * _sin((1200.-1000.*sin(1000.*t*sin(30.*t)))*t));
    
    float kick_blobb = s_crzy(10.*(1.-exp(-1000.*t))*exp(-30.*t) * _sin((300.-300.*t)*t));
    
	return vel * 2.*clamp(kick_body + kick_blobb + kick_click,-1.5,1.5);
}

float oddnoise(float t, float t_on, float vel)
{
    t = t - min(t, t_on); // reset time to Bon event
    
    float kick_click = 1.5*step(t,0.05) * _sin(t*5000. * _saw(t*1000.));
    
    kick_click = s_atan(8.*(1.-exp(-1000.*t))*exp(-80.*t) * _sin((1200.-1000.*sin(1000.*t*t))*t));
    
	return vel * 2.*clamp(kick_click,-1.5,1.5);
}

float shape(float x, float a)
{
    return mix(x,sin(a*PI*x),0.2);
}

float distsin(float t, float B, float Bon, float Boff, float note, int Bsyn)
{
    float Bprog = B-Bon;			// progress within Bar
    float Bproc = Bprog/(Boff-Bon); // relative progress
    float _t = SPB*(B - Bon); // reset time to Bon event
    float f = freqC1(note);

    float env = theta(B-Bon) * theta(Boff-B);
	float sound = clamp(1.1 * _sin(freqC1(note)*t), -0.999,0.999);

    if(Bsyn == -1) return 0.;
    
    if(Bsyn == 0)
    	return env * sound; // test reasons: just give out something simple

    else if(Bsyn == 76)
    {
        env = theta(Bprog) * exp(-(3.*pow(42./note,2.))*Bprog);

        f = f * (1. + .0004*Bprog*_sin(2.*Bprog));
        sound = clamp(shape(sound,10.),-1.,1.);
        sound = .5*sound + .5*TRISQ(t, f, 160, 0.4, 8. + sin(8.*B), 2.2, 0.2, 0.1);
    } 

    else if(Bsyn == 61)
    {   
        env = smoothstep(.0,.0002,Bprog) * smoothstep(.05, 0., B-Boff);

        float filterQ = 20.;
        float filterCO = 2000. + 1000. * env_ADSR(16. * Bprog,Boff-Bon,1.5,2.5,0.2,10.);

        sound = 0.9*macesaw(t, .5*f, filterCO, filterQ, 0.010, 0.020, 0., 0.)
               + .4 * macesaw(t, .499*f, filterCO, filterQ, 0.010, 0.020, 0., 0.);
    }    

    else if(Bsyn == 84)
    {   
        env = smoothstep(.0,.00001,Bprog) * smoothstep(.05, 0., B-Boff);
       
        float filterQ = 100.;
        float filterCO = 1500. + 1000. * smoothstep(0.,0.25,Bproc);

        sound = 0.9*macesaw(t, f, filterCO, filterQ, 0.010, 0.020, 0.3, 3.);
    }        
    
    else if(Bsyn == 75)
    {
        sound = TRISQ(t, freqC1(note), 160, 0.4 + .3 * sin(32.*B * (1.+sin(3.*B))), .1 + .1*sin(24.*B), 0.5, 0.2, 0.1);
    }
//    float QTRISQ(float t, float f, float QUANT, int MAXN, float MIX, float INR, float NDECAY, float RES, float RES_Q)
	else if(Bsyn == 100)
    {
        float dt = 1e-4;
        float b[4] = float[4](0.01, 0.01,0.08,0.2);
        float s[4] = float[4](0.,0.,0.,0.);
        sound = 0.;
        for(int n=0; n<4; n++)
        {
        	s[n] = _saw(_t*f*.5) + .3*_saw(_t*f*(1.-.005*_sin(2.2*B)));
        	_t = _t - dt;
            sound += b[n]*s[n];
        }
    }
    
    return clamp(env,0.,1.) * clamp(sound, -1., 1.);
}


float mainSynth(float time)
{

int NO_trks = 1;
int trk_sep[2] = int[2](0,1);
int trk_syn[1] = int[1](100);
float mod_on[1] = float[1](0.);
float mod_off[1] = float[1](8.);
int mod_ptn[1] = int[1](0);
float mod_transp[1] = float[1](12.);
float inv_NO_tracks = 1.;
float max_mod_off = 8.;
int drum_index = 6;
float drum_synths = 11.;
int NO_ptns = 2;
int ptn_sep[3] = int[3](0,39,39);
float note_on[39] = float[39](0.,.25,.5,.75,1.,1.25,1.5,1.75,2.,2.,2.25,2.5,2.5,2.75,3.,3.,3.25,3.5,3.75,4.,4.25,4.5,4.75,5.,5.25,5.5,5.75,6.,6.,6.25,6.5,6.5,6.75,7.,7.,7.25,7.5,7.5,7.75);
float note_off[39] = float[39](.25,.5,.75,1.,1.25,1.5,1.75,2.,2.25,2.25,2.5,2.75,2.75,3.,3.25,3.25,3.5,3.75,4.,4.25,4.5,4.75,5.,5.25,5.5,5.75,6.,6.25,6.25,6.5,6.75,6.75,7.,7.25,7.25,7.5,7.75,7.75,8.);
float note_pitch[39] = float[39](24.,27.,31.,32.,36.,38.,39.,43.,44.,24.,27.,31.,43.,32.,34.,42.,36.,38.,39.,24.,26.,27.,31.,32.,35.,39.,43.,48.,20.,22.,24.,47.,26.,44.,27.,29.,38.,31.,32.);
float note_vel[39] = float[39](1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.);

float max_release = 1.;
    
float global_norm = .6;
float track_norm[7] = float[7](1.,1.,.9,.9,1.,.9,1.7);
    
	float r = 0.;
    float d = 0.;

    //which beat are we at?
    float BT = mod(BPS * time, max_mod_off); // mod for looping
    if(BT > max_mod_off) return r;

    // drum / sidechaining parameters
    float amt_drum = 0.3;
    float r_sidechain = 1.;
    float amt_sidechain = 0.99;
    float dec_sidechain = 0.6;

    float Bon = 0.;
    float Boff = 0.;

    for(int trk = 0; trk < NO_trks; trk++)
    {
        int TLEN = trk_sep[trk+1] - trk_sep[trk];
        
        int _mod = TLEN;
        for(int i=0; i<TLEN; i++) if(BT < mod_off[(trk_sep[trk]+i)]) {_mod = i; break;}
		if(_mod == TLEN) continue;
        
        float B = BT - mod_on[trk_sep[trk]+_mod];
        
        int ptn = mod_ptn[trk_sep[trk]+_mod];
        int PLEN = ptn_sep[ptn+1] - ptn_sep[ptn];
        
        int _noteU = PLEN-1;
        for(int i=0; i<PLEN-1; i++) if(B < note_on[(ptn_sep[ptn]+i+1)]) {_noteU = i; break;}

        int _noteL = PLEN-1;
        for(int i=0; i<PLEN-1; i++) if(B <= note_off[(ptn_sep[ptn]+i)] + max_release) {_noteL = i; break;}
        
        for(int _note = _noteL; _note <= _noteU; _note++)
        {
            Bon    = note_on[(ptn_sep[ptn]+_note)];
            Boff   = note_off[(ptn_sep[ptn]+_note)];

            if(trk_syn[trk] == drum_index)
            {
                float Bdrum = mod(note_pitch[ptn_sep[ptn]+_note], drum_synths);
                float Bvel = 1.; //note_vel[(ptn_sep[ptn]+_note)] * pow(2.,mod_transp[_mod]/6.);

                float anticlick = 1.-exp(-1000.*(B-Bon));
                float _d = 0.;

            	if(Bdrum < .01) // Sidechain
                {
                    r_sidechain = anticlick - amt_sidechain * theta(B-Bon) * smoothstep(-dec_sidechain,0.,Bon-B);
                }
                else if(Bdrum < 1.01) // Kick1
                {
                    r_sidechain = anticlick - amt_sidechain * theta(B-Bon) * smoothstep(-dec_sidechain,0.,Bon-B);
                    r_sidechain *= 0.5;
                    _d = facekick(B*SPB, Bon*SPB, Bvel);
                }
                else if(Bdrum < 2.01) // Kick2
                {
                    r_sidechain = anticlick - amt_sidechain * theta(B-Bon) * smoothstep(-dec_sidechain,0.,Bon-B);
                    r_sidechain *= 0.5;
                    _d = hardkick(B*SPB, Bon*SPB, Bvel);
                }
                else if(Bdrum < 3.01) // Snare1
                {
                    _d = snare(B*SPB, Bon*SPB, Bvel);
                }
                else if(Bdrum < 4.01) // HiHat
                {
                    _d = hut(B*SPB, Bon*SPB, Bvel);
                }                
                else if(Bdrum < 5.01) // Shake
                {
                    _d = shake(B*SPB, Bon*SPB, Bvel);
                }
                else if(Bdrum < 6.01) // ...
                {
                }          
                d += track_norm[trk] * _d;
            }
            else
            {
                r += track_norm[trk] * distsin(time, B, Bon, Boff,
                                               note_pitch[(ptn_sep[ptn]+_note)] + mod_transp[_mod], trk_syn[trk]);
            }

        }
    }

    d *= global_norm;
    r *= global_norm;

    r_sidechain = 1.;
    amt_drum = .5;

    float snd = s_atan((1.-amt_drum) * r_sidechain * r + amt_drum * d);

    return s_atan(snd);
//    return sign(snd) * sqrt(abs(snd)); // eine von Matzes "besseren" Ideen
}

vec2 mainSound(float t)
{
    //maybe this works in enhancing the stereo feel
    float stereo_width = 0.1;
    float stereo_delay = 0.00001;
    
    //float comp_l = mainSynth(t) + stereo_width * mainSynth(t - stereo_delay);
    //float comp_r = mainSynth(t) + stereo_width * mainSynth(t + stereo_delay);
    
    //return vec2(comp_l * .99999, comp_r * .99999); 
    
    return vec2(mainSynth(t));
}

void main() 
{
   // compute time `t` based on the pixel we're about to write
   // the 512.0 means the texture is 512 pixels across so it's
   // using a 2 dimensional texture, 512 samples per row
   float t = iBlockOffset + ((gl_FragCoord.x-0.5) + (gl_FragCoord.y-0.5)*512.0)/iSampleRate;
    
//    t = mod(t, 4.5);
    
   // Get the 2 values for left and right channels
   vec2 y = iVolume * mainSound( t );

   // convert them from -1 to 1 to 0 to 65536
   vec2 v  = floor((0.5+0.5*y)*65536.0);

   // separate them into low and high bytes
   vec2 vl = mod(v,256.0)/255.0;
   vec2 vh = floor(v/256.0)/255.0;

   // write them out where 
   // RED   = channel 0 low byte
   // GREEN = channel 0 high byte
   // BLUE  = channel 1 low byte
   // ALPHA = channel 2 high byte
   gl_FragColor = vec4(vl.x,vh.x,vl.y,vh.y);
}


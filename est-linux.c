/* Electronic Strings by Team210
 * Copyright (C) 2018  Alexander Kraus <nr4@z10.info>
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

#include <X11/X.h>
#include <X11/Xlib.h>
#include <X11/Xatom.h>

#include <GL/gl.h>
#include <GL/glx.h>
#include "glext.h"

#include <alsa/asoundlib.h>
#define PCM_DEVICE "default"

#include <pthread.h>

#include <time.h>
#include <sys/time.h>

#include <dlfcn.h>

// OpenGL extensions
PFNGLGETPROGRAMIVPROC glGetProgramiv;
PFNGLGETSHADERIVPROC glGetShaderiv;
PFNGLGETSHADERINFOLOGPROC glGetShaderInfoLog;
PFNGLCREATESHADERPROC glCreateShader;
PFNGLCREATEPROGRAMPROC glCreateProgram;
PFNGLSHADERSOURCEPROC glShaderSource;
PFNGLCOMPILESHADERPROC glCompileShader;
PFNGLATTACHSHADERPROC glAttachShader;
PFNGLLINKPROGRAMPROC glLinkProgram;
PFNGLUSEPROGRAMPROC glUseProgram;
PFNGLGETUNIFORMLOCATIONPROC glGetUniformLocation;
PFNGLUNIFORM2FPROC glUniform2f;
PFNGLUNIFORM1FPROC glUniform1f;
PFNGLGENFRAMEBUFFERSPROC glGenFramebuffers;
PFNGLBINDFRAMEBUFFERPROC glBindFramebuffer;
PFNGLFRAMEBUFFERTEXTURE2DPROC glFramebufferTexture2D;
PFNGLNAMEDRENDERBUFFERSTORAGEEXTPROC glNamedRenderbufferStorageEXT;

// Shader globals
int w = 1366, h = 768; // TODO: add a way of configuring this in the future.
int lb_program, lb_progress_location, lb_resolution_location, lb_time_location;
int gfx_program, gfx_time_location, gfx_resolution_location;
int sfx_program, sfx_blockoffset_location, sfx_samplerate_location, sfx_volumelocation;

// Demo globals
double t_start = 0.;
int music_thread_handle;
pthread_t music_thread;
int sample_rate = 44100, channels = 2;
double duration1 = 312.*.43; //3 min running time
float *smusic1;
int music1_size;
#define texs 512
int block_size = texs*texs;
snd_pcm_t *snd_handle;

void MusicThread()
{
    snd_pcm_writei(snd_handle, smusic1, 4*music1_size);
}

//TODO: remove
void debug(int shader_handle)
{
    int compile_status = 0;
    glGetShaderiv(shader_handle, GL_COMPILE_STATUS, &compile_status);
    if(compile_status != GL_TRUE)
    {
        printf("FAILED.\n");
        int len = 12;
        glGetShaderiv(shader_handle, GL_INFO_LOG_LENGTH, &len);
        printf("log length: %d\n", len);
        GLchar CompileLog[1024];
        glGetShaderInfoLog(shader_handle, len, NULL, CompileLog);
        printf("error: %s\n", CompileLog);
    }
}

int main(int argc, char **args)
{
    XInitThreads();
    
    Display                 *display;
    Window                  root;
    GLint                   att[] = { GLX_RGBA, GLX_DEPTH_SIZE, 24, GLX_DOUBLEBUFFER, None };
    XVisualInfo             *vi;
    Colormap                cmap;
    XSetWindowAttributes    swa;
    Window                  win;
    GLXContext              glc;
    XWindowAttributes       gwa;
    XEvent                  xevent;

    display = XOpenDisplay(NULL);
    root = DefaultRootWindow(display);
    vi = glXChooseVisual(display, 0, att);
    cmap = XCreateColormap(display, root, vi->visual, AllocNone);
    swa.colormap = cmap;
    
    swa.event_mask = ExposureMask | KeyPressMask;
    win = XCreateWindow(display, root, 0, 0, w, h, 0, vi->depth, InputOutput, vi->visual, CWColormap | CWEventMask, &swa);
 
    Atom atoms[2] = { XInternAtom(display, "_NET_WM_STATE_FULLSCREEN", True), None };
    XChangeProperty(
        display, 
        win, 
        XInternAtom(display, "_NET_WM_STATE", True),
                    XA_ATOM, 32, PropModeReplace,(unsigned char*) atoms, 1
    );
    XMapWindow(display, win);
    glc = glXCreateContext(display, vi, NULL, GL_TRUE);
    glXMakeCurrent(display, win, glc);

    // OpenGL extensions
    void *gl = (void*)dlopen("libGL.so", RTLD_LAZY | RTLD_GLOBAL);
    glGetProgramiv = (PFNGLGETPROGRAMIVPROC) dlsym(gl, "glGetProgramiv");
    glGetShaderiv = (PFNGLGETSHADERIVPROC) dlsym(gl, "glGetShaderiv");
    glGetShaderInfoLog = (PFNGLGETSHADERINFOLOGPROC) dlsym(gl, "glGetShaderInfoLog");
    glCreateShader = (PFNGLCREATESHADERPROC) dlsym(gl, "glCreateShader");
    glCreateProgram = (PFNGLCREATEPROGRAMPROC) dlsym(gl, "glCreateProgram");
    glShaderSource = (PFNGLSHADERSOURCEPROC) dlsym(gl, "glShaderSource");
    glCompileShader = (PFNGLCOMPILESHADERPROC) dlsym(gl, "glCompileShader");
    glAttachShader = (PFNGLATTACHSHADERPROC) dlsym(gl, "glAttachShader");
    glLinkProgram = (PFNGLLINKPROGRAMPROC) dlsym(gl, "glLinkProgram");
    glUseProgram = (PFNGLUSEPROGRAMPROC) dlsym(gl, "glUseProgram");
    glGetUniformLocation = (PFNGLGETUNIFORMLOCATIONPROC) dlsym(gl, "glGetUniformLocation");
    glUniform2f = (PFNGLUNIFORM2FPROC) dlsym(gl, "glUniform2f");
    glUniform1f = (PFNGLUNIFORM1FPROC) dlsym(gl, "glUniform1f");
    glGenFramebuffers = (PFNGLGENFRAMEBUFFERSPROC) dlsym(gl, "glGenFramebuffers");
    glBindFramebuffer = (PFNGLBINDFRAMEBUFFERPROC) dlsym(gl, "glBindFramebuffer");
    glFramebufferTexture2D = (PFNGLFRAMEBUFFERTEXTURE2DPROC) dlsym(gl, "glFramebufferTexture2D");
    glNamedRenderbufferStorageEXT = (PFNGLNAMEDRENDERBUFFERSTORAGEEXTPROC) dlsym(gl, "glNamedRenderbufferStorage");
        
    // Init gfx
#undef VAR_IRESOLUTION
#undef VAR_ITIME
#include "gfx.h"
#ifndef VAR_IRESOLUTION
#define VAR_IRESOLUTION "iResolution"
#ifndef VAR_ITIME
#define VAR_ITIME "iTime"
    int gfx_size = strlen(gfx_frag),
        gfx_handle = glCreateShader(GL_FRAGMENT_SHADER);
    gfx_program = glCreateProgram();
    glShaderSource(gfx_handle, 1, (GLchar **)&gfx_frag, &gfx_size);
    glCompileShader(gfx_handle);
    glAttachShader(gfx_program, gfx_handle);
    glLinkProgram(gfx_program);
    glUseProgram(gfx_program);
    gfx_time_location =  glGetUniformLocation(gfx_program, VAR_ITIME);
    gfx_resolution_location = glGetUniformLocation(gfx_program, VAR_IRESOLUTION);
#endif
#endif
    
    debug(gfx_handle);
    
    // Init sfx
#undef VAR_IBLOCKOFFSET
#undef VAR_ISAMPLERATE
#undef VAR_IVOLUME
#include "sfx.h"
#ifndef VAR_IVOLUME
#define VAR_IVOLUME "iVolume"
#ifndef VAR_ISAMPLERATE
#define VAR_ISAMPLERATE "iSampleRate"
#ifndef VAR_IBLOCKOFFSET
#define VAR_IBLOCKOFFSET "iBlockOffset"
    int sfx_size = strlen(sfx_frag),
        sfx_handle = glCreateShader(GL_FRAGMENT_SHADER);
    sfx_program = glCreateProgram();
    glShaderSource(sfx_handle, 1, (GLchar **)&sfx_frag, &sfx_size);
    glCompileShader(sfx_handle);
    glAttachShader(sfx_program, sfx_handle);
    glLinkProgram(sfx_program);
    glUseProgram(sfx_program);
    sfx_samplerate_location = glGetUniformLocation(sfx_program, VAR_ISAMPLERATE);
    sfx_blockoffset_location = glGetUniformLocation(sfx_program, VAR_IBLOCKOFFSET);
    sfx_volumelocation = glGetUniformLocation(sfx_program, VAR_IVOLUME);
#endif
#endif
#endif
    
    int nblocks1 = sample_rate*duration1/block_size+1;
    music1_size = nblocks1*block_size; 
    smusic1 = (float*)malloc(4*music1_size);
    
    unsigned int snd_framebuffer;
    glGenFramebuffers(1, &snd_framebuffer);
    glBindFramebuffer(GL_FRAMEBUFFER, snd_framebuffer);
    glPixelStorei(GL_PACK_ALIGNMENT,  4);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 4);
    
    unsigned int snd_texture;
    glGenTextures(1, &snd_texture);
    glBindTexture(GL_TEXTURE_2D, snd_texture);
    glTexImage2D(GL_TEXTURE_2D, 0,GL_RGBA, texs, texs, 0, GL_RGBA, GL_BYTE, 0);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);

    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, snd_texture, 0);

    // Render sfx
    for(int i=0; i<nblocks1; ++i)
    {
        double tstart = (double)(i*block_size)/(double)sample_rate;
        
        glViewport(0,0,texs,texs);
        
        glUniform1f(sfx_volumelocation, 1.);
        glUniform1f(sfx_samplerate_location, (float)sample_rate);
        glUniform1f(sfx_blockoffset_location, (float)tstart);
        
        glBegin(GL_QUADS);
        glVertex3f(-1,-1,0);
        glVertex3f(-1,1,0);
        glVertex3f(1,1,0);
        glVertex3f(1,-1,0);
        glEnd();

        glFlush();

        glReadPixels(0, 0, texs, texs, GL_RGBA, GL_BYTE, smusic1+i*block_size);
    }
    
    // Reset everything for rendering gfx again
    glViewport(0, 0, w, h);
    glUseProgram(0);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);   
    glUseProgram(gfx_program);
    glUniform2f(gfx_resolution_location, w, h);
    
    // Get start time for relative time sync
    struct timeval tv;
    gettimeofday(&tv, NULL);
    t_start = (double)tv.tv_sec+(double)tv.tv_usec/1.e6;
  
    snd_pcm_hw_params_t *params;
    snd_pcm_uframes_t frames;

    snd_pcm_open(&snd_handle, PCM_DEVICE, SND_PCM_STREAM_PLAYBACK, 0);
    snd_pcm_set_params(snd_handle, SND_PCM_FORMAT_S16_LE, SND_PCM_ACCESS_RW_INTERLEAVED, channels, sample_rate, 0, (music1_size)*channels);
    
    music_thread_handle = pthread_create(&music_thread, NULL, (void*)&MusicThread, 0);
    
    int x_file_descriptor = ConnectionNumber(display);
    fd_set x_file_descriptor_set;
    
    // Main loop
    while(1)
    {
        // Exit demo if any key is pressed.
        while(XPending(display))
        {
            XNextEvent(display, &xevent);
        
            if(xevent.type == KeyPress) 
            {
                exit(0);
            }
        }
        
        FD_ZERO(&x_file_descriptor_set);
        FD_SET(x_file_descriptor, &x_file_descriptor_set);
        
        struct timeval t;
        t.tv_usec = 1.e6/60.;
        t.tv_sec = 0;
        
        int num_ready_fds = select(x_file_descriptor + 1, &x_file_descriptor_set, NULL, NULL, &t);
        if (num_ready_fds == 0)    
        {            
            struct timeval tv_now;
            gettimeofday(&tv_now, NULL);
            double t_now = (double)tv_now.tv_sec+(double)tv_now.tv_usec/1.e6;
            
            if(t_now-t_start > duration1) 
            {
                exit(0);
            }
            
            glUniform1f(gfx_time_location, t_now-t_start);
            
            glBegin(GL_QUADS);
            
            glVertex3f(-1,-1,0);
            glVertex3f(-1,1,0);
            glVertex3f(1,1,0);
            glVertex3f(1,-1,0);
            
            glEnd();

            glFlush();
            
            glXSwapBuffers(display, win);
        }
    }
    exit(0);
    return 0;
}



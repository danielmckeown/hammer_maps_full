import numpy as np
import ctypes
import math
import os.path
import struct
import array

def checklen(x):
    return len(np.array(x,ndmin=1));
def fcor(x):
    return np.array(x,dtype='f',ndmin=1)
def vfloat(x):
    return x.ctypes.data_as(ctypes.POINTER(ctypes.c_float));

def ok_scan(input,xmax=1.0e10,pos=0):
    if (pos==1):
        return (np.isnan(input)==False) & (abs(input)<=xmax) & (input > 0.);
    else:
        return (np.isnan(input)==False) & (abs(input)<=xmax);
    
    
    
    
def get_particle_hsml( x, y, z, vx, vy, vz, DesNgb=32, Hmax=0.):
    
    
    x=fcor(x); y=fcor(y); z=fcor(z); N=checklen(x); 

    vx=fcor(vx); vy=fcor(vy); vz=fcor(vz); N=checklen(x); 
    
    ok=(ok_scan(x) & ok_scan(y) & ok_scan(z));
    ok= ok & (ok_scan(vx) & ok_scan(vy) & ok_scan(vz));

    x=x[ok]; y=y[ok]; z=z[ok];
    vx=vx[ok]; vy=vy[ok]; vz=vz[ok];

    
    if(Hmax==0.):
        dx=np.max(x)-np.min(x); dy=np.max(y)-np.min(y); dz=np.max(z)-np.min(z); ddx=np.max([dx,dy,dz]);
        Hmax=5.*ddx*(np.float(N)**(-1./3.)); ## mean inter-particle spacing                                                                                                                                                                                   

  
  ##### now accessing the C code

    import inspect
    exec_call = "/Users/danmckeown/repos/FIRE_studio/firestudio/utils/stellar_utils/c_libraries/StellarHsml/starhsml.so"
    
    # calling program in the dir.
    h_routine=ctypes.cdll[exec_call];

    # These are python variables which we are casting as C variables
    h_out_cast=ctypes.c_float*N; H_OUT=h_out_cast();
    v_out_cast=ctypes.c_float*N; V_OUT=v_out_cast();
    v_mag_cast=ctypes.c_float*N; V_MAG=v_mag_cast();
    # First variable name needs to match last
    
    # These lines create the pointers that are passed to the C code  ( ctypes.byref(H_OUT), ctypes.byref(V_OUT))
    
    
    ## main call to the hsml-finding routine
    h_routine.stellarhsml( ctypes.c_int(N), \
    ## load the routine we need                                                                                                                                                                                                                               
    
        vfloat(x), vfloat(y), vfloat(z), \
        vfloat(vx), vfloat(vy), vfloat(vz), ctypes.c_int(DesNgb), \
        ctypes.c_float(Hmax), ctypes.byref(H_OUT), ctypes.byref(V_OUT),ctypes.byref(V_MAG)
    )                       
                          
    ## now put the output arrays into a useful format 
    
    # These two lines below get the pointers out
    h = np.ctypeslib.as_array(np.copy(H_OUT));                     
    vdisp = np.ctypeslib.as_array(np.copy(V_OUT));
    vmag = np.ctypeslib.as_array(np.copy(V_MAG));
    
    return h, vdisp,vmag 
                          
   

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "proto.h"
#include "ngbtree3d.h"

void allocate_3d(int Nsize);
void set_particle_pointer(int Nsize, float* x, float* y, float* z);
void free_memory_3d(int Nsize);

struct particle_3d 
{
  float Pos[3];
} **P3d;

float velocity_dispersion( float vx[], float vy[], float vz[], int n , float * abdx, float * abdy, float * abdz, float * final_abs,  float * vel_disp_radial   )
{ 

    int i;
    
    float total_vx,total_vy,total_vz ; 
    float total_radial_v;
    float Average_vx,Average_vy,Average_vz,Average_difference;
    float diff_vx[n] , abs_dif_vx[n] ,diff_vy[n], abs_dif_vy[n] , diff_vz[n] , abs_dif_vz[n];
    float radial_sq_v[n] , sq_diff_vx[n], sq_abs_dif_vx[n], sq_diff_vy[n], sq_abs_dif_vy[n], sq_diff_vz[n], sq_abs_dif_vz[n]  ;

    float total_diff_vx, total_diff_vy, total_diff_vz ;
    float vel_disp_sq_vx, vel_disp_sq_vy , vel_disp_sq_vz , vel_disp_radial_sq_v ;
    
    float radial_sq_root_abs_diff_v_mag[n] , radial_sq_abs_diff_v_mag_final[n];

    float total_abs_diff_vx_sq, total_abs_diff_vy_sq, total_abs_diff_vz_sq; 
    
    float total_abdx, total_abdy, total_abdz;
    
    float total_abs_radial_v_sq_root ;
    float abs_vel_radial_sq_root_av ;
    float radial_sq_v_mag[n];
    float total_abs_radial_v;
    float abs_mag_average;
    total_vx = 0;
    total_vy = 0;
    total_vz = 0;
    
    total_diff_vx = 0;
    total_diff_vy = 0;
    total_diff_vz = 0;
    total_abdx = 0;
    total_abdy = 0;
    total_abdz = 0;
    
    total_radial_v = 0;
    total_abs_radial_v = 0;
    total_abs_radial_v_sq_root = 0;
    

    for(i = 0; i < n; i++) {
       total_vx += vx[i];
       total_vy += vy[i];
       total_vz += vz[i];
    }
   

    Average_vx = total_vx /(float)n ;
    Average_vy = total_vy/(float)n ; 
    Average_vz = total_vz/(float)n ;

   

    for ( i = 0; i < n ; i++ ) {
      
      
        diff_vx[i] = vx[i] - Average_vx ; 
        sq_diff_vx[i] = diff_vx[ i ] * diff_vx[ i ] ;

        abs_dif_vx[i] = vx[0] - vx[i] ;
        sq_abs_dif_vx[i] = abs_dif_vx[i] * abs_dif_vx[i] ;

        diff_vy[ i ] = vy[i] - Average_vy; 
        sq_diff_vy[i] = diff_vy[ i ] * diff_vy[ i ] ;

        abs_dif_vy[i] = vy[0] - vy[i] ;
        sq_abs_dif_vy[i] = abs_dif_vy[i] * abs_dif_vy[i] ;
        
        diff_vz[ i ] = vz[i] - Average_vz; 
        sq_diff_vz[i] = diff_vz[ i ] * diff_vz[ i ] ;

        abs_dif_vz[i] = vz[0] - vz[i] ;
        sq_abs_dif_vz[i] = abs_dif_vz[i] * abs_dif_vz[i] ;
        
        radial_sq_v[i] = sq_diff_vx[i] + sq_diff_vy[i] + sq_diff_vz[i] ;
        
        radial_sq_v_mag[i] = sq_abs_dif_vx[i] + sq_abs_dif_vy[i] + sq_abs_dif_vz[i]; 
        
        radial_sq_root_abs_diff_v_mag[i] = sqrt(radial_sq_v_mag[i]) ;
   
   
   }

    
    for (int i = 0; i < n ; i++) {
     
     total_diff_vx += sq_diff_vx[i];

    }
    
    
    for (int i = 0; i < n ; i++) {
     
     total_abdx += sq_abs_dif_vx[i];

    }
    
    
    for (int i = 0; i < n ; i++) {
     
     total_abdy += sq_abs_dif_vy[i];

    }


    for (int i = 0; i < n ; i++) {
     
     total_abdz += sq_abs_dif_vz[i];

    }

    
    for (int i = 0; i < n ; i++) {
     
     total_diff_vy += sq_diff_vy[i];
    }

    for (int i = 0; i < n ; i++) {
     
     total_diff_vz += sq_diff_vz[i];
    
    }
    
    for (int i = 0; i < n ; i++) {
     
     total_radial_v += radial_sq_v[i];
    }

   
   for (int i = 0; i < n ; i++) {
     
     total_abs_radial_v  +=  radial_sq_v_mag[i];
    }



    for (int i = 0; i < n ; i++) {
     
     total_abs_radial_v_sq_root  += radial_sq_root_abs_diff_v_mag[i];
    }



    vel_disp_sq_vx = total_diff_vx / (float)n ;
    vel_disp_sq_vy = total_diff_vy / (float)n;
    vel_disp_sq_vz = total_diff_vz / (float)n;
    vel_disp_radial_sq_v = total_radial_v / (float)n  ;
    
    abs_mag_average = total_abs_radial_v / (float)n ;
    
    abs_vel_radial_sq_root_av  = total_abs_radial_v_sq_root / (float)31.0 ;
    
    
    *abdx = total_abdx / (float)n;
    *abdy = total_abdy / (float)n;
    *abdz = total_abdz / (float)n;
    *final_abs = abs_vel_radial_sq_root_av; 
   
    *vel_disp_radial = sqrt(vel_disp_radial_sq_v) ;
    
    
    return 0; 

}

///   THIS IS THE MAIN CALLING FUNTION  

// revised call for python calling //
//int stellarhsml(int argc,void *argv[])

// This gets called once for all the particles, so all x,y,z, vx , vy , vz are pulled in at once
// int N_in tells number of particles feeding the function

int stellarhsml(int N_in, float* x, float* y, float* z, float* vx, float* vy, float* vz, int DesNgb, float Hmax, float* H_OUT, float* V_MAG, float* V_OUT)
{
  float h_guess, h2, xyz[3], dummy[3], h_guess_0;
  int i, j, ind, ngbfound;
  allocate_3d(N_in);
  // allocate_3d(N_in) allocates the mememory necessary to store the particles I have in a tree
  // This sets up the memory for a tree that stores the neighbors of particles
  float *r2list;
  // list of distances squared to neighboring particles
  int *ngblist;
  // list of ids
  //float vdisp;
  

  float vx_ngb[DesNgb], vy_ngb[DesNgb], vz_ngb[DesNgb];
  // sets up 3 arrays,  which is the velocity for the nearest neighbors vx,vy,vz
  float abdx,abdy,abdz,final_abs, vel_disp_radial;
  // single value which is passed in and altered and replaced for each individual particle
  
  // comment out for now 
  //printf("N=%d\n",N_in); printf("Hmax=%g\n",Hmax); printf("DesNgb=%d\n",DesNgb);

  ngb3d_treeallocate(N_in, 2*N_in);
  // similar purpose for tree , more memory setup
  set_particle_pointer(N_in, x,y,z);
  // setting up a pointer for the tree
  ngb3d_treebuild((float **)&P3d[1], N_in, 0, dummy, dummy);
  // Hopkins Springel code
  h_guess = Hmax/150.0e0; h_guess_0=h_guess;
  // initial guess for smoothing len
  
  for(i=0;i<N_in;i++)
  {
	  xyz[0]=P3d[i+1]->Pos[0]+1.0e-10;
	  xyz[1]=P3d[i+1]->Pos[1]+1.0e-10;
	  xyz[2]=P3d[i+1]->Pos[2]+1.0e-10;
	  // looping through and getting x, y, z pos. of all the particles
    h2=ngb3d_treefind( xyz, DesNgb ,1.04*h_guess, &ngblist, &r2list, Hmax, &ngbfound); 
    // telling us the smoothing len guess,  &gnblist, &r2list are indices and number for nn, 
    // next in loop below it loops through list of neighbors, and it has stored the velocities of all the
    // neighbors to vx, vy, vz
    for (j=0;j<ngbfound;j++)
    {
      ind = ngblist[j];
      vx_ngb[j] = vx[ind];
      vy_ngb[j] = vy[ind];
      vz_ngb[j] = vz[ind];
      // printf("j = %d, ind = %d \n", j, ind ); fflush(stdout);
      // printf("x, y, z =  %g, %g, %g \n", x[ind], y[ind], z[ind] ); fflush(stdout);
      // printf("vx, vy, vz =  %g, %g, %g \n", vx[ind], vy[ind], vz[ind] ); fflush(stdout);
      // printf("vx, vy, vz =  %g, %g, %g \n", vx_ngb[j], vy_ngb[j], vz_ngb[j] ); fflush(stdout);
    }

    // vdisp = velocity_dispersion( &vx_ngb, &vy_ngb, &vz_ngb, ngbfound );
     // make sure to preseverve order
       velocity_dispersion( &vx_ngb, &vy_ngb, &vz_ngb, ngbfound,&abdx,&abdy, &abdz,  &final_abs , &vel_disp_radial );

    
    // ///   New Lines //   ///   ///    ///    ///   //   ///   ///    ////   ////    ////    ///// 
    //if (i < 5) {

        //printf("%.6f\n",  mag);
       // printf("%.6f\n",  final_abs);
        
        //printf("%.6f\n",  abdx);
        //printf("%.6f\n",  abdy);
        //printf("%.6f\n",  abdz);
    
    //}

    // comment out for now ///////
    
    //if(!(i%100000))
    if(i < 6)
    
    // why are there no spaces needed between entries, what does %d %g mean?
    {
    printf("i=%d hmax=%g h_guess=%g h=%g xyz=%g|%g|%g ngb=%d magnitude=%g vel_disp=%g  \n",
        i,Hmax,h_guess,sqrt(h2),xyz[0],xyz[1],xyz[2],ngbfound,final_abs,vel_disp_radial); fflush(stdout);
    }
      
      
      H_OUT[i] = sqrt(h2); 
      V_MAG[i] = final_abs;
      V_OUT[i] = vel_disp_radial;
      
      h_guess = H_OUT[i]; // use this value for next guess, should speed things up // 
      //if (h_guess>10.*h_guess_0) h_guess=2.*h_guess_0;
    } 

  ngb3d_treefree();
  free_memory_3d(N_in);
  printf("done\n");
  return 0;
}



void set_particle_pointer(int Nsize, float* x, float* y, float* z)
{
  int i;
  float *pos;
  for(i=1;i<=Nsize;i++)
    {
      P3d[i]->Pos[0] = x[i-1];
      P3d[i]->Pos[1] = y[i-1];
      P3d[i]->Pos[2] = z[i-1];
    }
}


void allocate_3d(int Nsize)
{
  printf("allocating memory...\n");
  int i;
  if(Nsize>0)
    {
      if(!(P3d=malloc(Nsize*sizeof(struct particle_3d *))))
	{
	  printf("failed to allocate memory. (A)\n");
	  exit(0);
	}
      P3d--;   /* start with offset 1 */
      if(!(P3d[1]=malloc(Nsize*sizeof(struct particle_3d))))
	{
	  printf("failed to allocate memory. (B)\n");
	  exit(0);
	}
      for(i=2;i<=Nsize;i++)   /* initiliaze pointer table */
	P3d[i]=P3d[i-1]+1;
    }
  printf("allocating memory...done\n");
}


void free_memory_3d(int Nsize)
{
  if(Nsize>0)
    {
      free(P3d[1]);
      P3d++;
      free(P3d);
    }
}



// Always recompile using make clean  and then make in the directory where the file is saved ....

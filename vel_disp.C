#include <stdio.h>
#include <math.h>

// function to return velocity dispersion
// in an array of size length_of_v, this just works on a single 32 nn pair , it is then incorporated for the larger 
// system of particles when put into the python wrapper


float  velocity_dispersion( float vx[], float vy[], float vz[], int n ,   float * final_abs, float * vel_disp_radial  )

{ 

  
    int i;
    
    float total_vx,total_vy,total_vz ; 
    float total_radial_v;
    float Average_vx,Average_vy,Average_vz,Average_difference; 
    float diff_vx[n] , diff_vy[n], diff_vz[n]  ;
    float abs_dif_vx[n] , abs_dif_vy[n] , abs_dif_vz[n];
    float radial_sq_v[n] , sq_diff_vx[n], sq_abs_dif_vx[n], sq_diff_vy[n], sq_abs_dif_vy[n], sq_diff_vz[n], sq_abs_dif_vz[n]  ;
    float total_diff_vx, total_diff_vy, total_diff_vz ;
    
    float vel_disp_sq_vx, vel_disp_sq_vy , vel_disp_sq_vz , vel_disp_radial_sq_v ;
    
    float radial_sq_root_abs_diff_v_mag[n] , radial_sq_abs_diff_v_mag_final[n];

    
    float total_abs_diff_vx_sq, total_abs_diff_vy_sq, total_abs_diff_vz_sq , total_abs_radial_v_sq_root ;
    float abs_vel_radial_sq_root_av ;
    total_vx = 0;
    total_vy = 0;
    total_vz = 0;
   
   
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
        
        abs_dif_vx[i] = vx[0] - vx[i] ;
        // square the components to take the magnitude
        sq_abs_dif_vx[i] = abs_dif_vx[i] * abs_dif_vx[i] ;

        sq_diff_vx[i] = diff_vx[ i ] * diff_vx[ i ] ;
      
        diff_vy[ i ] = vy[i] - Average_vy; 
        
        abs_dif_vy[i] = vy[0] - vy[i] ;
        
        sq_abs_dif_vy[i] = abs_dif_vy[i] * abs_dif_vy[i] ;
        

        sq_diff_vy[i] = diff_vy[ i ] * diff_vy[ i ] ;
      
        diff_vz[ i ] = vz[i] - Average_vz; 
      
        abs_dif_vz[i] = vz[0] - vz[i] ;

        sq_abs_dif_vz[i] = abs_dif_vz[i] * abs_dif_vz[i] ;
        

        sq_diff_vz[i] = diff_vz[ i ] * diff_vz[ i ] ;
   
        radial_sq_v[i] =  sq_diff_vx[i] + sq_diff_vy[i] + sq_diff_vz[i] ;

        
        radial_sq_root_abs_diff_v_mag[i] = sqrt (sq_abs_dif_vx[i] + sq_abs_dif_vy[i] + sq_abs_dif_vz[i] ) ;

        //radial_sq_abs_diff_v_mag_final[i] = radial_sq_abs_diff_v_mag[i] * radial_sq_abs_diff_v_mag[i] ;

       
        printf("%.6f\n", radial_sq_v[i]);
        
        printf("%.6f\n", radial_sq_root_abs_diff_v_mag[i]);

        //printf("%.6f\n", radial_sq_abs_diff_v_mag_final[i] );

   }

    for (int i = 0; i < n ; i++) {
     
     total_diff_vx += sq_diff_vx[i];

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

    printf("%.6f\n", total_radial_v );

    
    for (int i = 0; i < n ; i++) {
     
     total_abs_radial_v_sq_root  += radial_sq_root_abs_diff_v_mag[i];
    }




    vel_disp_sq_vx = total_diff_vx / (float)n ;
    vel_disp_sq_vy = total_diff_vy / (float)n;
    vel_disp_sq_vz = total_diff_vz / (float)n;
    
    

    printf("%.6f\n",vel_disp_radial_sq_v );
    printf("%.6f\n", total_abs_radial_v_sq_root/ (float)n );

       
    abs_vel_radial_sq_root_av  = total_abs_radial_v_sq_root / (float)n ;
    
    *final_abs = abs_vel_radial_sq_root_av ; 
    
    vel_disp_radial_sq_v = total_radial_v / (float)n  ;

    *vel_disp_radial = sqrt(vel_disp_radial_sq_v) ;
    
    return 0;
    
   
}



int main (void) {


float vx[] = { 7,59, 46 , 12, 96, 8. , 46, 76 , 68, 87,107,78, 91,34, 107, 60, 110 ,70, 86 ,115 ,42 ,45, 99 ,56,59 ,39,51,121,45,33,118 ,113 };
float vy[] = { 49,31,80,43 ,46.,32,18,16 ,41, 25,36 ,21 ,61 ,80 ,50 ,43,73,25,108,49,86 ,22 ,30,121,37 ,88 ,95,56,123 ,80 ,53,83};
float vz[] = {73,73, 25, 68, 91, 34, 9,64,107,44, 22, 44, 91, 85, 99, 106, 17, 61, 42, 22, 69, 74, 50, 104, 89, 61, 41, 28, 33, 127, 35, 97};

float  final_abs, vel_disp_radial;

// keep all this stuff later 


// float radial_vel_disp;

 
 velocity_dispersion(vx, vy, vz,  32. , &final_abs , &vel_disp_radial  ) ;
 
  printf("%.6f\n",  final_abs * final_abs);
  
  printf( "%.6f\n",   2 * vel_disp_radial * vel_disp_radial );
 
 
  // printf("%.6f\n", abs_vel_radial_sq_av * abs_vel_radial_sq_av );
  return 0;
}




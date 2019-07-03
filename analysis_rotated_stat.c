#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_histogram.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include <math.h>


#define TAU 20
#define DIAM 20
#define DIST 6
#define GAIN 8
#define SCELTA 1
#define N_RIP 115
#define N_RIP0 1
#define N_BIN 20
#define N_BIN_CALIB 20
#define SUBTRACT_DIPOLE 0

#define OLD_METHOD 0
#define ALONG_CELL_AXIS 1
#define ALONG_DIPOLE_AXIS 0

#define N_TRIALS 100


struct vector {

    double x;
    double y;

};

struct vector_int {

    int x;
    int y;

};

struct dipole_term {

    double xx;
    double xy;
    double yx;
    double yy;

};

struct quadrupole_term {

    double xxx;
    double xxy;
    double xyx;
    double xyy;
    double yxx;
    double yxy;
    double yyx;
    double yyy;

};



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// DECLARE FUNCTIONS
//
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void projection_cell_axis(struct vector *, struct vector * , struct vector * , struct vector * ,   struct vector * , struct vector * , double *, double *,  int );


void project_traction(struct vector * , struct vector * , struct vector * , struct vector * , int );

void get_polar_coords(struct vector * ,  double * ,  double * , int );


void project_position(struct vector * , struct vector * , struct vector * , struct vector * , struct vector * ,  int );

void get_dipole(struct vector * , struct vector  , struct vector * ,  struct dipole_term *, int,  int , int );


void get_quadrupole(struct vector * , struct vector , struct vector * , struct dipole_term * , struct quadrupole_term * , int , int, int );


void get_dipole_axis(double * , double * , struct dipole_term * , double * , double * , double *, double * );
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




int main()
{
    FILE *fptr;
    FILE *out_dip ,*out_quad, *out_omega, *out_tract, *out_max_eigenv,  *out_angle_eigenv,  *out_angle_eigenv_axis, *out_angle_cell_speed, *out_speed_cell, *out_angle_cell_speed_axis, *out_monop, *out_cell_pos;

    struct vector *cell_traject, *data_vel, *dipole_center;

    struct vector *cellpos, *dist, *cellvel, *vel, *vel_neg[N_RIP], *vel_pos[N_RIP], cell_initial_pos, u, u_perp;
    struct vector *pos_dip_axis, *vel_dip_axis;
    struct vector mediapos, meanvel_neg, variancevel_neg, meanvel_pos, variancevel_pos;
    struct vector  *data_pos, *pos;
    struct vector_int sizes_neg[N_RIP], sizes_pos[N_RIP], count_pos, count_neg;

    struct dipole_term dipole;
    struct quadrupole_term quadrupole;

    int sizes[N_RIP];
    int size, sizenew, sizec, half, countereffective;
    int screeing;

    double dipolemix;
    double *velX, *velY;


    double meanmagnitudestress, variancemagnitudestress, meanvelx,meanvely, variancevelx, variancevely;
    double *cellanglepos, *cellanglevel, *cellanglevel_axis;

    double meanmagnitudepos, variancemagnitudepos;
    double screening, screening_quad;

    double *anglepos, *anglevel, *magnvel, *magnpos;
    double trace, determ,  max_eigenv, min_eigenv, angle_max_eigenv, axis_max_eigenv, angle_min_eigenv, axis_min_eigenv;

    int a, k, n, i;
    int dummy_int;
    double Omega;
    double dummy_float;

    char file_name_centre[550], file_name_pos[550], file_name_vel[550], file_name_dip[550], file_name_quad[550],  file_name_omega[550], file_name_tract[550], file_name_max_eigenv[550], file_name_angle_eigenv[550], file_name_angle_eigenv_axis[550], file_name_angle_cell_speed[550], file_name_speed_cell[550], file_name_angle_cell_speed_axis[550], file_name_monop[550], file_name_cell_pos[550];



    /*** output files ***/
    sprintf(file_name_cell_pos, "results-analysis/cell_pos_new.dat");
    out_cell_pos=fopen(file_name_cell_pos,"w");

    sprintf(file_name_omega, "results-analysis/omega.dat");
    out_omega = fopen(file_name_omega,"w");

    sprintf(file_name_tract, "results-analysis/magnitude_traction.dat");
    out_tract = fopen(file_name_tract,"w");

    sprintf(file_name_max_eigenv, "results-analysis/max_eigen.dat");
    out_max_eigenv = fopen(file_name_max_eigenv,"w");

    sprintf(file_name_angle_eigenv, "results-analysis/angles_eigen.dat");
    out_angle_eigenv = fopen(file_name_angle_eigenv,"w");

    sprintf(file_name_angle_eigenv_axis, "results-analysis/axis_eigen.dat");
    out_angle_eigenv_axis = fopen(file_name_angle_eigenv_axis,"w");

    sprintf(file_name_speed_cell, "results-analysis/cell_speed.dat");
    out_speed_cell = fopen(file_name_speed_cell,"w");

    sprintf(file_name_angle_cell_speed, "results-analysis/angles_cell_speed.dat");
    out_angle_cell_speed = fopen(file_name_angle_cell_speed,"w");

    sprintf(file_name_angle_cell_speed_axis, "results-analysis/axis_cell_speed.dat");
    out_angle_cell_speed_axis = fopen(file_name_angle_cell_speed_axis,"w");

    
    int aa;
    
    for( aa = 90; aa < N_TRIALS; aa++){
        
        sprintf(file_name_dip, "results-analysis/dipole_trial%d.dat", aa);
        out_dip=fopen(file_name_dip,"w");
        
        sprintf(file_name_quad, "results-analysis/cell_quadrdupole_trial%d.dat", aa);
        out_quad=fopen(file_name_quad,"w");

        sprintf(file_name_monop, "results-analysis/monopole_trial%d.dat", aa);
        out_monop = fopen(file_name_monop,"w");

        
        
    /********* cell pos *****************/

    sprintf(file_name_centre, "cell-pos.dat"); //read positions files

    fptr = fopen(file_name_centre, "r" );
    printf("try to read file %s \n",file_name_centre);

    if(fptr==NULL)  printf("file not found\n");
    else  printf("reading file %s \n",file_name_centre);



    sizenew = 0;
    while ( fscanf( fptr, "%lf %lf ", &dummy_float, &dummy_float) != EOF )  sizenew++;
    rewind( fptr );


    cell_traject = (struct vector *) malloc( sizenew * sizeof( struct vector ) );

    cellpos =(struct vector *)malloc(sizenew *sizeof(struct vector));
    cellvel =(struct vector *)malloc(sizenew *sizeof(struct vector));

    cellanglepos =(double*)malloc(sizenew *sizeof(double));
    cellanglevel =(double*)malloc(sizenew *sizeof(double));
    cellanglevel_axis =(double*)malloc(sizenew *sizeof(double));


      if ( cell_traject == NULL ) printf( "Error alocating memory\n" ), exit( 0 );


      for ( a = 0; a < sizenew; a++ ) fscanf( fptr, "%lf %lf", &cell_traject[a].x, &cell_traject[a].y );

      fclose(fptr);


    cell_initial_pos.x = cell_traject[0].x;
    cell_initial_pos.y = cell_traject[0].y;

    double angle;

    angle = 0.0 ;
    //angle = 2.8842689386187046;

    u.x = cos(angle);
    u.y = sin(angle);

    u_perp.x = sin(angle);
    u_perp.y = -cos(angle);

    projection_cell_axis( &cell_traject[0], &cell_initial_pos, &u,  &u_perp,  cellpos,
                         cellvel, cellanglepos, cellanglevel,  sizenew);

    
    for ( a = 0; a < sizenew; a++ ){

        if( cellanglevel[a] > 0.5*M_PI) cellanglevel_axis[a] = cellanglevel[a]-M_PI;
        else if( cellanglevel[a] < -0.5*M_PI) cellanglevel_axis[a] = cellanglevel[a]+ M_PI;
        else cellanglevel_axis[a] = cellanglevel[a];

    }



    sprintf(file_name_centre, "cell-pos.dat"); //read positions files

    /********* apro file in sequenza *****************/
    fptr = fopen(file_name_centre, "r" );
    printf("try to read file %s \n",file_name_centre);

    if(fptr==NULL)  printf("file not found\n");
    else  printf("reading file %s \n",file_name_centre);

    sizenew = 0;
    while ( fscanf( fptr, "%lf %lf ", &dummy_float, &dummy_float) != EOF )  sizenew++;
    rewind( fptr );


    dipole_center = (struct vector *) malloc( sizenew * sizeof( struct vector ) );

    dist =(struct vector *)malloc(sizenew *sizeof(struct vector));


    if ( dipole_center == NULL ) printf( "Error alocating memory\n" ), exit( 0 );


    for ( a = 0; a < sizenew; a++ ) fscanf( fptr, "%lf %lf", &dipole_center[a].x, &dipole_center[a].y );

    fclose(fptr);

    for ( a = 0; a < sizenew; a++ ){

        dist[a].x =  (dipole_center[a].x-cell_initial_pos.x)*u.x + (dipole_center[a].y-cell_initial_pos.y)*u.y ;
        dist[a].y =   (dipole_center[a].x-cell_initial_pos.x)*u_perp.x + (dipole_center[a].y-cell_initial_pos.y)*u_perp.y ;

    }


    for(k=N_RIP0;k<N_RIP;k++){

        sizes[k] = 0;
        sizes_neg[k].x = 0;
        sizes_neg[k].y = 0;
        sizes_pos[k].x = 0;
        sizes_pos[k].y = 0;

    }

// begins the repetitions
  for(k=N_RIP0; k<N_RIP; k++){

      /****** file of traction vector field ******/
      sprintf(file_name_pos, "pos-%d.dat",k);
      sprintf(file_name_vel, "vel-%d.dat",k);

      /************ open position file  *************/
      fptr = fopen(file_name_pos, "r" );
      printf("try to read file %s \n",file_name_pos);

      if(fptr==NULL)  printf("file not found\n");
      else  printf("reading file positions %s \n",file_name_pos);
      while ( fscanf( fptr, "%lf %lf", &dummy_float, &dummy_float) != EOF) sizes[k]++;

      rewind( fptr );

      data_pos = (struct vector *)malloc( sizes[k] * sizeof( struct vector) );

      if( data_pos == NULL ) printf( "Error alocating memory\n" ), exit( 0 );
      printf( "realiz %d size %d \n", k, sizes[k] );

      for ( i = 0; i < sizes[k]; i++ ) fscanf( fptr, "%lf %lf", &data_pos[i].x, &data_pos[i].y);
      //att correct here
      fclose(fptr);

      /************ open velocity (deformation rate) file  *************/
      fptr = fopen(file_name_vel, "r" );
      if(fptr==NULL)  printf("file not found\n");
      else  printf("reading file %s \n",file_name_vel);

      sizec = 0;
      while ( fscanf( fptr, "%lf %lf", &dummy_float, &dummy_float) != EOF )  sizec++;

      rewind( fptr );
      data_vel = (struct vector *) malloc( sizes[k] * sizeof( struct vector ) );

      if ( data_vel == NULL ) printf( "Error alocating memory\n" ), exit( 0 );

      for ( i = 0; i < sizes[k]; i++ ) fscanf( fptr, "%lf %lf", &data_vel[i].x, &data_vel[i].y);
      //att correct here

      fclose(fptr);

      /*******************************************/
      /*                                         */
      /******************************************/


      pos =(struct vector *)malloc(sizes[k] *sizeof(struct vector));
      vel =(struct vector*)malloc(sizes[k] *sizeof(struct vector));

      pos_dip_axis =(struct vector *)malloc(sizes[k] *sizeof(struct vector));
      vel_dip_axis =(struct vector*)malloc(sizes[k] *sizeof(struct vector));

      velX =(double*)malloc(sizes[k] *sizeof(double));
      velY =(double*)malloc(sizes[k] *sizeof(double));

      anglepos =(double*)malloc(sizes[k] *sizeof(double));
      anglevel =(double*)malloc(sizes[k] *sizeof(double));

      magnpos =(double*)malloc(sizes[k] *sizeof(double));
      magnvel =(double*)malloc(sizes[k] *sizeof(double));

      for( i = 0; i <  sizes[k]; i++){
          
          vel[i].x = 0;
          vel[i].y = 0;
          
          pos[i].x = 0;
          pos[i].y = 0;
          
          pos_dip_axis[i].x = 0 ;
          pos_dip_axis[i].y = 0 ;
          
          vel_dip_axis[i].x = 0 ;
          vel_dip_axis[i].y = 0 ;
          
          magnpos[i] = 0;
          magnvel[i] = 0;
          
          anglepos[i]= 0;
          anglevel[i]= 0;
          
          velX[i] = 0 ;
          velY[i] = 0 ;
          
      }

      project_traction( &data_vel[0], &u, &u_perp, &vel[0],   sizes[k]);

      project_position( &data_pos[0], &cell_initial_pos, &u, &u_perp, &pos[0], sizes[k]);

      get_polar_coords(&vel[0], &anglevel[0], &magnvel[0] , sizes[k]);

      get_polar_coords(&pos[0], &anglepos[0], &magnpos[0] , sizes[k]);


      for( i = 0; i < sizes[k]; i++){

              if(vel[i].x >= 0.0) sizes_pos[k].x++;
              if(vel[i].x < 0.0) sizes_neg[k].x++;
              if(vel[i].y >= 0.0) sizes_pos[k].y++;
              if(vel[i].y < 0.0) sizes_neg[k].y++;

      }

      vel_neg[k] =(struct vector *)malloc(sizes[k] *sizeof(struct vector));
      vel_pos[k] =(struct vector *)malloc(sizes[k] *sizeof(struct vector));

      count_pos.x =0;
      count_pos.y =0;
      count_neg.x =0;
      count_neg.y =0;

      for( i = 0; i < sizes[k]; i++){

              if(vel[i].x >= 0.0){

                  vel_pos[k][count_pos.x].x = vel[i].x  ;
                  count_pos.x ++;
              }

              if(vel[i].x < 0.0){
                  vel_neg[k][count_neg.x].x = vel[i].x  ;
                  count_neg.x ++;
              }

              if(vel[i].y >= 0.0){
                  vel_pos[k][count_pos.y].y = vel[i].y  ;
                  count_pos.y ++;
              }

              if(vel[i].y < 0.0){
                  vel_neg[k][count_neg.y].y = vel[i].y  ;
                  count_neg.y ++;
              }
      }



      /***** media, varianza, max e min  ******/

        meanmagnitudestress = gsl_stats_mean(magnvel, 1, sizes[k]);
        variancemagnitudestress = gsl_stats_variance(magnvel, 1, sizes[k]);

        meanmagnitudepos = gsl_stats_mean(magnpos, 1, sizes[k]);
        variancemagnitudepos = gsl_stats_variance(magnpos, 1, sizes[k]);
      
      
      int countereffective = 0 ;
      double dist_x, dist_y, check_dist;
      
        for( i = 0; i < sizes[k]; i++){
            
            dist_x= pos[i].x - cell_traject[k].x ;
            dist_y= pos[i].y - cell_traject[k].y;

            check_dist = sqrt( (dist_x)*(dist_x) + (dist_y)*(dist_y) );
            
            if( (int)check_dist  < 11600 + aa*4){
                
                velX[i] = vel[i].x;
                velY[i] = vel[i].y;

                countereffective++;
                
            }
            
            else{
                
                velX[i] = 0 ;
                velY[i] = 0 ;
            }



        }

        meanvelx = ( gsl_stats_mean(velX, 1, countereffective ) ) ;
        meanvely = ( gsl_stats_mean(velY, 1, countereffective ) ) ;

        variancevelx = gsl_stats_variance(velX, 1, countereffective ) ;
        variancevely = gsl_stats_variance(velY, 1, countereffective) ;


        countereffective = 0;


      if(OLD_METHOD == 1){


          get_dipole(&data_pos[0],   cell_traject[k],  &data_vel[0],  &dipole,  screening = 0,  sizes[k] ,aa );
          get_quadrupole(&data_pos[0],   cell_traject[k],  &data_vel[0],  &dipole, &quadrupole,  screening = 0,  sizes[k], aa );


      }


      if(ALONG_CELL_AXIS == 1){

          get_dipole(&pos[0],   cellpos[k],  &vel[0],  &dipole,  screening = 0,  sizes[k] , aa);
          get_quadrupole(&pos[0],   cellpos[k],  &vel[0],  &dipole, &quadrupole,  screening = 0,  sizes[k], aa );

      }



      if(ALONG_DIPOLE_AXIS == 1){

          get_dipole(&data_pos[0],   cell_traject[k],  &data_vel[0],  &dipole,  screening = 0,  sizes[k] , aa);

          get_dipole_axis(&trace, &determ, &dipole, &max_eigenv,  &min_eigenv, &angle_min_eigenv, &angle_max_eigenv);

          struct vector u_dip, u_dip_perp;

          u_dip.x = cos(angle_max_eigenv);
          u_dip.y = sin(angle_max_eigenv);

          u_dip_perp.x = sin(angle_max_eigenv);
          u_dip_perp.y = -cos(angle_max_eigenv);

          project_traction( &data_vel[0], &u_dip, &u_dip_perp, &vel_dip_axis[0],  sizes[k]);

          project_position( &data_pos[0], &dipole_center[k], &u_dip, &u_dip_perp, &pos_dip_axis[0],  sizes[k]);

          get_dipole(&pos_dip_axis[0],   dipole_center[k],  &vel_dip_axis[0],  &dipole,  screening = 0,  sizes[k], aa);

          get_quadrupole(&pos_dip_axis[0],   dipole_center[k],  &vel_dip_axis[0],  &dipole, &quadrupole,  screening = 0,  sizes[k], aa );

      }

      get_dipole_axis(&trace, &determ, &dipole, &max_eigenv,  &min_eigenv, &angle_min_eigenv, &angle_max_eigenv);


      Omega = (0.5*(dipole.xy - dipole.yx)); // asymmetric part of the dipole tensor



      if( angle_max_eigenv > 0.5*M_PI){

          axis_max_eigenv = angle_max_eigenv-M_PI;
        }

      else if( angle_max_eigenv < -0.5*M_PI){

          axis_max_eigenv = angle_max_eigenv+ M_PI;

      }

      else axis_max_eigenv = angle_max_eigenv;


      printf( "realiz %d countereff %d \n", k, countereffective);

      if(aa==1){

      fprintf(out_omega,"%d\t %lf\n",k, Omega );

      fprintf(out_max_eigenv,"%d\t %lf\t %lf\t %d\n",k, max_eigenv, min_eigenv, sizes[k]);

      fprintf(out_angle_eigenv,"%d\t %lf\t %lf\n",k, angle_max_eigenv, angle_min_eigenv);
      fprintf(out_angle_eigenv_axis,"%d\t %lf\n",k, axis_max_eigenv);

      fprintf(out_angle_cell_speed,"%d\t %lf\t %lf\n", k, cellanglevel[k], cos(cellanglevel_axis[k]-axis_max_eigenv) );

      fprintf(out_angle_cell_speed_axis,"%d\t  %lf\n", k, sqrt(cellvel[k].x*cellvel[k].x + cellvel[k].y*cellvel[k].y)*cellanglevel_axis[k]);

      fprintf(out_speed_cell,"%d\t %lf\t %lf\t %lf\n", k, cellvel[k].x, cellvel[k].y, cellvel[k].x*cellvel[k].x + cellvel[k].y*cellvel[k].y);

      fprintf(out_cell_pos,"%lf\t %lf\n", cellpos[k].x, cellpos[k].y);
          
          
      fprintf(out_tract,"%d\t %lf\t %lf\t %lf\t %lf\n",k, meanmagnitudestress,variancemagnitudestress , meanmagnitudepos, variancemagnitudepos );
          

      }

      fprintf(out_quad,"%d\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\n",k, quadrupole.xxx, quadrupole.xxy, quadrupole.xyx, quadrupole.xyy,
              quadrupole.yxx, quadrupole.yxy, quadrupole.yyx, quadrupole.yyy);
      
      fprintf(out_dip,"%d\t %lf\t  %lf\t %lf\t %lf\t %lf\t %lf\n",k, dipole.xx, dipole.xy , dipole.yx , dipole.yy, trace, determ);

      fprintf(out_monop,"%d\t %lf\t %lf\t %lf\t %lf\n", k, meanvelx, variancevelx,  meanvely, variancevely);


      
      free(data_pos);
      free(data_vel);

      free(pos);
      free(vel);

      free(anglepos);
      free(anglevel);

      free(magnpos);
      free(magnvel);



    } // k --- N_RIP2

        
        fclose(out_quad);
        fclose(out_dip);
        fclose(out_monop);

        if(aa > 1){
            fclose(out_cell_pos);

        }

    }


  return(0);
}



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
//
//      FUNCTIONS
//
//
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  project on cell axis
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void projection_cell_axis(struct vector * cell_traject, struct vector * cell_initial_pos, struct vector * u, struct vector * u_perp,   struct vector * cellpos, struct vector * cellvel, double * cellanglepos, double * cellanglevel,  int sizenew){

    int a;

for ( a = 0; a < sizenew; a++ ){

    cellpos[a].x =  (cell_traject[a].x-cell_initial_pos->x)*u->x + (cell_traject[a].y-cell_initial_pos->y)*u->y ;
    cellpos[a].y =   (cell_traject[a].x-cell_initial_pos->x)*u_perp->x + (cell_traject[a].y-cell_initial_pos->y)*u_perp->y ;

    cellanglepos[a]= atan2(cellpos[a].y, cellpos[a].x);

    if( a <= sizenew-2){

        cellvel[a].x = (cell_traject[a+1].x - cell_traject[a].x)*u->x  + (cell_traject[a+1].y - cell_traject[a].y)*u->y ;
        cellvel[a].y = (cell_traject[a+1].y - cell_traject[a].y)*u_perp->y  + (cell_traject[a+1].x - cell_traject[a].x)*u_perp->x   ;

        cellanglevel[a]= atan2(cellvel[a].y, cellvel[a].x);

    }

    else{

        cellvel[a].x = 0.0 ;
        cellvel[a].y = 0.0 ;
        cellanglevel[a] = 0.0;
    }

}


}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  polar coords
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void get_polar_coords(struct vector * v,  double * angle,  double * magnitude, int siz){

    int i;

    for( i = 0; i < siz; i++){

    angle[i] = atan2(v[i].y, v[i].x);
    magnitude[i] = sqrt(v[i].x*v[i].x + v[i].y*v[i].y);

    }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  project traction on a given axis
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



void project_traction(struct vector * data_vel, struct vector * u, struct vector * u_perp, struct vector * vel,
                                  int siz){

    int i;

    for( i = 0; i < siz; i++){

        vel[i].x =  data_vel[i].x*u->x ; //+ data_vel[i].y*u->y;
        vel[i].y =  data_vel[i].y*u_perp->y ; // + data_vel[i].x*u_perp->x ;
        //printf("%lf \n", u->x);

    }

}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  project position on a given axis
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


void project_position(struct vector * data_pos, struct vector * new_origin, struct vector * u, struct vector * u_perp, struct vector * pos,  int siz){

    int i;

    for( i = 0; i < siz; i++){

        pos[i].y = (data_pos[i].x)*u->x ; //+ (data_pos[i].x-new_origin->y)*u->y; // att this is correct
        pos[i].x = (data_pos[i].y)*u_perp->y ; // + (data_pos[i].y-new_origin->x)*u_perp->x ; // att this is correct
        //    pos[i].y = data_pos[i].x ; // att this is correct
        //    pos[i].x = data_pos[i].y ; // att this is correct


    }

}



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//     Dipole
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void get_dipole(struct vector * pos, struct vector cellposs, struct vector * vel, struct dipole_term * dipole, int screening, int siz, int aa ){


    int i , countereffective;
    double screen, check_dist, dipolemix, dist_x, dist_y;

    dipole->xx = 0 ;
    dipole->xy = 0 ;
    dipole->yx = 0 ;
    dipole->yy = 0 ;
    dipolemix = 0 ;

    countereffective = 0 ;
    
    
    for( i = 0; i < siz; i++){

        dist_x = pos[i].x - cellposs.x ;
        dist_y = pos[i].y - cellposs.y ;

        if( screening == 1 ){

            screen = (dist_x)*(dist_x) + (dist_y)*(dist_y);

            if(isnan(screen)==1) printf("problem with the screening computation %lf\n", screen ), exit(0);

        }


        else screen = 1. ;

        check_dist = sqrt((dist_x)*(dist_x) + (dist_y)*(dist_y) );

        if( (int)check_dist  < 600 + aa*4){
            
            countereffective++;
            
        }
        
        else{
            
            dist_x=0;
            dist_y=0;
        }


        dipole->xx += (dist_x*vel[i].x)/screen;
        dipole->xy += (dist_x*vel[i].y)/screen;
        dipole->yx += (dist_y*vel[i].x)/screen;
        dipole->yy += (dist_y*vel[i].y)/screen;

        }

    if (countereffective > 0) {

        dipole->xx /= (double)countereffective;
        dipole->yy /= (double)countereffective;

        dipolemix = (0.5*(dipole->xy + dipole->yx))/(double)countereffective;   //  symmetric part of the dipole tensor

        dipole->xy = dipolemix;
        dipole->yx = dipolemix;

    }

}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//     Quadrupole
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void get_quadrupole(struct vector * pos, struct vector cellposs, struct vector * vel,
                                    struct dipole_term * dipole, struct quadrupole_term * quadrupole, int screening, int siz , int aa){

    int i, countereffective;
;
    double screen, check_dist, screening_quad, dist_x, dist_y;

    quadrupole->xxx = 0 ;
    quadrupole->xxy = 0 ;
    quadrupole->xyx = 0 ;
    quadrupole->xyy= 0 ;

    quadrupole->yxx = 0 ;
    quadrupole->yxy = 0 ;
    quadrupole->yyx = 0 ;
    quadrupole->yyy = 0 ;
    
    countereffective = 0 ;

    for( i = 0; i < siz; i++){

        dist_x = pos[i].x - cellposs.x;
        dist_y = pos[i].y - cellposs.y;
        
        if(screening==1) {

            screen = (dist_x)*(dist_x) + (dist_y)*(dist_y);
            if(isnan(screen)==1) printf("problem %lf\n", screen ), exit(0);

        }

        else screen = 1. ;

        check_dist = sqrt((dist_x)*(dist_x) + (dist_y)*(dist_y));

        if( (int)check_dist  < 600 + aa*4){
            
            countereffective++;
            
        }
        
        else{
            
            dist_x=0;
            dist_y=0;
        }

        screening_quad = screen*screen;

        quadrupole->xxx += (dist_x*dist_x*vel[i].x)/screening_quad;
        quadrupole->xxy += (dist_x*dist_x*vel[i].y)/screening_quad;
        quadrupole->xyx += (dist_x*dist_y*vel[i].x)/screening_quad;
        quadrupole->xyy += (dist_x*dist_y*vel[i].y)/screening_quad;

        quadrupole->yxx += (dist_y*dist_x*vel[i].x)/screening_quad;
        quadrupole->yxy += (dist_y*dist_x*vel[i].y)/screening_quad;
        quadrupole->yyx += (dist_y*dist_y*vel[i].x)/screening_quad;
        quadrupole->yyy += (dist_y*dist_y*vel[i].y)/screening_quad;

        if(SUBTRACT_DIPOLE == 1){

            quadrupole->xxx  -= dist_x*dipole->xx;
            quadrupole->xxy  -= dist_x*dipole->xy;
            quadrupole->xyx  -= dist_x*dipole->yx;
            quadrupole->xyy  -= dist_x*dipole->yy;
            quadrupole->yxx  -= dist_y*dipole->xx;
            quadrupole->yxy  -= dist_y*dipole->xy;
            quadrupole->yyx  -= dist_y*dipole->yx;
            quadrupole->yyy  -= dist_x*dipole->yy;

        }

    }

    if (countereffective > 0) {


        quadrupole->xxx  /= (double)countereffective;
        quadrupole->xxy  /= (double)countereffective;
        quadrupole->xyx  /= (double)countereffective;
        quadrupole->xyy  /= (double)countereffective;

        quadrupole->yxx  /= (double)countereffective;
        quadrupole->yxy  /= (double)countereffective;
        quadrupole->yyx   /= (double)countereffective;
        quadrupole->yyy   /= (double)countereffective;



    }
 }




//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//     Get dipole axis
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void get_dipole_axis(double * trace, double * determ, struct dipole_term * dipole, double * max_eigenv, double * min_eigenv,
                      double * angle_min_eigenv,
                     double * angle_max_eigenv){

    double memo_plus, memo_minus, u_max_x, u_max_y,  u_min_x, u_min_y;

    *trace = dipole->xx + dipole->yy;
    *determ = dipole->xx*dipole->yy - dipole->xy*dipole->yx;

    memo_plus = *trace + sqrt(*trace*(*trace)-4.0*(*determ) );
    memo_minus = *trace - sqrt(*trace*(*trace)-4.0*(*determ) );

    if ( fabs(memo_plus) >  fabs(memo_minus) ) {

        *max_eigenv = 0.5*(memo_plus);
        *min_eigenv = 0.5*(memo_minus);
    }

    else{

        *max_eigenv = 0.5*(memo_minus);
        *min_eigenv = 0.5*(memo_plus);
    }



    if( dipole->xy != 0 ){

    //mod_u_max = sqrt(1+ ((max_eigenv- dipole.xx)*(max_eigenv- dipole.xx))/(dipole.xy*dipole.xy));
    //mod_u_min = sqrt(1+ ((min_eigenv- dipole.xx)*(min_eigenv- dipole.xx))/(dipole.xy*dipole.xy));

    u_max_x = 1.;
    u_max_y = (*max_eigenv - dipole->xx)/(dipole->xy);

    u_min_x = 1.;
    u_min_y =(*min_eigenv- dipole->xx)/(dipole->xy);

    }

    else{
        u_max_x = 1.;
        u_max_y = 0.;

        u_min_x = 0.;
        u_min_y = 1. ;

    }

    *angle_min_eigenv = atan2(u_min_y , u_min_x);
    *angle_max_eigenv  = atan2(u_max_y , u_max_x);


}

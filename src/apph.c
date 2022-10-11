//=============================================================================
// in cartesian coordinates t,x,y,z for x,y,z in [-1,1]
//
// NOTES: ==> uses the global variables of AdS4D
//        ==> hierarchy between Lmin and Lmax must be in-sync
//        ==> assumed time structure ... n,nm1,np1
//        ==> the current AH structure is specified by c_AH (0 indexed)
//
// The following routine searches for an apparent horizon over levels
// AH_Lmin to AH_Lmax.
//
// The hypersurface describing the AH is given by
//
// r = AH_R(chi,phi)
//
// An iterative flow-method is used, to absolute tolerance AH_tol, maximum
// iteration AH_max_iter. 
//
// If AH_tol<0 then the tolerance is dynamically adjusted (via change in R of min(dx)*-AH_tol)
// 
// If (use_R_ic), then existing value of R is used as an initial guess, else
// AH_rsteps values of R=AH_r0 to R=AH_r1 are used, with the one giving the
// smallest expansion used to start.
//
// AH_theta(AH_Nchi,AH_Nphi) is a real work array used to hold the null expansion.
// AH_own(AH_Nchi,AH_Nphi) and AH_lev(AH_Nchi,AH_Nphi) are integer work arrays describing the
// ownership of points.
// 
// returns (true) if found, and if so M & J are filled in to reflect the
// effect mass and angular momentum via the dynamical horizons approach i
// (not implemented yet)
//=============================================================================

#include <stdlib.h>
#include <stdio.h>
#include <pamr.h>
#include <amrd.h>
#include <math.h>
#include <m_util_r8.h>
#include <bbhutil.h>
#include <mpi.h>
#include "AdS4D.h"
#include "apph.h"

#define UNOWNED -1e10
real AH_ct[MAX_BHS];

//=============================================================================
// some utility routines called by find_apph below
//=============================================================================

//-----------------------------------------------------------------------------
// initialize ownership array ... i.e., those points this node is
// responsible for
//-----------------------------------------------------------------------------
int fill_own(int Lmax, int ltrace, int *first)
{
   int i,np,L,valid;

   np=AH_Nchi[c_AH]*AH_Nphi[c_AH];

   for (i=0; i<np; i++) AH_own[c_AH][i]=-1;

   for (L=Lmax; L>=AH_Lmin[c_AH]; L--)
   {
      valid=PAMR_init_s_iter(L,PAMR_AMRH,1);  // loop over ALL grids!

      while(valid)
      {
         ldptr();

         ah_fill_own_(AH_R[c_AH],AH_xc[c_AH],AH_own[c_AH],AH_lev[c_AH],bbox,&dx,&dy,&dz,&g_rank,&L,&AH_Nchi[c_AH],&AH_Nphi[c_AH], &axisym);
         valid=PAMR_next_g();
      }
   }
   *first=0;

   for (i=0; i<np ; i++) if (AH_own[c_AH][i]==-1) 
   {
      if (my_rank==0 && ltrace)
      printf("find_apph: error ... point %i,%i is not owned; at AH_R=%f\n,chi=%f\n",i%AH_Nchi[c_AH]+1,i/AH_Nchi[c_AH]+1,AH_R[c_AH][i],(i%AH_Nchi[c_AH])*M_PI/(AH_Nchi[c_AH]-1));
      return 0;        
   }

   return 1;
}

//-----------------------------------------------------------------------------
// computes theta and AH metric components, assuming ownership calculated ... alters AH_w1,AH_w3,AH_w4
//
// is_ex set to 1 if any point couldn't be calculated due to closeness
// of excision zone ... returns theta average
//-----------------------------------------------------------------------------
#define USE_SMOOTH_A 0 //HB//
#define MAX_TRACE 5000
real fill_theta_ahmetric(double *AH_theta0, real eps0, real *area, real *c_equat, real *c_polar, real *c_polar2, int *is_ex,
                        real *ief_bh_r0,real *a_rot0, int *kerrads_background, real *ct)
{
   int i,j,np,valid,dvtrace=0,i0,j0,is_int;
   static int num_trace=0;
   int size_name;
   char *name;
   int AH_shape[3],rank;
   real AH_bbox[6],resid,da[4],area_owned[4],area_global[4];  // elements 2,3,4 for circumferences.

   int k;
   int quit_cleanly=0;
   real tmp=1.0;

   // outputs AH_*_iter gfns
   dvtrace=1;

   np=AH_Nchi[c_AH]*AH_Nphi[c_AH];

   area_owned[0]=area_owned[1]=area_owned[2]=area_owned[3]=0;
   *is_ex=0;

   for (i=0; i<np; i++) 
   {
      AH_theta0[i]=UNOWNED;
      AH_x0[c_AH][i]=0;
      AH_y0[c_AH][i]=0;
      AH_z0[c_AH][i]=0;
      AH_g0_xx[c_AH][i]=0;
      AH_g0_xy[c_AH][i]=0;
      AH_g0_xz[c_AH][i]=0;
      AH_g0_yy[c_AH][i]=0;
      AH_g0_yz[c_AH][i]=0;
      AH_g0_zz[c_AH][i]=0;
      AH_g0_chichi[c_AH][i]=0;
      AH_g0_chiphi[c_AH][i]=0;
      AH_g0_phiphi[c_AH][i]=0;
      AH_kretsch[c_AH][i]=0;
      AH_riemanncube[c_AH][i]=0;
      AH_ahr[c_AH][i]=0;
      AH_dch[c_AH][i]=0;
      AH_dph[c_AH][i]=0;
      AH_da0[c_AH][i]=0;
      AH_dcq[c_AH][i]=0;
      AH_dcp[c_AH][i]=0;
      AH_dcp2[c_AH][i]=0;
      if (AH_own[c_AH][i]==my_rank)
      {
         i0=i%AH_Nchi[c_AH]+1;
         j0=i/AH_Nchi[c_AH]+1;
         valid=PAMR_init_s_iter(AH_lev[c_AH][i],PAMR_AMRH,0); 
         while(valid)
         {
            ldptr_bbox();
            ah_is_int_(&is_int,AH_R[c_AH],AH_xc[c_AH],&i0,&j0,bbox,&dx,&dy,&dz,&AH_Nchi[c_AH],&AH_Nphi[c_AH],&axisym);
            if (is_int)
            {
               ldptr(); 

               // compute full theta
               if ((*ct)!=0)
               {
                  //(NOTE: for t>t0, have cycled time sequence np1,n,nm1 to time sequence n,nm1,np1,
                  // so here, time level n is the most advanced time level)   
                  calc_exp_metric0_(AH_R[c_AH],AH_xc[c_AH],AH_theta0,
                         &i0,&j0,&AH_Nchi[c_AH],
                         &AH_Nphi[c_AH],theta,f,&da[0],&da[1],&da[2],&da[3],
                         AH_x0[c_AH],AH_y0[c_AH],AH_z0[c_AH],
                         AH_g0_xx[c_AH],AH_g0_xy[c_AH],AH_g0_xz[c_AH],
                         AH_g0_yy[c_AH],AH_g0_yz[c_AH],AH_g0_zz[c_AH],
                         AH_g0_chichi[c_AH],AH_g0_chiphi[c_AH],AH_g0_phiphi[c_AH],
                         AH_kretsch[c_AH],
                         AH_riemanncube[c_AH],
                         AH_ahr[c_AH],AH_dch[c_AH],AH_dph[c_AH],
                         AH_da0[c_AH],AH_dcq[c_AH],AH_dcp[c_AH],AH_dcp2[c_AH],
                         gb_tt_n,gb_tt_nm1,gb_tt_np1,
                         gb_tx_n,gb_tx_nm1,gb_tx_np1,
                         gb_ty_n,gb_ty_nm1,gb_ty_np1,
                         gb_tz_n,gb_tz_nm1,gb_tz_np1,
                         gb_xx_n,gb_xx_nm1,gb_xx_np1,
                         gb_xy_n,gb_xy_nm1,gb_xy_np1,
                         gb_xz_n,gb_xz_nm1,gb_xz_np1,
                         gb_yy_n,gb_yy_nm1,gb_yy_np1,
                         gb_yz_n,gb_yz_nm1,gb_yz_np1,
                         gb_zz_n,gb_zz_nm1,gb_zz_np1,
                         Hb_t_n,Hb_t_nm1,Hb_t_np1,
                         Hb_x_n,Hb_x_nm1,Hb_x_np1,
                         Hb_y_n,Hb_y_nm1,Hb_y_np1,
                         Hb_z_n,Hb_z_nm1,Hb_z_np1,
                         phi1_n,phi1_nm1,phi1_np1,
                         kretsch_n,
                         riemanncube_n,
                         &AdS_L,x,y,z,&dt,chr,&AMRD_ex,&AMRD_do_ex,&Nx,&Ny,&Nz,&axisym,
                         ief_bh_r0,a_rot0,kerrads_background);
               }
               else
               {
                  //(NOTE: for t=t0, have *not* cycled time sequence, so still np1,n,nm1,
                  // so here, time level np1 is the most advanced time level) 
                  calc_exp_metric0_(AH_R[c_AH],AH_xc[c_AH],AH_theta0,
                         &i0,&j0,&AH_Nchi[c_AH],
                         &AH_Nphi[c_AH],theta,f,&da[0],&da[1],&da[2],&da[3],
                         AH_x0[c_AH],AH_y0[c_AH],AH_z0[c_AH],
                         AH_g0_xx[c_AH],AH_g0_xy[c_AH],AH_g0_xz[c_AH],
                         AH_g0_yy[c_AH],AH_g0_yz[c_AH],AH_g0_zz[c_AH],
                         AH_g0_chichi[c_AH],AH_g0_chiphi[c_AH],AH_g0_phiphi[c_AH],
                         AH_kretsch[c_AH],
                         AH_riemanncube[c_AH],
                         AH_ahr[c_AH],AH_dch[c_AH],AH_dph[c_AH],
                         AH_da0[c_AH],AH_dcq[c_AH],AH_dcp[c_AH],AH_dcp2[c_AH],
                         gb_tt_np1,gb_tt_n,gb_tt_nm1,
                         gb_tx_np1,gb_tx_n,gb_tx_nm1,
                         gb_ty_np1,gb_ty_n,gb_ty_nm1,
                         gb_tz_np1,gb_tz_n,gb_tz_nm1,
                         gb_xx_np1,gb_xx_n,gb_xx_nm1,
                         gb_xy_np1,gb_xy_n,gb_xy_nm1,
                         gb_xz_np1,gb_xz_n,gb_xz_nm1,
                         gb_yy_np1,gb_yy_n,gb_yy_nm1,
                         gb_yz_np1,gb_yz_n,gb_yz_nm1,
                         gb_zz_np1,gb_zz_n,gb_zz_nm1,
                         Hb_t_np1,Hb_t_n,Hb_t_nm1,
                         Hb_x_np1,Hb_x_n,Hb_x_nm1,
                         Hb_y_np1,Hb_y_n,Hb_y_nm1,
                         Hb_z_np1,Hb_z_n,Hb_z_nm1,
                         phi1_np1,phi1_n,phi1_nm1,
                         kretsch_np1,
                         riemanncube_np1,
                         &AdS_L,x,y,z,&dt,chr,&AMRD_ex,&AMRD_do_ex,&Nx,&Ny,&Nz,&axisym,
                         ief_bh_r0,a_rot0,kerrads_background);
               }

               area_owned[0]+=da[0];
               area_owned[1]+=da[1];
               area_owned[2]+=da[2];
               area_owned[3]+=da[3];
               valid=0;
            }
            else valid=PAMR_next_g();
         }
         // to quit "cleanly"
         if (AH_theta0[i]==UNOWNED) quit_cleanly=1;

      }

   }

   if (quit_cleanly==1) AMRD_stop("fill_theta_ahmetric: error ... point 'unowned'",""); //quitting cleanly (cf. above)

   // for each i point on the AH surface, save max{AH_theta0[i]}_allprocessors into AH_w1[i],
   MPI_Allreduce(AH_theta0,AH_w1[c_AH],np,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
   // copy AH_w1 into AH_theta0
   for (i=0; i<np; i++) {AH_theta0[i]=AH_w1[c_AH][i];}

   // and sum{area_owned[i]} into area_global[i]
   MPI_Allreduce(area_owned,area_global,4,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

   *area    =area_global[0];
   *c_equat =area_global[1];
   *c_polar =area_global[2];
   *c_polar2=area_global[3];

//   if (my_rank==0)
//   {
//      printf("area=%lf,c_equat=%lf,c_polar=%lf,c_polar2=%lf\n",*area,*c_equat,*c_polar,*c_polar2);
//   }

   if (*area<0) { *is_ex=1; *area=0; }

   if (output_moreAHquant_sdf||output_moreAHquant_ascii)
   {
      if (output_moreAHquant_ascii)
      { 
        // save cartesian coordinates on the horizon, accounting for all processors
        MPI_Allreduce(AH_x0[c_AH],AH_w3[c_AH],np,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
        MPI_Allreduce(AH_x0[c_AH],AH_w4[c_AH],np,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
        for (i=0; i<np; i++) {AH_x0[c_AH][i]=AH_w3[c_AH][i]+AH_w4[c_AH][i];}
        MPI_Allreduce(AH_y0[c_AH],AH_w3[c_AH],np,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
        MPI_Allreduce(AH_y0[c_AH],AH_w4[c_AH],np,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
        for (i=0; i<np; i++) {AH_y0[c_AH][i]=AH_w3[c_AH][i]+AH_w4[c_AH][i];}
        MPI_Allreduce(AH_z0[c_AH],AH_w3[c_AH],np,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
        MPI_Allreduce(AH_z0[c_AH],AH_w4[c_AH],np,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
        for (i=0; i<np; i++) {AH_z0[c_AH][i]=AH_w3[c_AH][i]+AH_w4[c_AH][i];} 
      }

      if (output_metricAH_cart_sdf||output_metricAH_cart_ascii)
      {
        // for each i point on the AH surface, save sum{AH_g0_xx[c_AH][i]}_allprocessors,etc.
        MPI_Allreduce(AH_g0_xx[c_AH],AH_w3[c_AH],np,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
        MPI_Allreduce(AH_g0_xx[c_AH],AH_w4[c_AH],np,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
        for (i=0; i<np; i++) {AH_g0_xx[c_AH][i]=AH_w3[c_AH][i]+AH_w4[c_AH][i];}
        MPI_Allreduce(AH_g0_xy[c_AH],AH_w3[c_AH],np,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
        MPI_Allreduce(AH_g0_xy[c_AH],AH_w4[c_AH],np,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
        for (i=0; i<np; i++) {AH_g0_xy[c_AH][i]=AH_w3[c_AH][i]+AH_w4[c_AH][i];}
        MPI_Allreduce(AH_g0_xz[c_AH],AH_w3[c_AH],np,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
        MPI_Allreduce(AH_g0_xz[c_AH],AH_w4[c_AH],np,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
        for (i=0; i<np; i++) {AH_g0_xz[c_AH][i]=AH_w3[c_AH][i]+AH_w4[c_AH][i];}
        MPI_Allreduce(AH_g0_yy[c_AH],AH_w3[c_AH],np,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
        MPI_Allreduce(AH_g0_yy[c_AH],AH_w4[c_AH],np,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
        for (i=0; i<np; i++) {AH_g0_yy[c_AH][i]=AH_w3[c_AH][i]+AH_w4[c_AH][i];}
        MPI_Allreduce(AH_g0_yz[c_AH],AH_w3[c_AH],np,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
        MPI_Allreduce(AH_g0_yz[c_AH],AH_w4[c_AH],np,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
        for (i=0; i<np; i++) {AH_g0_yz[c_AH][i]=AH_w3[c_AH][i]+AH_w4[c_AH][i];}
        MPI_Allreduce(AH_g0_zz[c_AH],AH_w3[c_AH],np,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
        MPI_Allreduce(AH_g0_zz[c_AH],AH_w4[c_AH],np,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
        for (i=0; i<np; i++) {AH_g0_zz[c_AH][i]=AH_w3[c_AH][i]+AH_w4[c_AH][i];}
      }
      
      if (output_metricAH_sph_sdf||output_metricAH_sph_ascii)
      {
        // save components of induced metric on the horizon, accounting for all processors
        MPI_Allreduce(AH_g0_chichi[c_AH],AH_w3[c_AH],np,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
        MPI_Allreduce(AH_g0_chichi[c_AH],AH_w4[c_AH],np,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
        for (i=0; i<np; i++) {AH_g0_chichi[c_AH][i]=AH_w3[c_AH][i]+AH_w4[c_AH][i];}
        MPI_Allreduce(AH_g0_chiphi[c_AH],AH_w3[c_AH],np,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
        MPI_Allreduce(AH_g0_chiphi[c_AH],AH_w4[c_AH],np,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
        for (i=0; i<np; i++) {AH_g0_chiphi[c_AH][i]=AH_w3[c_AH][i]+AH_w4[c_AH][i];}
        MPI_Allreduce(AH_g0_phiphi[c_AH],AH_w3[c_AH],np,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
        MPI_Allreduce(AH_g0_phiphi[c_AH],AH_w4[c_AH],np,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
        for (i=0; i<np; i++) {AH_g0_phiphi[c_AH][i]=AH_w3[c_AH][i]+AH_w4[c_AH][i];}
      }

      if (output_kretschAH_sdf||output_kretschAH_ascii)
      {
        // save Kretschmann scalar on the horizon, accounting for all processors
        MPI_Allreduce(AH_kretsch[c_AH],AH_w3[c_AH],np,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
        MPI_Allreduce(AH_kretsch[c_AH],AH_w4[c_AH],np,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
        for (i=0; i<np; i++) {AH_kretsch[c_AH][i]=AH_w3[c_AH][i]+AH_w4[c_AH][i];}
      }

      if (output_riemanncubeAH_sdf||output_riemanncubeAH_ascii)
      {
        // save relative Kretschmann scalar on the horizon, accounting for all processors
        MPI_Allreduce(AH_riemanncube[c_AH],AH_w3[c_AH],np,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
        MPI_Allreduce(AH_riemanncube[c_AH],AH_w4[c_AH],np,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
        for (i=0; i<np; i++) {AH_riemanncube[c_AH][i]=AH_w3[c_AH][i]+AH_w4[c_AH][i];}
      }

      if (output_diagnosticAH_ascii)
      { 
        // save components of diagnostics
        MPI_Allreduce(AH_ahr[c_AH],AH_w3[c_AH],np,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
        MPI_Allreduce(AH_ahr[c_AH],AH_w4[c_AH],np,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
        for (i=0; i<np; i++) {AH_ahr[c_AH][i]=AH_w3[c_AH][i]+AH_w4[c_AH][i];}
        MPI_Allreduce(AH_dch[c_AH],AH_w3[c_AH],np,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
        MPI_Allreduce(AH_dch[c_AH],AH_w4[c_AH],np,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
        for (i=0; i<np; i++) {AH_dch[c_AH][i]=AH_w3[c_AH][i]+AH_w4[c_AH][i];}
        MPI_Allreduce(AH_dph[c_AH],AH_w3[c_AH],np,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
        MPI_Allreduce(AH_dph[c_AH],AH_w4[c_AH],np,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
        for (i=0; i<np; i++) {AH_dph[c_AH][i]=AH_w3[c_AH][i]+AH_w4[c_AH][i];}
        MPI_Allreduce(AH_da0[c_AH],AH_w3[c_AH],np,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
        MPI_Allreduce(AH_da0[c_AH],AH_w4[c_AH],np,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
        for (i=0; i<np; i++) {AH_da0[c_AH][i]=AH_w3[c_AH][i]+AH_w4[c_AH][i];}
        MPI_Allreduce(AH_dcq[c_AH],AH_w3[c_AH],np,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
        MPI_Allreduce(AH_dcq[c_AH],AH_w4[c_AH],np,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
        for (i=0; i<np; i++) {AH_dcq[c_AH][i]=AH_w3[c_AH][i]+AH_w4[c_AH][i];}
        MPI_Allreduce(AH_dcp[c_AH],AH_w3[c_AH],np,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
        MPI_Allreduce(AH_dcp[c_AH],AH_w4[c_AH],np,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
        for (i=0; i<np; i++) {AH_dcp[c_AH][i]=AH_w3[c_AH][i]+AH_w4[c_AH][i];}
        MPI_Allreduce(AH_dcp2[c_AH],AH_w3[c_AH],np,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
        MPI_Allreduce(AH_dcp2[c_AH],AH_w4[c_AH],np,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
        for (i=0; i<np; i++) {AH_dcp2[c_AH][i]=AH_w3[c_AH][i]+AH_w4[c_AH][i];} 
      }
    }



   //--------------------------------------------------------------------------
   //HB//
   // regularity is essential at the chi=0,PI poles, and we need to do it before
   // smoothing as calc_exp_metric0 does not fill in the axis, and afterwards
   // again to make sure that theta and hence R is always exactly regular
   //--------------------------------------------------------------------------
  reg_ah_r_(AH_theta0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
  if (output_moreAHquant_sdf||output_moreAHquant_ascii)
  {
    if (output_metricAH_cart_sdf||output_metricAH_cart_ascii)
    {
     reg_ah_r_(AH_g0_xx[c_AH],&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
     reg_ah_r_(AH_g0_xy[c_AH],&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
     reg_ah_r_(AH_g0_xz[c_AH],&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
     reg_ah_r_(AH_g0_yy[c_AH],&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
     reg_ah_r_(AH_g0_yz[c_AH],&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
     reg_ah_r_(AH_g0_zz[c_AH],&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
    }
    if (output_metricAH_sph_sdf||output_metricAH_sph_ascii)
    {
     reg_ah_r_(AH_g0_chichi[c_AH],&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
     reg_ah_r_(AH_g0_chiphi[c_AH],&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
     reg_ah_r_(AH_g0_phiphi[c_AH],&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
    }
    if (output_kretschAH_sdf||output_kretschAH_ascii)
    {
      reg_ah_r_(AH_kretsch[c_AH],&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
    }
    if (output_riemanncubeAH_sdf||output_riemanncubeAH_ascii)
    {
      reg_ah_r_(AH_riemanncube[c_AH],&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
    }
  }

   // AH smoothing
   if (eps0>0 && eps0<1)
   {
      if (USE_SMOOTH_A) 
      {
        smooth_ah_r_(AH_theta0,AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
        if (output_moreAHquant_sdf||output_moreAHquant_ascii)
        {
          if (output_metricAH_cart_sdf||output_metricAH_cart_ascii)
          {
            smooth_ah_r_(AH_g0_xx[c_AH],AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
            smooth_ah_r_(AH_g0_xy[c_AH],AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
            smooth_ah_r_(AH_g0_xz[c_AH],AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
            smooth_ah_r_(AH_g0_yy[c_AH],AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
            smooth_ah_r_(AH_g0_yz[c_AH],AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
            smooth_ah_r_(AH_g0_zz[c_AH],AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
          }
          if (output_metricAH_sph_sdf||output_metricAH_sph_ascii)
          {
            smooth_ah_r_(AH_g0_chichi[c_AH],AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
            smooth_ah_r_(AH_g0_chiphi[c_AH],AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
            smooth_ah_r_(AH_g0_phiphi[c_AH],AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
          }
          if (output_kretschAH_sdf||output_kretschAH_ascii)
          {
            smooth_ah_r_(AH_kretsch[c_AH],AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
          }
          if (output_riemanncubeAH_sdf||output_riemanncubeAH_ascii)
          {
            smooth_ah_r_(AH_riemanncube[c_AH],AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
          }
        }
      }
      else 
      {
        smooth_ah_r_b_(AH_theta0,AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
        if (output_moreAHquant_sdf||output_moreAHquant_ascii)
        {
          if (output_metricAH_cart_sdf||output_metricAH_cart_ascii)
          {
            smooth_ah_r_b_(AH_g0_xx[c_AH],AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
            smooth_ah_r_b_(AH_g0_xy[c_AH],AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
            smooth_ah_r_b_(AH_g0_xz[c_AH],AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
            smooth_ah_r_b_(AH_g0_yy[c_AH],AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
            smooth_ah_r_b_(AH_g0_yz[c_AH],AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
            smooth_ah_r_b_(AH_g0_zz[c_AH],AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
          }
          if (output_metricAH_sph_sdf||output_metricAH_sph_ascii)
          {
            smooth_ah_r_b_(AH_g0_chichi[c_AH],AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
            smooth_ah_r_b_(AH_g0_chiphi[c_AH],AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
            smooth_ah_r_b_(AH_g0_phiphi[c_AH],AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
          }
          if (output_kretschAH_sdf||output_kretschAH_ascii)
          {
            smooth_ah_r_b_(AH_kretsch[c_AH],AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
          }
          if (output_riemanncubeAH_sdf||output_riemanncubeAH_ascii)
          {
            smooth_ah_r_b_(AH_riemanncube[c_AH],AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
          }
        }
      }
   }
   else if (eps0>1 && eps0<2) 
   {
      eps0-=1;
      if (USE_SMOOTH_A) 
      {
        smooth_ah_r_(AH_theta0,AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
        if (output_moreAHquant_sdf||output_moreAHquant_ascii)
        {
          if (output_metricAH_cart_sdf||output_metricAH_cart_ascii)
          {
            smooth_ah_r_(AH_g0_xx[c_AH],AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
            smooth_ah_r_(AH_g0_xy[c_AH],AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
            smooth_ah_r_(AH_g0_xz[c_AH],AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
            smooth_ah_r_(AH_g0_yy[c_AH],AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
            smooth_ah_r_(AH_g0_yz[c_AH],AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
            smooth_ah_r_(AH_g0_zz[c_AH],AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
          }
          if (output_metricAH_sph_sdf||output_metricAH_sph_ascii)
          {
            smooth_ah_r_(AH_g0_chichi[c_AH],AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
            smooth_ah_r_(AH_g0_chiphi[c_AH],AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
            smooth_ah_r_(AH_g0_phiphi[c_AH],AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
          }
          if (output_kretschAH_sdf||output_kretschAH_ascii)
          {
            smooth_ah_r_(AH_kretsch[c_AH],AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
          }
          if (output_riemanncubeAH_sdf||output_riemanncubeAH_ascii)
          {
            smooth_ah_r_(AH_riemanncube[c_AH],AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
          }
        }
      }
      else 
      {
        smooth_ah_r_b_(AH_theta0,AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
        smooth_ah_r_(AH_theta0,AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
        reg_ah_r_(AH_theta0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]); 
        smooth_ah_r_(AH_theta0,AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
        reg_ah_r_(AH_theta0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]); 

        if (output_moreAHquant_sdf||output_moreAHquant_ascii)
        {
          if (output_metricAH_cart_sdf||output_metricAH_cart_ascii)
          {
            smooth_ah_r_b_(AH_g0_xx[c_AH],AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
            smooth_ah_r_(AH_g0_xx[c_AH],AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
            reg_ah_r_(AH_g0_xx[c_AH],&AH_Nchi[c_AH],&AH_Nphi[c_AH]); 
            smooth_ah_r_(AH_g0_xx[c_AH],AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
            reg_ah_r_(AH_g0_xx[c_AH],&AH_Nchi[c_AH],&AH_Nphi[c_AH]);  
  
            smooth_ah_r_b_(AH_g0_xy[c_AH],AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
            smooth_ah_r_(AH_g0_xy[c_AH],AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
            reg_ah_r_(AH_g0_xy[c_AH],&AH_Nchi[c_AH],&AH_Nphi[c_AH]); 
            smooth_ah_r_(AH_g0_xy[c_AH],AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
            reg_ah_r_(AH_g0_xy[c_AH],&AH_Nchi[c_AH],&AH_Nphi[c_AH]);  
  
            smooth_ah_r_b_(AH_g0_xz[c_AH],AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
            smooth_ah_r_(AH_g0_xz[c_AH],AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
            reg_ah_r_(AH_g0_xz[c_AH],&AH_Nchi[c_AH],&AH_Nphi[c_AH]); 
            smooth_ah_r_(AH_g0_xz[c_AH],AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
            reg_ah_r_(AH_g0_xz[c_AH],&AH_Nchi[c_AH],&AH_Nphi[c_AH]);  
  
            smooth_ah_r_b_(AH_g0_yy[c_AH],AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
            smooth_ah_r_(AH_g0_yy[c_AH],AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
            reg_ah_r_(AH_g0_yy[c_AH],&AH_Nchi[c_AH],&AH_Nphi[c_AH]); 
            smooth_ah_r_(AH_g0_yy[c_AH],AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
            reg_ah_r_(AH_g0_yy[c_AH],&AH_Nchi[c_AH],&AH_Nphi[c_AH]);  
  
            smooth_ah_r_b_(AH_g0_yz[c_AH],AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
            smooth_ah_r_(AH_g0_yz[c_AH],AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
            reg_ah_r_(AH_g0_yz[c_AH],&AH_Nchi[c_AH],&AH_Nphi[c_AH]); 
            smooth_ah_r_(AH_g0_yz[c_AH],AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
            reg_ah_r_(AH_g0_yz[c_AH],&AH_Nchi[c_AH],&AH_Nphi[c_AH]);  
  
            smooth_ah_r_b_(AH_g0_zz[c_AH],AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
            smooth_ah_r_(AH_g0_zz[c_AH],AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
            reg_ah_r_(AH_g0_zz[c_AH],&AH_Nchi[c_AH],&AH_Nphi[c_AH]); 
            smooth_ah_r_(AH_g0_zz[c_AH],AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
            reg_ah_r_(AH_g0_zz[c_AH],&AH_Nchi[c_AH],&AH_Nphi[c_AH]);   
          }
          if (output_metricAH_sph_sdf||output_metricAH_sph_ascii)
          {
            smooth_ah_r_b_(AH_g0_chichi[c_AH],AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
            smooth_ah_r_(AH_g0_chichi[c_AH],AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
            reg_ah_r_(AH_g0_chichi[c_AH],&AH_Nchi[c_AH],&AH_Nphi[c_AH]); 
            smooth_ah_r_(AH_g0_chichi[c_AH],AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
            reg_ah_r_(AH_g0_chichi[c_AH],&AH_Nchi[c_AH],&AH_Nphi[c_AH]);  
  
            smooth_ah_r_b_(AH_g0_chiphi[c_AH],AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
            smooth_ah_r_(AH_g0_chiphi[c_AH],AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
            reg_ah_r_(AH_g0_chiphi[c_AH],&AH_Nchi[c_AH],&AH_Nphi[c_AH]); 
            smooth_ah_r_(AH_g0_chiphi[c_AH],AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
            reg_ah_r_(AH_g0_chiphi[c_AH],&AH_Nchi[c_AH],&AH_Nphi[c_AH]);  
  
            smooth_ah_r_b_(AH_g0_phiphi[c_AH],AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
            smooth_ah_r_(AH_g0_phiphi[c_AH],AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
            reg_ah_r_(AH_g0_phiphi[c_AH],&AH_Nchi[c_AH],&AH_Nphi[c_AH]); 
            smooth_ah_r_(AH_g0_phiphi[c_AH],AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
            reg_ah_r_(AH_g0_phiphi[c_AH],&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
          }
          if (output_kretschAH_sdf||output_kretschAH_ascii)
          {
            smooth_ah_r_b_(AH_kretsch[c_AH],AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
            smooth_ah_r_(AH_kretsch[c_AH],AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
            reg_ah_r_(AH_kretsch[c_AH],&AH_Nchi[c_AH],&AH_Nphi[c_AH]); 
            smooth_ah_r_(AH_kretsch[c_AH],AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
            reg_ah_r_(AH_kretsch[c_AH],&AH_Nchi[c_AH],&AH_Nphi[c_AH]);  
          }
          if (output_riemanncubeAH_sdf||output_riemanncubeAH_ascii)
          {
            smooth_ah_r_b_(AH_riemanncube[c_AH],AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
            smooth_ah_r_(AH_riemanncube[c_AH],AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
            reg_ah_r_(AH_riemanncube[c_AH],&AH_Nchi[c_AH],&AH_Nphi[c_AH]); 
            smooth_ah_r_(AH_riemanncube[c_AH],AH_w1[c_AH],&eps0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
            reg_ah_r_(AH_riemanncube[c_AH],&AH_Nchi[c_AH],&AH_Nphi[c_AH]);  
          }
        }

        eps0+=1;
      }
   }

   reg_ah_r_(AH_theta0,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
   if (output_moreAHquant_sdf||output_moreAHquant_ascii)
   {
      if (output_metricAH_cart_sdf||output_metricAH_cart_ascii)
      {
        reg_ah_r_(AH_g0_xx[c_AH],&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
        reg_ah_r_(AH_g0_xy[c_AH],&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
        reg_ah_r_(AH_g0_xz[c_AH],&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
        reg_ah_r_(AH_g0_yy[c_AH],&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
        reg_ah_r_(AH_g0_yz[c_AH],&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
        reg_ah_r_(AH_g0_zz[c_AH],&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
      }
      if (output_metricAH_sph_sdf||output_metricAH_sph_ascii)
      {
        reg_ah_r_(AH_g0_chichi[c_AH],&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
        reg_ah_r_(AH_g0_chiphi[c_AH],&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
        reg_ah_r_(AH_g0_phiphi[c_AH],&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
      }
      if (output_kretschAH_sdf||output_kretschAH_ascii)
      {
        reg_ah_r_(AH_kretsch[c_AH],&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
      }
      if (output_riemanncubeAH_sdf||output_riemanncubeAH_ascii)
      {
        reg_ah_r_(AH_riemanncube[c_AH],&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
      }
    }

   // compute resid as the root mean square of AH_theta values on the AH surface
   for (i=0, resid=0; i<np; i++) resid+=AH_theta0[i]*AH_theta0[i];
   resid=sqrt(resid/np);   

   if (!my_rank && dvtrace && num_trace<MAX_TRACE)
   {
      if (!((int)AH_ct[c_AH] % dvtrace))
      {
         num_trace++;
         AH_shape[0]=AH_Nchi[c_AH];
         AH_shape[1]=AH_Nphi[c_AH];
         AH_bbox[0]=0;
         AH_bbox[1]=M_PI;
         AH_bbox[2]=0;
         AH_bbox[3]=2*M_PI;
         rank=2; 
   
         size_name = snprintf(NULL, 0, "%sAH_%i_R_iter",AMRD_save_tag,c_AH+1) + 1;
         if ( (name = malloc(size_name) )==NULL) AMRD_stop("fill_theta_ahmetric...out of memory",""); 
         sprintf(name,"%sAH_%i_R_iter",AMRD_save_tag,c_AH+1);
         gft_out_bbox(name,AH_ct[c_AH],AH_shape,rank,AH_bbox,AH_R[c_AH]);
         free(name); name=NULL;

         size_name = snprintf(NULL, 0, "%sAH_%i_theta_iter",AMRD_save_tag,c_AH+1) + 1;
         if ( (name = malloc(size_name) )==NULL) AMRD_stop("fill_theta_ahmetric...out of memory",""); 
         sprintf(name,"%sAH_%i_theta_iter",AMRD_save_tag,c_AH+1);
         gft_out_bbox(name,AH_ct[c_AH],AH_shape,rank,AH_bbox,AH_theta0);
         free(name); name=NULL;


         //sprintf(name,"%sAH_%i_g0_xx_iter",AMRD_save_tag,c_AH+1);
         //gft_out_bbox(name,AH_ct[c_AH],AH_shape,rank,AH_bbox,AH_g0_xx[c_AH]);
         //sprintf(name,"%sAH_%i_g0_xy_iter",AMRD_save_tag,c_AH+1);
         //gft_out_bbox(name,AH_ct[c_AH],AH_shape,rank,AH_bbox,AH_g0_xy[c_AH]);
         //sprintf(name,"%sAH_%i_g0_xz_iter",AMRD_save_tag,c_AH+1);
         //gft_out_bbox(name,AH_ct[c_AH],AH_shape,rank,AH_bbox,AH_g0_xz[c_AH]);
         //sprintf(name,"%sAH_%i_g0_yy_iter",AMRD_save_tag,c_AH+1);
         //gft_out_bbox(name,AH_ct[c_AH],AH_shape,rank,AH_bbox,AH_g0_yy[c_AH]);
         //sprintf(name,"%sAH_%i_g0_yz_iter",AMRD_save_tag,c_AH+1);
         //gft_out_bbox(name,AH_ct[c_AH],AH_shape,rank,AH_bbox,AH_g0_yz[c_AH]);
         //sprintf(name,"%sAH_%i_g0_zz_iter",AMRD_save_tag,c_AH+1);
         //gft_out_bbox(name,AH_ct[c_AH],AH_shape,rank,AH_bbox,AH_g0_zz[c_AH]);
         //sprintf(name,"%sAH_%i_g0_chichi_iter",AMRD_save_tag,c_AH+1);
         //gft_out_bbox(name,AH_ct[c_AH],AH_shape,rank,AH_bbox,AH_g0_chichi[c_AH]);
         //sprintf(name,"%sAH_%i_g0_chiphi_iter",AMRD_save_tag,c_AH+1);
         //gft_out_bbox(name,AH_ct[c_AH],AH_shape,rank,AH_bbox,AH_g0_chiphi[c_AH]);
         //sprintf(name,"%sAH_%i_g0_phiphi_iter",AMRD_save_tag,c_AH+1);
         //gft_out_bbox(name,AH_ct[c_AH],AH_shape,rank,AH_bbox,AH_g0_phiphi[c_AH]);

      }

      AH_ct[c_AH]++;
   }

   return resid;
}



//-----------------------------------------------------------------------------
// driver routine
//-----------------------------------------------------------------------------
#define SMOOTH_R 0
#define SMOOTH_R_AFT 1
int find_apph(real *M, real *J, real *area, real *c_equat, real *c_polar, real *c_polar2, int use_R_ic, real *AH_min_resid, 
              real *ief_bh_r0,real *a_rot0, int *kerrads_background, int *AH_analytic_kerrads, real *ct)
{
   int iter,i,j,l,np,Lmax,Lmax_AH,is_ex;
   real resid,prev_resid,min_resid,min_R,c_R;
   int ltrace=2,numinc,skip_rest;
   real eps0,resid0,resid1,tol0,dx0[3],prev_AH_xc[2],eps1;
   static real min_tol,max_tol,first_c=1;
   int first=1;
   static int first_tol=1;

   *M=*J=0;
   *AH_min_resid=1e10;

   if (first_c) { for (l=0; l<MAX_BHS; l++) AH_ct[l]=0; first_c=0; }

   Lmax=PAMR_get_max_lev(PAMR_AMRH);
   if (Lmax>AH_Lmax[c_AH]) Lmax=AH_Lmax[c_AH];
   if (AH_Lmin[c_AH]>Lmax) return 0;

   np=AH_Nchi[c_AH]*AH_Nphi[c_AH];

   if (ltrace && my_rank==0)
   {
      printf("\n=========================================================================\n\n");
      printf("Searching for AH[%i] between L=%i and %i\n\n",c_AH+1,AH_Lmin[c_AH],Lmax);
   }

   eps0=eps1=AH_eps[c_AH];
   tol0=AH_tol[c_AH];
   if (eps1>1) eps1=eps1-1;

   // for found_AH:=use_R_ic=0, figure out an initial guess for AH_R
   if ((!use_R_ic)&&(!(*AH_analytic_kerrads)))
   {
      if (AH_rsteps[c_AH]==1) // for one radius iteration
      {
        for (i=0; i<np; i++) AH_R[c_AH][i]=AH_r0[c_AH];
      }
      else                                                                    // for several radii iterations
      {
         min_R=AH_r0[c_AH]; 
         min_resid=1e10;
         for (iter=1; iter<=(AH_rsteps[c_AH]+1); iter++)
         {
            c_R=AH_r0[c_AH]+(iter-1)*(AH_r1[c_AH]-AH_r0[c_AH])/(AH_rsteps[c_AH]-1);

            for (i=0; i<np; i++) AH_R[c_AH][i]=c_R;

            // compute initial theta values
            if (!fill_own(Lmax,ltrace,&first)) return 0;
            resid=fill_theta_ahmetric(AH_theta[c_AH],eps0,area,c_equat,c_polar,c_polar2,&is_ex,
                                      ief_bh_r0,a_rot0,kerrads_background,ct);

            if (is_ex)
            {
               printf("   is_ex on initial iteration\n");
               return 0;
            }
            if (resid<min_resid && !isnan(resid)) 
            {
               min_R=c_R;
               min_resid=resid;
            }
            if (ltrace>1 && my_rank==0)
               printf("   init iter %i: AH_r0=%lf ... |Theta|=%lf, Theta(mid)=%lf\n",iter,c_R,resid,AH_theta[c_AH][(np+1)/2]);
         }
         for (i=0; i<np; i++) AH_R[c_AH][i]=min_R;
      }
   }

   // save AH with minimum tol over all iterations in w2
   for (i=0; i<np; i++) AH_w2[c_AH][i]=AH_R[c_AH][i];
   for (i=0; i<2; i++) prev_AH_xc[i]=AH_xc[c_AH][i];

   int restarts=0;

   prev_resid=resid=tol0+1;
   iter=0;
   numinc=0;

   // (MAIN LOOP) while not yet AH_max_iter, or while residual is above tolerance
   while (iter < AH_max_iter[c_AH] && resid > tol0 && (numinc < AH_maxinc[c_AH] || AH_lambda[c_AH]>AH_lambda_min[c_AH]))   
   {
      iter++;
      skip_rest=0;

      if (numinc==AH_maxinc[c_AH])
      {
          AH_lambda[c_AH]*=(0.75);
          if (AH_lambda[c_AH]<AH_lambda_min[c_AH]) AH_lambda[c_AH]=AH_lambda_min[c_AH];
          if (my_rank==0 && ltrace) printf ("|Theta| increased for %i steps or out of bounds ... \n"
             "restoring AH shape to initial shape and decreasing lambda to %lf\n",numinc,AH_lambda[c_AH]);
          for (i=0; i<np; i++) AH_R[c_AH][i]=AH_w2[c_AH][i];
          for (i=0; i<2; i++) AH_xc[c_AH][i]=prev_AH_xc[i];
          numinc=0;
      }

      // do the adjustment within the iteration to keep changes small ... otherwise transformation is invalid
      adjust_ah_xc_(AH_R[c_AH],AH_xc[c_AH],&AH_Nchi[c_AH],&AH_Nphi[c_AH],&dx,&dy,&dz,&axisym);

      if (!fill_own(Lmax,ltrace,&first))
      {
         skip_rest=1;
         numinc=AH_maxinc[c_AH];
      }
      
      if (!skip_rest)
      {
      
         // compute theta values (and metric components at AH, OPTIONAL)
         resid=fill_theta_ahmetric(AH_theta[c_AH],eps0,area,c_equat,c_polar,c_polar2,&is_ex,
                                      ief_bh_r0,a_rot0,kerrads_background,ct);
      
//         if (my_rank==0)
//         {
//            printf("iter=%i:\n\n",iter);
//            for (i=0; i<AH_Nchi[c_AH]; i++)
//               for (j=0; j<AH_Nphi[c_AH]; j++) 
//                  printf("  i,j=%i,%i: theta=%lf, R=%lf\n",i,j,(AH_theta[c_AH])[i+j*AH_Nchi[c_AH]],(AH_R[c_AH])[i+j*AH_Nchi[c_AH]]);
//         }

         if (is_ex)
         {
            if (ltrace>1 && my_rank==0) printf("is_ex iter %i: AH_R(mid)[%i]=%lf ... |Theta|=%lf, Theta(mid)=%lf \n",iter,(np+1)/2,AH_R[c_AH][(np+1)/2],resid,AH_theta[c_AH][(np+1)/2]);
            if (ltrace>0 && my_rank==0) printf("AH too close to excision boundary (iter=%i,min_resid=%lf)\n",iter,*AH_min_resid);

            resid=tol0+10; goto fin;
         }

         // update AH_R[i], for each i point on the AH surface, by a fraction AH_lambda of AH_theta
         //(NOTE: AH_R decreases when theta positive, and AH_R increases when theta negative)
            for (i=0; i<np; i++)
            {
                  //printf("pre AH finder update:  i,j=%i,%i: AH_theta=%lf, R=%lf\n",i,j,(AH_theta[c_AH])[i+j*AH_Nchi[c_AH]],(AH_R[c_AH])[i+j*AH_Nchi[c_AH]]);
                  AH_R[c_AH][i]+=(-AH_lambda[c_AH]*AH_theta[c_AH][i]);
                  //printf("post AH finder update: AH_R[c_AH][i]=%lf\n",AH_R[c_AH][i+j*AH_Nchi[c_AH]]);
            }

         if (SMOOTH_R && eps1>0)
         {
            if (AH_xc[c_AH][1]==0) reg_ah_r_(AH_R[c_AH],&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
            smooth_ah_r_(AH_R[c_AH],AH_w1[c_AH],&eps1,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
         }

//LOOKING//
         /*if (iter>2 && resid>prev_resid) numinc++; else*/ numinc=0;

         prev_resid=resid;

         // print out AH_R(middle), |Theta|, Theta(middle), (xc,yc,zc)
         if (ltrace>1 && my_rank==0) printf("   iter %i: AH_R(mid)[%i]=%lf ... |Theta|=%lf, Theta(mid)=%lf, (xc,yc,zc)=(%lf,%lf,%lf) \n",iter,(np+1)/2,AH_R[c_AH][(np+1)/2],resid,AH_theta[c_AH][(np+1)/2],AH_xc[c_AH][0],AH_xc[c_AH][1],AH_xc[c_AH][2]);

         // checks for not-allowed negative AH_R 
         for (i=0; i<np; i++)
         {
            if (AH_R[c_AH][i]<0)
            {
              if (ltrace>0 && my_rank==0) printf("AH flowed to negative AH_R (iter=%i,min_resid=%lf)\n",iter,*AH_min_resid); resid=tol0+10; goto fin;
            }
         }

         if (resid < (*AH_min_resid) && !isnan(resid)) *AH_min_resid=resid;
      }
   }

   //compute Schwarzschild-AdS mass with given horizon area
   *M=(sqrt((*area)/4/M_PI))*(1+((*area)/4/M_PI)/AdS_L/AdS_L)/2;

   if (ltrace && my_rank==0)
   {
      if (resid>tol0) 
      {
         printf("\n ... failed to find an AH in %i iterations\n",iter);
         if (iter<AH_max_iter[c_AH]) printf("     (stopped because of increase in residual ... min_resid=%lf)\n",*AH_min_resid);
      }
      else 
      {
         printf("\n ... found an AH (to within %lf) in %i iterations ... \n",
                tol0,iter);
         printf("     horizon area approximation (calculating the area density and computing the integral over the apparent horizon, e.g., in Mathematica, is slightly more precise): %5.6lf,\n areal horizon radius (non-compactified): %5.6lf , horizon Schwarzschild-AdS mass with given area: %5.6lf \n",
              *area,sqrt((*area)/4/M_PI),*M);
         printf("     equat circum (x=0): %5.3lf,  polar circum 1 (y=0, i.e. phi=PI/2 && phi=3*PI/2): %5.3lf, and polar circum 2 (z=0, i.e. phi=0 && phi=PI): %5.3lf\n",
              *c_equat,*c_polar,*c_polar2);
      }

      printf("\n=========================================================================\n");
   }

   if (SMOOTH_R_AFT && eps1>0 && resid<=tol0)
   {
      if (AH_xc[c_AH][1]==0) reg_ah_r_(AH_R[c_AH],&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
      smooth_ah_r_(AH_R[c_AH],AH_w1[c_AH],&eps1,&AH_Nchi[c_AH],&AH_Nphi[c_AH]);
   }

   AH_ct[c_AH]+=100;

fin:

   if (resid>tol0)
   {
      for (i=0; i<np; i++) AH_R[c_AH][i]=AH_w2[c_AH][i];
      for (i=0; i<2; i++) AH_xc[c_AH][i]=prev_AH_xc[i];
      return 0; 
   }
   else return 1;
}

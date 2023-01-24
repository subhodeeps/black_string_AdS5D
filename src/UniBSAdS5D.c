//============================================================================
// Application interface functions for simulations of 
// perturbations of the uniform AdS-black string metric in 5 spacetime dimensions
// with a preserved SO(3) symmetry
//
// In coordinates (t,x=xi,y,omega_1,z), where
// x=xi in [0,1] is the radial coordinate along which we approach the AdS boundary.
// Each hypersurface of the uniform black string at fixed x is Schwarzschild-UniBSAdS5D. 
// (t,y,omega_1,z) are compactified Cartesian coordinates on each of these hypersurfaces.
// They are obtained by taking the uncompactified, horizon-penetrating 
// (e.g., Kerr-Schild or Gullstrand-Painleve-like) spherical coordinates (t,r,theta,phi), 
// using a compactified radial coordinate rho defined by r=2*rho/(1-rho^2),
// and then defining Cartesian coordinate y,omega_1,z as
// y=rho*cos(theta), omega_1= rho*sin(theta)*cos(phi), z=rho*sin(theta)*sin(phi)
//
// We only study perturbationswe that preserve an SO(3) symmetry, 
// i.e., that preserve the spherical symmetry of the xi-fixed slices.
// Thus, we can use the modified Cartoon method to reduce the simulation 
// on the hypersurface S at omega_1=z=0 and y in [0,1].
//=============================================================================

#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <pamr.h>
#include <amrd.h>
#include <math.h>
#include <m_util_r8.h>
#include <bbhutil.h>
#include <mpi.h>
#include "UniBSAdS5D.h"

const int app_dim=2; //It should match AMRD_dim!

int cp_version;

real *g_norms;

real *x,*y;
int shape[app_dim],ghost_width[2*app_dim],Nx,Ny,phys_bdy[2*app_dim],size,g_rank;
real base_bbox[2*app_dim],bbox[2*app_dim],dx,dy,dt,dx_Lc,dy_Lc;
int g_L;

//=============================================================================
// Background parameters
//=============================================================================

int unibs_background;
real AdS_L;
real ief_bh_r0;

//=============================================================================
// Excision parameters
//=============================================================================

int ex_reset_ybuf;
real ex_ybuf;
real ex_y,ex_yc;

//=============================================================================
// Enable/disable parameters
//=============================================================================

int background;

//=============================================================================
// Gauge parameters
//=============================================================================
int gauge_t;
real x1_t,x2_t,x3_t,x4_t,y3_t,y4_t,xi1_t,xi2_t,cbulk_t;
int gauge_i;
real x1_i,x2_i,x3_i,x4_i,y3_i,y4_i,xi1_i,xi2_i,cbulk_i;

//=============================================================================
// Pointers for grid functions
//=============================================================================

real *gb_tt,*gb_tt_n,*gb_tt_np1,*gb_tt_nm1; 
real *gb_tx,*gb_tx_n,*gb_tx_np1,*gb_tx_nm1; 
real *gb_ty,*gb_ty_n,*gb_ty_np1,*gb_ty_nm1;
real *gb_xx,*gb_xx_n,*gb_xx_np1,*gb_xx_nm1; 
real *gb_xy,*gb_xy_n,*gb_xy_np1,*gb_xy_nm1; 
real *gb_yy,*gb_yy_n,*gb_yy_np1,*gb_yy_nm1; 
real *gb_zz,*gb_zz_n,*gb_zz_np1,*gb_zz_nm1;

real *Hb_t,*Hb_t_n,*Hb_t_np1,*Hb_t_nm1;
real *Hb_x,*Hb_x_n,*Hb_x_np1,*Hb_x_nm1;
real *Hb_y,*Hb_y_n,*Hb_y_np1,*Hb_y_nm1;

real *gb_tt_t,*gb_tt_t_n;
real *gb_tx_t,*gb_tx_t_n;
real *gb_ty_t,*gb_ty_t_n;
real *gb_xx_t,*gb_xx_t_n;
real *gb_xy_t,*gb_xy_t_n;
real *gb_yy_t,*gb_yy_t_n;
real *gb_zz_t,*gb_zz_t_n;
real *Hb_t_t,*Hb_t_t_n;
real *Hb_x_t,*Hb_x_t_n;
real *Hb_y_t,*Hb_y_t_n;

real *mask,*mask_mg;
real *chr,*chr_mg;
real *gb_res;

real *efe_all_ires;
real *efe_tt_ires,*efe_tx_ires,*efe_ty_ires;
real *efe_xx_ires,*efe_xy_ires,*efe_yy_ires;
real *efe_zz_ires;

real *cl_res;

real *hb_t_res,*hb_i_res;
real *Hb_t_0,*Hb_x_0,*Hb_y_0;

//=============================================================================
// Grid function numbers: a grid function number (gfn) is an index,
// starting at 1, into an array of pointers to actual grid function data
//=============================================================================

int gb_tt_gfn,gb_tt_n_gfn,gb_tt_np1_gfn,gb_tt_nm1_gfn; 
int gb_tx_gfn,gb_tx_n_gfn,gb_tx_np1_gfn,gb_tx_nm1_gfn; 
int gb_ty_gfn,gb_ty_n_gfn,gb_ty_np1_gfn,gb_ty_nm1_gfn; 
int gb_xx_gfn,gb_xx_n_gfn,gb_xx_np1_gfn,gb_xx_nm1_gfn; 
int gb_xy_gfn,gb_xy_n_gfn,gb_xy_np1_gfn,gb_xy_nm1_gfn; 
int gb_yy_gfn,gb_yy_n_gfn,gb_yy_np1_gfn,gb_yy_nm1_gfn; 
int gb_zz_gfn,gb_zz_n_gfn,gb_zz_np1_gfn,gb_zz_nm1_gfn; 

int Hb_t_gfn,Hb_t_n_gfn,Hb_t_np1_gfn,Hb_t_nm1_gfn;
int Hb_x_gfn,Hb_x_n_gfn,Hb_x_np1_gfn,Hb_x_nm1_gfn;
int Hb_y_gfn,Hb_y_n_gfn,Hb_y_np1_gfn,Hb_y_nm1_gfn;

int gb_tt_t_gfn,gb_tt_t_n_gfn;
int gb_tx_t_gfn,gb_tx_t_n_gfn;
int gb_ty_t_gfn,gb_ty_t_n_gfn;
int gb_xx_t_gfn,gb_xx_t_n_gfn;
int gb_xy_t_gfn,gb_xy_t_n_gfn;
int gb_yy_t_gfn,gb_yy_t_n_gfn;
int gb_zz_t_gfn,gb_zz_t_n_gfn;
int Hb_t_t_gfn,Hb_t_t_n_gfn;
int Hb_x_t_gfn,Hb_x_t_n_gfn;
int Hb_y_t_gfn,Hb_y_t_n_gfn;

int mask_gfn,mask_mg_gfn;
int chr_gfn,chr_mg_gfn;

int gb_res_gfn;
int efe_all_ires_gfn;
int efe_tt_ires_gfn,efe_tx_ires_gfn,efe_ty_ires_gfn;
int efe_xx_ires_gfn,efe_xy_ires_gfn,efe_yy_ires_gfn;
int efe_zz_ires_gfn;

int cl_res_gfn;

int hb_t_res_gfn,hb_i_res_gfn;
int Hb_t_0_gfn,Hb_x_0_gfn,Hb_y_0_gfn;


//=============================================================================
// call after variables have been defined
//=============================================================================
void set_gfns(void)
{
	if ((gb_tt_gfn     = PAMR_get_gfn("gb_tt",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
   if ((gb_tt_nm1_gfn = PAMR_get_gfn("gb_tt",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
   if ((gb_tt_n_gfn   = PAMR_get_gfn("gb_tt",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
   if ((gb_tt_np1_gfn = PAMR_get_gfn("gb_tt",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

   if ((gb_tx_gfn     = PAMR_get_gfn("gb_tx",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
   if ((gb_tx_nm1_gfn = PAMR_get_gfn("gb_tx",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
   if ((gb_tx_n_gfn   = PAMR_get_gfn("gb_tx",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
   if ((gb_tx_np1_gfn = PAMR_get_gfn("gb_tx",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

   if ((gb_ty_gfn     = PAMR_get_gfn("gb_ty",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
   if ((gb_ty_nm1_gfn = PAMR_get_gfn("gb_ty",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
   if ((gb_ty_n_gfn   = PAMR_get_gfn("gb_ty",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
   if ((gb_ty_np1_gfn = PAMR_get_gfn("gb_ty",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

   if ((gb_xx_gfn     = PAMR_get_gfn("gb_xx",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
   if ((gb_xx_nm1_gfn = PAMR_get_gfn("gb_xx",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
   if ((gb_xx_n_gfn   = PAMR_get_gfn("gb_xx",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
   if ((gb_xx_np1_gfn = PAMR_get_gfn("gb_xx",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

   if ((gb_xy_gfn     = PAMR_get_gfn("gb_xy",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
   if ((gb_xy_nm1_gfn = PAMR_get_gfn("gb_xy",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
   if ((gb_xy_n_gfn   = PAMR_get_gfn("gb_xy",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
   if ((gb_xy_np1_gfn = PAMR_get_gfn("gb_xy",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

   if ((gb_yy_gfn     = PAMR_get_gfn("gb_yy",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
   if ((gb_yy_nm1_gfn = PAMR_get_gfn("gb_yy",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
   if ((gb_yy_n_gfn   = PAMR_get_gfn("gb_yy",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
   if ((gb_yy_np1_gfn = PAMR_get_gfn("gb_yy",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

   if ((gb_zz_gfn     = PAMR_get_gfn("gb_zz",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
   if ((gb_zz_nm1_gfn = PAMR_get_gfn("gb_zz",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
   if ((gb_zz_n_gfn   = PAMR_get_gfn("gb_zz",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
   if ((gb_zz_np1_gfn = PAMR_get_gfn("gb_zz",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

   if ((Hb_t_gfn      = PAMR_get_gfn("Hb_t",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
   if ((Hb_t_nm1_gfn  = PAMR_get_gfn("Hb_t",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
   if ((Hb_t_n_gfn    = PAMR_get_gfn("Hb_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
   if ((Hb_t_np1_gfn  = PAMR_get_gfn("Hb_t",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

   if ((Hb_x_gfn      = PAMR_get_gfn("Hb_x",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
   if ((Hb_x_nm1_gfn  = PAMR_get_gfn("Hb_x",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
   if ((Hb_x_n_gfn    = PAMR_get_gfn("Hb_x",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
   if ((Hb_x_np1_gfn  = PAMR_get_gfn("Hb_x",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

   if ((Hb_y_gfn      = PAMR_get_gfn("Hb_y",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
   if ((Hb_y_nm1_gfn  = PAMR_get_gfn("Hb_y",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
   if ((Hb_y_n_gfn    = PAMR_get_gfn("Hb_y",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
   if ((Hb_y_np1_gfn  = PAMR_get_gfn("Hb_y",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

   if ((gb_tt_t_gfn = PAMR_get_gfn("gb_tt_t",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
   if ((gb_tt_t_n_gfn = PAMR_get_gfn("gb_tt_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);

   if ((gb_tx_t_gfn = PAMR_get_gfn("gb_tx_t",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
   if ((gb_tx_t_n_gfn = PAMR_get_gfn("gb_tx_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);

   if ((gb_ty_t_gfn = PAMR_get_gfn("gb_ty_t",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
   if ((gb_ty_t_n_gfn = PAMR_get_gfn("gb_ty_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);

   if ((gb_xx_t_gfn = PAMR_get_gfn("gb_xx_t",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
   if ((gb_xx_t_n_gfn = PAMR_get_gfn("gb_xx_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);

   if ((gb_xy_t_gfn = PAMR_get_gfn("gb_xy_t",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
   if ((gb_xy_t_n_gfn = PAMR_get_gfn("gb_xy_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);

   if ((gb_yy_t_gfn = PAMR_get_gfn("gb_yy_t",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
   if ((gb_yy_t_n_gfn = PAMR_get_gfn("gb_yy_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);

   if ((gb_zz_t_gfn = PAMR_get_gfn("gb_zz_t",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
   if ((gb_zz_t_n_gfn = PAMR_get_gfn("gb_zz_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);

   if ((Hb_t_t_gfn  = PAMR_get_gfn("Hb_t_t",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
   if ((Hb_t_t_n_gfn  = PAMR_get_gfn("Hb_t_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);

   if ((Hb_x_t_gfn  = PAMR_get_gfn("Hb_x_t",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
   if ((Hb_x_t_n_gfn  = PAMR_get_gfn("Hb_x_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);

   if ((Hb_y_t_gfn  = PAMR_get_gfn("Hb_y_t",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
   if ((Hb_y_t_n_gfn  = PAMR_get_gfn("Hb_y_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);

   if ((mask_mg_gfn = PAMR_get_gfn("cmask",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
   if ((mask_gfn    = PAMR_get_gfn("cmask",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
   if ((chr_gfn     = PAMR_get_gfn("chr",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
   if ((chr_mg_gfn  = PAMR_get_gfn("chr",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);

   if ((gb_res_gfn    = PAMR_get_gfn("gb_res",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
   if ((efe_all_ires_gfn   = PAMR_get_gfn("efe_all_ires",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
   if ((efe_tt_ires_gfn    = PAMR_get_gfn("efe_tt_ires",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
   if ((efe_tx_ires_gfn    = PAMR_get_gfn("efe_tx_ires",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
   if ((efe_ty_ires_gfn    = PAMR_get_gfn("efe_ty_ires",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
   if ((efe_xx_ires_gfn    = PAMR_get_gfn("efe_xx_ires",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
   if ((efe_xy_ires_gfn    = PAMR_get_gfn("efe_xy_ires",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
   if ((efe_yy_ires_gfn    = PAMR_get_gfn("efe_yy_ires",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
   if ((efe_zz_ires_gfn    = PAMR_get_gfn("efe_zz_ires",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

   if ((cl_res_gfn   = PAMR_get_gfn("cl_res",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

   if ((hb_t_res_gfn  = PAMR_get_gfn("hb_t_res",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
   if ((hb_i_res_gfn  = PAMR_get_gfn("hb_i_res",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

   if ((Hb_t_0_gfn  = PAMR_get_gfn("Hb_t_0",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
   if ((Hb_x_0_gfn  = PAMR_get_gfn("Hb_x_0",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
   if ((Hb_y_0_gfn  = PAMR_get_gfn("Hb_y_0",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

   g_norms=AMRD_get_global_norms();
}


//=============================================================================
// call with valid iter to set up globals:
//=============================================================================
void ldptr_bbox(void)
{
	real dx0[2];
  	static int first=1;	

	if (first) 
  	{
      first=0; 
      set_gfns();
      PAMR_get_global_bbox(base_bbox);
      if (PAMR_get_max_lev(PAMR_AMRH)>1) PAMR_get_dxdt(2,dx0,&dt); else PAMR_get_dxdt(1,dx0,&dt);
      dx_Lc=dx0[0];
      dy_Lc=dx0[1];
  	}	
  	PAMR_get_g_rank(&g_rank);
  	PAMR_get_g_shape(shape);
  	PAMR_get_g_bbox(bbox);
  	PAMR_get_g_ghost_width(ghost_width);
  	PAMR_get_g_level(&g_L);
  	PAMR_get_dxdt(g_L,dx0,&dt);
  	dx=dx0[0];
  	dy=dx0[1];
  	if ((bbox[0]-base_bbox[0])<dx/2) phys_bdy[0]=1; else phys_bdy[0]=0;
  	if ((base_bbox[1]-bbox[1])<dx/2) phys_bdy[1]=1; else phys_bdy[1]=0;
  	if ((bbox[2]-base_bbox[2])<dy/2) phys_bdy[2]=1; else phys_bdy[2]=0;
  	if ((base_bbox[3]-bbox[3])<dy/2) phys_bdy[3]=1; else phys_bdy[3]=0;
  	Nx=shape[0];
  	Ny=shape[1];
  	size=Nx*Ny;
}

void ldptr(void)
{

	real *x0[2],*gfs[PAMR_MAX_GFNS];

   ldptr_bbox();
   PAMR_get_g_x(x0);   
   x=x0[0];
   y=x0[1];  
   PAMR_get_g_gfs(gfs);   
   gb_tt     = gfs[gb_tt_gfn-1];
   gb_tt_n   = gfs[gb_tt_n_gfn-1];
   gb_tt_np1 = gfs[gb_tt_np1_gfn-1];
   gb_tt_nm1 = gfs[gb_tt_nm1_gfn-1];   
   gb_tx     = gfs[gb_tx_gfn-1];
   gb_tx_n   = gfs[gb_tx_n_gfn-1];
   gb_tx_np1 = gfs[gb_tx_np1_gfn-1];
   gb_tx_nm1 = gfs[gb_tx_nm1_gfn-1];   
   gb_ty     = gfs[gb_ty_gfn-1];
   gb_ty_n   = gfs[gb_ty_n_gfn-1];
   gb_ty_np1 = gfs[gb_ty_np1_gfn-1];
   gb_ty_nm1 = gfs[gb_ty_nm1_gfn-1];
   gb_xx     = gfs[gb_xx_gfn-1];
   gb_xx_n   = gfs[gb_xx_n_gfn-1];
   gb_xx_np1 = gfs[gb_xx_np1_gfn-1];
   gb_xx_nm1 = gfs[gb_xx_nm1_gfn-1];   
   gb_xy     = gfs[gb_xy_gfn-1];
   gb_xy_n   = gfs[gb_xy_n_gfn-1];
   gb_xy_np1 = gfs[gb_xy_np1_gfn-1];
   gb_xy_nm1 = gfs[gb_xy_nm1_gfn-1];
   gb_yy     = gfs[gb_yy_gfn-1];
   gb_yy_n   = gfs[gb_yy_n_gfn-1];
   gb_yy_np1 = gfs[gb_yy_np1_gfn-1];
   gb_yy_nm1 = gfs[gb_yy_nm1_gfn-1];
   gb_zz     = gfs[gb_zz_gfn-1];
   gb_zz_n   = gfs[gb_zz_n_gfn-1];
   gb_zz_np1 = gfs[gb_zz_np1_gfn-1];
   gb_zz_nm1 = gfs[gb_zz_nm1_gfn-1];
   Hb_t      = gfs[Hb_t_gfn-1];
   Hb_t_n    = gfs[Hb_t_n_gfn-1];
   Hb_t_nm1  = gfs[Hb_t_nm1_gfn-1];
   Hb_t_np1  = gfs[Hb_t_np1_gfn-1];    
   Hb_x      = gfs[Hb_x_gfn-1];
   Hb_x_n    = gfs[Hb_x_n_gfn-1];
   Hb_x_nm1  = gfs[Hb_x_nm1_gfn-1];
   Hb_x_np1  = gfs[Hb_x_np1_gfn-1];    
   Hb_y      = gfs[Hb_y_gfn-1];
   Hb_y_n    = gfs[Hb_y_n_gfn-1];
   Hb_y_nm1  = gfs[Hb_y_nm1_gfn-1];
   Hb_y_np1  = gfs[Hb_y_np1_gfn-1];
   gb_tt_t = gfs[gb_tt_t_gfn-1];
   gb_tt_t_n = gfs[gb_tt_t_n_gfn-1];   
   gb_tx_t = gfs[gb_tx_t_gfn-1];
   gb_tx_t_n = gfs[gb_tx_t_n_gfn-1];   
   gb_ty_t = gfs[gb_ty_t_gfn-1];
   gb_ty_t_n = gfs[gb_ty_t_n_gfn-1];
   gb_xx_t = gfs[gb_xx_t_gfn-1];
   gb_xx_t_n = gfs[gb_xx_t_n_gfn-1];   
   gb_xy_t = gfs[gb_xy_t_gfn-1];
   gb_xy_t_n = gfs[gb_xy_t_n_gfn-1];
   gb_yy_t = gfs[gb_yy_t_gfn-1];
   gb_yy_t_n = gfs[gb_yy_t_n_gfn-1];
   gb_zz_t = gfs[gb_zz_t_gfn-1];
   gb_zz_t_n = gfs[gb_zz_t_n_gfn-1];   
   Hb_t_t  = gfs[Hb_t_t_gfn-1];
   Hb_t_t_n  = gfs[Hb_t_t_n_gfn-1];
   Hb_x_t  = gfs[Hb_x_t_gfn-1];
   Hb_x_t_n  = gfs[Hb_x_t_n_gfn-1];
   Hb_y_t  = gfs[Hb_y_t_gfn-1];
   Hb_y_t_n  = gfs[Hb_y_t_n_gfn-1];
   gb_res    = gfs[gb_res_gfn-1];
   efe_all_ires  = gfs[efe_all_ires_gfn-1];
   efe_tt_ires  = gfs[efe_tt_ires_gfn-1];
   efe_tx_ires  = gfs[efe_tx_ires_gfn-1];
   efe_ty_ires  = gfs[efe_ty_ires_gfn-1];
   efe_xx_ires  = gfs[efe_xx_ires_gfn-1];
   efe_xy_ires  = gfs[efe_xy_ires_gfn-1];
   efe_yy_ires  = gfs[efe_yy_ires_gfn-1];
   efe_zz_ires  = gfs[efe_zz_ires_gfn-1];
   cl_res   = gfs[cl_res_gfn-1];
   hb_t_res  = gfs[hb_t_res_gfn-1];
   hb_i_res  = gfs[hb_i_res_gfn-1];
   Hb_t_0  = gfs[Hb_t_0_gfn-1];
   Hb_x_0  = gfs[Hb_x_0_gfn-1];
   Hb_y_0  = gfs[Hb_y_0_gfn-1]; 
}

//=============================================================================
// PAMR_get_dxdt() only works with AMR hierarchy levels ... here we use
// lambda for dt, but this only works if rhosp=rhotm
//=============================================================================
void ldptr_mg(void)
{
	real lambda;

	ldptr();

	dx=x[1]-x[0]; dy=y[1]-y[0];
	PAMR_get_lambda(&lambda);
	dt=lambda*dx;
}

//=============================================================================
// utility routines
//=============================================================================
real norm_l2(real *f, real *cmask, real *chr)
{
   int i;
   real norm=0;
   int sum=0;

   for (i=0; i<Nx*Ny; i++) 
      if (cmask[i]==AMRD_CMASK_ON && (chr[i]!=AMRD_ex)) { sum++; norm+=f[i]*f[i]; }

   if (!sum) sum=1;
   return (sqrt(norm/sum));
}

void print_background_info(real AdS_L, real ief_bh_r0)
{
   real rh,M0,rhoh,yh;

   rh=-pow(AdS_L,2)
            /(pow(3,(1.0/3.0)))
            /(pow((9*pow(AdS_L,2)*(ief_bh_r0/2)+sqrt(3.0)*sqrt(pow(AdS_L,6)+27*pow(AdS_L,4)*pow((ief_bh_r0/2),2))),(1.0/3.0)))
            +(pow((9*pow(AdS_L,2)*(ief_bh_r0/2)+sqrt(3.0)*sqrt(pow(AdS_L,6)+27*pow(AdS_L,4)*pow((ief_bh_r0/2),2))),(1.0/3.0)))
            /(pow(3,(2.0/3.0)));
   M0=ief_bh_r0/2;
   rhoh=(-1 + sqrt(1 + pow(rh,2)))/rh;
   yh=rhoh; //y coordinate of the horizon; recall that we set the other Cartesian coordinates on Schwarzschild spacetime to 0 in this code
   printf("\n==============Information about background spacetime===================\n"
          "\nUniform black string AdS5D initial data with Schwarzschild-AdS4D slices\n"
          "Schwarzschild-AdS4D parameters:\n"
          "r0/L=%lf\n"
          "Mass parameter M0 = r0/2 = rh*(1+rh^2/L^2)/2 = %lf\n"
          "Horizon radius in uncompactified Schwarzschild coordinates: rh/L=%lf\n"
          "Horizon radius in compactified Schwarzschild coordinates: rhoh=yh=%lf )\n\n" 
          ,ief_bh_r0/AdS_L,M0,rh/AdS_L,yh);
}

real set_excision_radius(real AdS_L, real ief_bh_r0, real ex_ybuf)
{
   real rh,M0,rhoh,yh,ex_y;

   rh=-pow(AdS_L,2)
            /(pow(3,(1.0/3.0)))
            /(pow((9*pow(AdS_L,2)*(ief_bh_r0/2)+sqrt(3.0)*sqrt(pow(AdS_L,6)+27*pow(AdS_L,4)*pow((ief_bh_r0/2),2))),(1.0/3.0)))
            +(pow((9*pow(AdS_L,2)*(ief_bh_r0/2)+sqrt(3.0)*sqrt(pow(AdS_L,6)+27*pow(AdS_L,4)*pow((ief_bh_r0/2),2))),(1.0/3.0)))
            /(pow(3,(2.0/3.0)));
   M0=ief_bh_r0/2;
   rhoh=(-1 + sqrt(1 + pow(rh,2)))/rh;
   yh=rhoh; //y coordinate of the horizon; recall that we set the other Cartesian coordinates on Schwarzschild spacetime to 0 in this code
   ex_y=yh*(1-ex_ybuf);

   return ex_y;
}

void print_excision_info(real ex_yc, real ex_ybuf, real ex_y)
{
   printf("\n==============Information about excision===================\n"
          "Cylindrical excision along x inside the black string horizon\n"
          "Excision center is at yc=%lf\n"
          "Excision buffer (i.e. size of the evolved region within the horizon) is ex_ybuf=%lf\n\n"
          ,ex_yc,ex_ybuf);
   printf("The apparent horizon finder is not implemented yet, so we use fixed excision inside the background black string horizon\n"
          "The excision surface is at the y-coordinate value ex_y=%lf\n\n",ex_y);
}

// to be callable from fortran
void check_nan_(real *x, int *is_nan)
{
    if (isnan(*x)) *is_nan=1; else *is_nan=0;
}

//=============================================================================
// Routines required by amrd:
//=============================================================================

//=============================================================================
// Returns 0 to use default mechanism, or is expected to calculate
// the correct initial hierarchy and return 1:
//=============================================================================
int UniBSAdS5D_id(void)
{
   return 0;
}

//=============================================================================
// Sets custom parameters, variables, etc. Split up into two segments,
// one called before the pamr context is initialized and standard
// parameters are read, and the other afterwards
//=============================================================================
void UniBSAdS5D_var_pre_init(char *pfile)
{

	AMRD_echo_params=1; AMRD_int_param(pfile,"echo_params",&AMRD_echo_params,1);
   cp_version=UNIBSADS5D_CP_VERSION; AMRD_int_param(pfile,"cp_version",&cp_version,1);

	return;
}

void UniBSAdS5D_var_post_init(char *pfile)
{

	if (my_rank==0)
   {
      printf("===================================================================\n");
      printf("Reading UniBSAdS5D parameters:\n\n");
   }

//=============================================================================
// Read background parameters
//=============================================================================

   unibs_background=1; AMRD_int_param(pfile,"unibs_background",&unibs_background,1);
   if (unibs_background==2) AMRD_stop("Schwarzschild slices in Gullstrand-Painleve coords not yet implemented\n","");
   AdS_L=1.0; AMRD_real_param(pfile,"AdS_L",&AdS_L,1);
   ief_bh_r0=0.0; AMRD_real_param(pfile,"ief_bh_r0",&ief_bh_r0,1);

//=============================================================================
// Read excision parameters
//=============================================================================

   ex_reset_ybuf=0; AMRD_int_param(pfile,"ex_reset_ybuf",&ex_reset_ybuf,1);
   if (!AMRD_cp_restart || ex_reset_ybuf)
   {
   	AMRD_real_param(pfile,"ex_ybuf",&ex_ybuf,1);
   	if (ex_ybuf<0 || ex_ybuf>1 ) printf("WARNING ... ex_ybuf=%lf is outside of standard range\n",ex_ybuf);
   }

   // initialise y-coordinate of the centre of excised region
   // and coordinate-distance of excision boundary from centre of excised region
   if (!AMRD_cp_restart)
   {
      ex_yc=0;
      ex_y=0;
   }

//=============================================================================
// Read enable/disable parameters
//=============================================================================

   background=0; AMRD_int_param(pfile,"background",&background,1);

//=============================================================================
// Read gauge parameters
//=============================================================================
   gauge_t=0; AMRD_int_param(pfile,"gauge_t",&gauge_t,1);
   x1_t=0; AMRD_real_param(pfile,"x1_t",&x1_t,1);
   x2_t=0; AMRD_real_param(pfile,"x2_t",&x2_t,1);
   x3_t=0; AMRD_real_param(pfile,"x3_t",&x3_t,1);
   x4_t=0; AMRD_real_param(pfile,"x4_t",&x4_t,1);
   y3_t=0; AMRD_real_param(pfile,"y3_t",&x3_t,1);
   y4_t=0; AMRD_real_param(pfile,"y4_t",&x4_t,1);
   xi1_t=0; AMRD_real_param(pfile,"xi1_t",&xi1_t,1);
   xi2_t=0; AMRD_real_param(pfile,"xi2_t",&xi2_t,1);
   cbulk_t=0; AMRD_real_param(pfile,"cbulk_t",&cbulk_t,1);
   gauge_i=0; AMRD_int_param(pfile,"gauge_i",&gauge_t,1);
   x1_i=0; AMRD_real_param(pfile,"x1_i",&x1_i,1);
   x2_i=0; AMRD_real_param(pfile,"x2_i",&x2_i,1);
   x3_i=0; AMRD_real_param(pfile,"x3_i",&x3_i,1);
   x4_i=0; AMRD_real_param(pfile,"x4_i",&x4_i,1);
   y3_i=0; AMRD_real_param(pfile,"y3_i",&x3_i,1);
   y4_i=0; AMRD_real_param(pfile,"y4_i",&x4_i,1);
   xi1_i=0; AMRD_real_param(pfile,"xi1_i",&xi1_i,1);
   xi2_i=0; AMRD_real_param(pfile,"xi2_i",&xi2_i,1);
   cbulk_i=0; AMRD_real_param(pfile,"cbulk_i",&cbulk_i,1);

//=============================================================================
// Print info
//=============================================================================
   if (my_rank==0) print_background_info(AdS_L,ief_bh_r0);

//=============================================================================
// Apply excision
//=============================================================================
   if (AMRD_do_ex==0) AMRD_stop("require excision to be on","");
   ex_y=set_excision_radius(AdS_L,ief_bh_r0,ex_ybuf);
   if (my_rank==0) print_excision_info(ex_yc,ex_ybuf,ex_y);
   PAMR_excision_on("chr",&UniBSAdS5D_fill_ex_mask,AMRD_ex,1);

   if (my_rank==0) printf("===================================================================\n");
   return;

}

//=============================================================================
// Sets all variables to their 'zero' values:
//=============================================================================
void UniBSAdS5D_AMRH_var_clear(void)
{

   ldptr();
   init_gbhb_(gb_tt_n,gb_tx_n,gb_ty_n,
              gb_xx_n,gb_xy_n,gb_yy_n,
              gb_zz_n,
              Hb_t_n,Hb_x_n,Hb_y_n,
              &AdS_L,x,y,chr,&AMRD_ex,&Nx,&Ny);

	return;
}

//=============================================================================
// Initial data for free fields: (at tn=2) ... following vars also used in 
// t0_cnst_data
//=============================================================================
void UniBSAdS5D_free_data(void)
{
   //specify initial values of gb's. 
   //These are "free" since we are not solving the constraints of GR for now.
   return;
}

//=============================================================================
// Initialize any "elliptic_vars_t0" post construction of MGH, but before
// the start of vcycling.
//=============================================================================
void UniBSAdS5D_elliptic_vars_t0_init(void)
{

}

//=============================================================================
// Initial constraint data --- called after each MG iteration.
//=============================================================================
void UniBSAdS5D_t0_cnst_data(void)
{
	return;
}

//=============================================================================
// Calculations prior to saving info to disk.
//
// NOTE: at this point, the time sequence is: n,nm1,np1
//=============================================================================

void UniBSAdS5D_pre_io_calc(void)
{
	return;
}

//=============================================================================
// Returns some norm of the residual for the evolution variables ... just
// use the value from the most recent evolution step
//=============================================================================
#define LIN_ZERO_BND 1
int lin_zero_bnd_all=1;
real UniBSAdS5D_evo_residual(void)
{
   real l2norm=0,l2norm_gb,l2norm_hb_t,l2norm_hb_i;
   int is_nan;

   ldptr();    
   if (LIN_ZERO_BND) 
   {
      lin_zero_bnd_res_(gb_res,phys_bdy,&lin_zero_bnd_all,&Nx,&Ny,&app_dim);
      lin_zero_bnd_res_(hb_t_res,phys_bdy,&lin_zero_bnd_all,&Nx,&Ny,&app_dim);
      lin_zero_bnd_res_(hb_i_res,phys_bdy,&lin_zero_bnd_all,&Nx,&Ny,&app_dim);
   }
   l2norm_gb=norm_l2(gb_res,mask,chr);
   l2norm_hb_t=norm_l2(hb_t_res,mask,chr);
   l2norm_hb_i=norm_l2(hb_i_res,mask,chr);
   if (!background) l2norm+=(l2norm_gb+l2norm_hb_t+l2norm_hb_i);
   check_nan_(&l2norm,&is_nan);    
   if (is_nan)
   {
      printf("\nl2norm_gb=%lf, l2norm_hb_t=%lf, l2norm_hb_i=%lf\n",
                l2norm_gb,l2norm_hb_t,l2norm_hb_i);
      printf("[Nx,Ny]=[%i,%i],L=%i\n",Nx,Ny,g_L);
      AMRD_stop("l2norm is nan ... stopping","");
      l2norm=0;
   }

   return l2norm;
}

//=============================================================================
// Performs 1 iteration of the evolution equations 
//
// NOTE: at this point, the time sequence is: np1,n,nm1
//=============================================================================
void UniBSAdS5D_evolve(int iter)
{
   int ltrace=0;
   real ct;

   ldptr();    
   ct=PAMR_get_time(g_L);  
   if (ltrace) printf("UniBSAdS5D_evolve: iter=%i , time=%lf, lev=%i, rank=%i\n",iter,ct,g_L,my_rank); 

   if (!background)
   {
      hb_t_evo_(hb_t_res,
                gb_tt_np1,gb_tt_n,gb_tt_nm1,
                gb_tx_np1,gb_tx_n,gb_tx_nm1,
                gb_ty_np1,gb_ty_n,gb_ty_nm1,
                gb_xx_np1,gb_xx_n,gb_xx_nm1,
                gb_xy_np1,gb_xy_n,gb_xy_nm1,
                gb_yy_np1,gb_yy_n,gb_yy_nm1,
                gb_zz_np1,gb_zz_n,gb_zz_nm1,
                Hb_t_np1,Hb_t_n,Hb_t_nm1,
                Hb_x_np1,Hb_x_n,Hb_x_nm1,
                Hb_y_np1,Hb_y_n,Hb_y_nm1,
                &AdS_L,x,y,&dt,chr,&AMRD_ex,
                phys_bdy,ghost_width,&Nx,&Ny,&app_dim,
                Hb_t_0,Hb_x_0,Hb_y_0,
                &gauge_t,&ct,
                &x1_t,&x2_t,&x3_t,&x4_t,&y3_t,&y4_t,&xi1_t,&xi2_t,
                &cbulk_t,&unibs_background);  
      hb_i_evo_(hb_i_res,
                gb_tt_np1,gb_tt_n,gb_tt_nm1,
                gb_tx_np1,gb_tx_n,gb_tx_nm1,
                gb_ty_np1,gb_ty_n,gb_ty_nm1,
                gb_xx_np1,gb_xx_n,gb_xx_nm1,
                gb_xy_np1,gb_xy_n,gb_xy_nm1,
                gb_yy_np1,gb_yy_n,gb_yy_nm1,
                gb_zz_np1,gb_zz_n,gb_zz_nm1,
                Hb_t_np1,Hb_t_n,Hb_t_nm1,
                Hb_x_np1,Hb_x_n,Hb_x_nm1,
                Hb_y_np1,Hb_y_n,Hb_y_nm1,
                &AdS_L,x,y,&dt,chr,&AMRD_ex,
                phys_bdy,ghost_width,&Nx,&Ny,&app_dim,
                Hb_t_0,Hb_x_0,Hb_y_0,
                &gauge_i,&ct,
                &x1_i,&x2_i,&x3_i,&x4_i,&y3_i,&y4_i,&xi1_i,&xi2_i,
                &cbulk_i,&unibs_background);
//      g_evo_opt_(gb_res,phi1_res,cl_res,
//                 gb_tt_np1,gb_tt_n,gb_tt_nm1,
//                 gb_tx_np1,gb_tx_n,gb_tx_nm1,
//                 gb_ty_np1,gb_ty_n,gb_ty_nm1,
//                 gb_xx_np1,gb_xx_n,gb_xx_nm1,
//                 gb_xy_np1,gb_xy_n,gb_xy_nm1,     
//                 gb_yy_np1,gb_yy_n,gb_yy_nm1,
//                 gb_zz_np1,gb_zz_n,gb_zz_nm1,
//                 Hb_t_np1,Hb_t_n,Hb_t_nm1,
//                 Hb_x_np1,Hb_x_n,Hb_x_nm1,
//                 Hb_y_np1,Hb_y_n,Hb_y_nm1,
//                 &AdS_L,x,y,&dt,chr,&AMRD_ex,
//                 phys_bdy,ghost_width,&Nx,&Ny,&app_dim,
//                 &background,&kappa_cd,&rho_cd,
//                 &interptype,&i_shift,&regtype,
//                 &diss_kmax,tfunction,
//                 &ief_bh_r0,&a_rot0,&unibs_background);
   }

	return;
}

//=============================================================================
// sets excision mask (NO ITERATOR, SO DON'T LOAD POINTERS!!!)
// when using a radial coordinate for the code,
// the only excised region is that inside the horizon
//=============================================================================
void UniBSAdS5D_fill_ex_mask(real *mask, int dim, int *shape, real *bbox, real excised)
{
	real x,y,dx,dy,yp;
	int i,j,ind;

	dx=(bbox[1]-bbox[0])/(shape[0]-1);
   dy=(bbox[3]-bbox[2])/(shape[1]-1);

   for (i=0; i<shape[0]; i++)
   {
      x=bbox[0]+i*dx;
      for (j=0; j<shape[1]; j++)
      {
         y=bbox[2]+j*dy;
         ind=i+shape[0]*j;  
         mask[ind]=excised-1;

         yp=y-ex_yc;
         if (yp<ex_y)
         {
            mask[ind]=excised;
	      }
      }
   }
}

//=============================================================================
//=============================================================================
void UniBSAdS5D_fill_bh_bboxes(real *bbox, int *num, int max_num)
{   
    *num=1;
}

//=============================================================================
// The following routine searches for AH's, manages excision, 
// and if t==0, determines past time level information using evolve
//
// NOTE: at this point, the time sequence is: n,nm1,np1 (unless t=0)
//=============================================================================
void UniBSAdS5D_pre_tstep(int L)
{
	return;
}

//=============================================================================
// The following routine prints diagnostic quantities
//
// NOTE: at this point, the time sequence is: n,nm1,np1 (unless t=0)
//=============================================================================
void UniBSAdS5D_post_tstep(int L)
{
	return;
}

#define ACTION_RELAX 1
#define ACTION_LOP 3
#define ACTION_RESIDUAL 2
//=============================================================================
// Returns some norm of the residual for the MG variables, *AND* 
// stores the point-wise residual in "f_res" for each MG variable "f" (for
// computing new RHS's)
//=============================================================================
real UniBSAdS5D_MG_residual(void)
{
   real norm=0;

   return norm;
}


//=============================================================================
// Performs 1 relaxation sweep of the MG equations, and returns an estimate
// of the norm of the residual.
//=============================================================================
real UniBSAdS5D_MG_relax(void)
{
   real norm=0;

   return norm;
}

//=============================================================================
// Computes the differential operator L acting on the current grid,
// in a region specified by the grid function "mask". Stores the result
// in "f_lop" for each MG variable "f"
//=============================================================================
void UniBSAdS5D_L_op(void)
{ 
    return;
}

//=============================================================================
//=============================================================================
void UniBSAdS5D_scale_tre(void)
{
}

//=============================================================================
// post-regrid initialization of constant functions
//=============================================================================
void UniBSAdS5D_post_regrid(void)
{
	return;
}

//=============================================================================
//check-pointing
//=============================================================================
#define CP_DATA_SIZE 50000
void UniBSAdS5D_copy_block(char *p, char **q, int n, int dir, int *tot_size)
{
   char *p0;
   int n0; 
   if (n==0) return;   
   if ((*tot_size+n) > (CP_DATA_SIZE))
      AMRD_stop("UniBSAdS5D_copy_block: error ... CP_DATA_SIZE too small\n","");
   *tot_size+=n;   
   n0=n;
   p0=p;   
   if (dir==AMRD_CP_SAVE) while(n0--) *(*q)++=*p0++;
   else while(n0--) *p0++=*(*q)++; 
   return;
}

void UniBSAdS5D_cp(int dir, char *data)
{
   int size=0;    
   if (dir==AMRD_CP_SAVE)
   {
      cp_version=UNIBSADS5D_CP_VERSION;
      UniBSAdS5D_copy_block((char *)&cp_version,&data,sizeof(int),dir,&size);
   }
}

//=============================================================================
int main(int argc, char **argv)
{
   amrd_set_app_user_cp_hook(UniBSAdS5D_cp,CP_DATA_SIZE);
   amrd_set_app_pre_tstep_hook(UniBSAdS5D_pre_tstep);
   amrd_set_elliptic_vars_t0_init(UniBSAdS5D_elliptic_vars_t0_init);
   amrd(argc,argv,&UniBSAdS5D_id,&UniBSAdS5D_var_pre_init,
      &UniBSAdS5D_var_post_init, &UniBSAdS5D_AMRH_var_clear,
      &UniBSAdS5D_free_data, &UniBSAdS5D_t0_cnst_data,
      &UniBSAdS5D_evo_residual, &UniBSAdS5D_MG_residual,
      &UniBSAdS5D_evolve, &UniBSAdS5D_MG_relax, &UniBSAdS5D_L_op, 
      &UniBSAdS5D_pre_io_calc, &UniBSAdS5D_scale_tre, 
      &UniBSAdS5D_post_regrid, &UniBSAdS5D_post_tstep,
      &UniBSAdS5D_fill_ex_mask, &UniBSAdS5D_fill_bh_bboxes);
}











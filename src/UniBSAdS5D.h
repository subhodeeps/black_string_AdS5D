#ifndef _UNIBSADS5D_H
#define _UNIBSADS5D_H

/*============================================================================= */
/* global variables and prototypes for UniBSAdS5D                               */
/* "extern" is short for external. It means that the declared variable          */
/* is defined (i.e., memory is allocated)                                       */
/* outside of the current file (in UniBSAdS5D.c in this code).                  */
/* It is necessary to declare these variables here if we want to use them in    */
/* files other than the one where they are declared.                            */
/*============================================================================= */

#define UNIBSADS5D_CP_VERSION 1

//extern real AdS_L;
//
//extern real *x,*y;
//extern int Nx,Ny;
//extern real ex;

//extern int cp_version;
//
//extern int gb_tt_gfn,gb_tt_n_gfn,gb_tt_np1_gfn,gb_tt_nm1_gfn; 
//extern int gb_tx_gfn,gb_tx_n_gfn,gb_tx_np1_gfn,gb_tx_nm1_gfn; 
//extern int gb_ty_gfn,gb_ty_n_gfn,gb_ty_np1_gfn,gb_ty_nm1_gfn; 
//extern int gb_xx_gfn,gb_xx_n_gfn,gb_xx_np1_gfn,gb_xx_nm1_gfn; 
//extern int gb_xy_gfn,gb_xy_n_gfn,gb_xy_np1_gfn,gb_xy_nm1_gfn; 
//extern int gb_yy_gfn,gb_yy_n_gfn,gb_yy_np1_gfn,gb_yy_nm1_gfn; 
//extern int gb_zz_gfn,gb_zz_n_gfn,gb_zz_np1_gfn,gb_zz_nm1_gfn; 
//
//extern int Hb_t_gfn,Hb_t_n_gfn,Hb_t_np1_gfn,Hb_t_nm1_gfn;
//extern int Hb_x_gfn,Hb_x_n_gfn,Hb_x_np1_gfn,Hb_x_nm1_gfn;
//extern int Hb_y_gfn,Hb_y_n_gfn,Hb_y_np1_gfn,Hb_y_nm1_gfn;
//
//extern real *gb_tt_t,*gb_tt_t_n;
//extern real *gb_tx_t,*gb_tx_t_n;
//extern real *gb_ty_t,*gb_ty_t_n;
//extern real *gb_xx_t,*gb_xx_t_n;
//extern real *gb_xy_t,*gb_xy_t_n;
//extern real *gb_yy_t,*gb_yy_t_n;
//extern real *gb_zz_t,*gb_zz_t_n;
//
//extern real *Hb_t_t,*Hb_t_t_n;
//extern real *Hb_x_t,*Hb_x_t_n;
//extern real *Hb_y_t,*Hb_y_t_n;
//
//extern real *mask,*mask_mg,*chr,*chr_mg,AMRD_ex;
//
//extern real *efe_all_ires;
//extern real *efe_tt_ires,*efe_tx_ires,*efe_ty_ires;
//extern real *efe_xx_ires,*efe_xy_ires,*efe_yy_ires;
//extern real *efe_zz_ires;
//
//extern real *cl_res;

/* prototypes for the various fortran functions we use */

void init_gbhb_(real *gb_tt, real *gb_tx, real *gb_ty,
                real *gb_xx, real *gb_xy, real *gb_yy, 
                real *gb_zz, 
                real *Hb_t, real *Hb_x, real *Hb_y, 
                real *AdS_L, real *x, real *y, real *chr, real *ex, int *Nx, int *Ny);

void lin_zero_bnd_res_(real *f, int *phys_bdy, int *all, int *Nx, int *Ny,const int *app_dim);

void hb_t_evo_(real *res,
               real *gb_tt_np1, real *gb_tt_n, real *gb_tt_nm1,
               real *gb_tx_np1, real *gb_tx_n, real *gb_tx_nm1,
               real *gb_ty_np1, real *gb_ty_n, real *gb_ty_nm1,
               real *gb_xx_np1, real *gb_xx_n, real *gb_xx_nm1,
               real *gb_xy_np1, real *gb_xy_n, real *gb_xy_nm1,
               real *gb_yy_np1, real *gb_yy_n, real *gb_yy_nm1,
               real *gb_zz_np1, real *gb_zz_n, real *gb_zz_nm1,
               real *Hb_t_np1, real *Hb_t_n, real *Hb_t_nm1,
               real *Hb_x_np1, real *Hb_x_n, real *Hb_x_nm1,
               real *Hb_y_np1, real *Hb_y_n, real *Hb_y_nm1,
               real *AdS_L, real *x, real *y, real *dt, real *chr, real *ex, 
               int *phys_bdy, int *ghost_width, int *Nx, int *Ny, const int *app_dim,
               real *Hb_t_0, real *Hb_x_0, real *Hb_y_0,
               int *gauge, real *t_n, 
               real *x1_t, real *x2_t, real *x3_t, real *x4_t, real *y3_t, real *y4_t, real *xi1_t, real *xi2_t, 
               real *cbulk, int *unibs_background);

void hb_i_evo_(real *res,
               real *gb_tt_np1, real *gb_tt_n, real *gb_tt_nm1,
               real *gb_tx_np1, real *gb_tx_n, real *gb_tx_nm1,
               real *gb_ty_np1, real *gb_ty_n, real *gb_ty_nm1,
               real *gb_xx_np1, real *gb_xx_n, real *gb_xx_nm1,
               real *gb_xy_np1, real *gb_xy_n, real *gb_xy_nm1,
               real *gb_yy_np1, real *gb_yy_n, real *gb_yy_nm1,
               real *gb_zz_np1, real *gb_zz_n, real *gb_zz_nm1,
               real *Hb_t_np1, real *Hb_t_n, real *Hb_t_nm1,
               real *Hb_x_np1, real *Hb_x_n, real *Hb_x_nm1,
               real *Hb_y_np1, real *Hb_y_n, real *Hb_y_nm1,
               real *AdS_L, real *x, real *y, real *dt, real *chr, real *ex, 
               int *phys_bdy, int *ghost_width, int *Nx, int *Ny, const int *app_dim,
               real *Hb_t_0, real *Hb_x_0, real *Hb_y_0,
               int *gauge, real *t_n, 
               real *x1_t, real *x2_t, real *x3_t, real *x4_t, real *y3_t, real *y4_t, real *xi1_t, real *xi2_t, 
               real *cbulk, int *unibs_background);

void UniBSAdS5D_fill_ex_mask(real *mask, int dim, int *shape, real *bbox, real excised);

#endif
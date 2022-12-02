#ifndef _UNIBSADS5D_H
#define _UNIBSADS5D_H

/*============================================================================= */
/* global variables and prototypes for UniBSAdS5D                               */
/* "extern" is short for external. It means that the declared variable          */
/* is defined (i.e., memory is allocated)                                       */
/* outside of the current file (in UniBSAdS5D.c in this code).                  */
/* It is needed for variables that are used in the declarations of functions in */
/* this file.                                                                   */
/*============================================================================= */


#define UNIBSADS5D_CP_VERSION 1
extern int cp_version;

extern int gb_tt_gfn,gb_tt_n_gfn,gb_tt_np1_gfn,gb_tt_nm1_gfn; 
extern int gb_tx_gfn,gb_tx_n_gfn,gb_tx_np1_gfn,gb_tx_nm1_gfn; 
extern int gb_ty_gfn,gb_ty_n_gfn,gb_ty_np1_gfn,gb_ty_nm1_gfn; 
extern int gb_xx_gfn,gb_xx_n_gfn,gb_xx_np1_gfn,gb_xx_nm1_gfn; 
extern int gb_xy_gfn,gb_xy_n_gfn,gb_xy_np1_gfn,gb_xy_nm1_gfn; 
extern int gb_yy_gfn,gb_yy_n_gfn,gb_yy_np1_gfn,gb_yy_nm1_gfn; 
extern int gb_zz_gfn,gb_zz_n_gfn,gb_zz_np1_gfn,gb_zz_nm1_gfn; 

extern int Hb_t_gfn,Hb_t_n_gfn,Hb_t_np1_gfn,Hb_t_nm1_gfn;
extern int Hb_x_gfn,Hb_x_n_gfn,Hb_x_np1_gfn,Hb_x_nm1_gfn;
extern int Hb_y_gfn,Hb_y_n_gfn,Hb_y_np1_gfn,Hb_y_nm1_gfn;

extern real *gb_tt_t,*gb_tt_t_n;
extern real *gb_tx_t,*gb_tx_t_n;
extern real *gb_ty_t,*gb_ty_t_n;
extern real *gb_xx_t,*gb_xx_t_n;
extern real *gb_xy_t,*gb_xy_t_n;
extern real *gb_yy_t,*gb_yy_t_n;
extern real *gb_zz_t,*gb_zz_t_n;

extern real *Hb_t_t,*Hb_t_t_n;
extern real *Hb_x_t,*Hb_x_t_n;
extern real *Hb_y_t,*Hb_y_t_n;

extern real *mask,*mask_mg,*chr,*chr_mg,AMRD_ex;

extern real *efe_all_ires;
extern real *efe_tt_ires,*efe_tx_ires,*efe_ty_ires;
extern real *efe_xx_ires,*efe_xy_ires,*efe_yy_ires;
extern real *efe_zz_ires;

extern real *cl_res;


#endif
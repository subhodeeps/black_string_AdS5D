#ifndef _APPH_H
#define _APPH_H

//-----------------------------------------------------------------------------
// global variables and prototypes for AH-finder
//-----------------------------------------------------------------------------

extern int AH_Nchi[MAX_BHS],AH_Nphi[MAX_BHS],AH_Lmin[MAX_BHS],AH_Lmax[MAX_BHS],AH_find_best_fit[MAX_BHS];
extern int AH_max_iter[MAX_BHS],AH_freq[MAX_BHS],AH_freq_aft[MAX_BHS],AH_rsteps[MAX_BHS],AH_maxinc[MAX_BHS];
extern real AH_tol[MAX_BHS],AH_tol_aft[MAX_BHS],AH_r0[MAX_BHS],AH_r1[MAX_BHS],AH_lambda[MAX_BHS],AH_lambda_min[MAX_BHS];
extern real AH_eps[MAX_BHS],AH_tol_scale[MAX_BHS],AH_reset_scale[MAX_BHS];
extern real AH_xc[MAX_BHS][3],AH_max_tol_inc[MAX_BHS],AH_tmin[MAX_BHS],AH_omt_scale[MAX_BHS];
extern int use_AH_new_smooth,use_AH_new;
extern int c_AH;

extern real *AH_theta[MAX_BHS],*AH_R[MAX_BHS],*AH_w1[MAX_BHS],*AH_w2[MAX_BHS],*AH_w3[MAX_BHS],*AH_w4[MAX_BHS];
extern real *AH_x0[MAX_BHS],*AH_y0[MAX_BHS],*AH_z0[MAX_BHS];
extern real *AH_g0_tt[MAX_BHS];
extern real *AH_g0_tx[MAX_BHS];
extern real *AH_g0_ty[MAX_BHS];
extern real *AH_g0_tz[MAX_BHS];
extern real *AH_g0_xx[MAX_BHS];
extern real *AH_g0_xy[MAX_BHS];
extern real *AH_g0_xz[MAX_BHS];
extern real *AH_g0_yy[MAX_BHS];
extern real *AH_g0_yz[MAX_BHS];
extern real *AH_g0_zz[MAX_BHS];
extern real *AH_g0_chichi[MAX_BHS];
extern real *AH_g0_chiphi[MAX_BHS];
extern real *AH_g0_phiphi[MAX_BHS];
extern real *AH_kretsch[MAX_BHS];
extern real *AH_riemanncube[MAX_BHS];
extern real *AH_ahr[MAX_BHS],*AH_dph[MAX_BHS],*AH_dch[MAX_BHS];
extern real *AH_da0[MAX_BHS],*AH_dcq[MAX_BHS],*AH_dcp[MAX_BHS],*AH_dcp2[MAX_BHS];

extern int *AH_lev[MAX_BHS],*AH_own[MAX_BHS];

extern int AH_count[MAX_BHS],found_AH[MAX_BHS],freq0[MAX_BHS];
extern int found_count_AH[MAX_BHS];

//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------

int find_apph(real *M, real *J, real *area, real *c_equat, real *c_polar, real *c_polar2s, int use_R_ic, real *AH_min_resid, 
    real *ief_bh_r0,real *a_rot0,int *kerrads_background, int *AH_analytic_kerrads, real *ct);

//-----------------------------------------------------------------------------
// related fotran routines
//-----------------------------------------------------------------------------

void is_inside_(int *a_in_b, int *a_int_b, real *r_a, real *xc_a, real *r_b, real *xc_b, int *dim);

void ah_is_int_(int *is_int, real *AH_R, real *AH_xc, int *i, int *j, real *bbox, real *dx, real *dy, real *dz,
               int *AH_Nchi, int *AH_Nphi, int *axisym);

void ah_fill_own_(real *AH_R, real *AH_xc, int *AH_own, int *AH_lev, real *bbox, real *dx, real *dy, real *dz, 
                 int *rank, int *AdS_L, int *AH_Nchi, int *AH_Nphi, int *axisym);

void ah_fill_f_(real *AH_R, real *AH_xc, real *f, int *is, int *ie, int *js, int *je, int *ks, int *ke,
                real *x, real *y, real *z, int *AH_Nchi, int *AH_Nphi, int *Nx, int *Ny, int *Nz, int *axisym);

void calc_exp_metric0_(real *AH_R, real *AH_xc, real *AH_theta, int *i0, int *j0, int *AH_Nchi,int *AH_Nphi, 
                real *theta, real *f, real *da, real *d_ceq, real *d_cp, real *d_cp2,
                real *AH_x0, real *AH_y0, real *AH_z0,
                real *AH_g0_xx,real *AH_g0_xy,real *AH_g0_xz,
                real *AH_g0_yy,real *AH_g0_yz,real *AH_g0_zz,
                real *AH_g0_chichi, real *AH_g0_chiphi, real *AH_g0_phiphi, 
                real *AH_kretsch,
                real *AH_riemanncube,
                real *AH_ahr,real *AH_dch,real *AH_dph,
                real *AH_da0, real *AH_dcq, real *AH_dcp, real *AH_dcp2,
                real *gb_tt_np1, real *gb_tt_n, real *gb_tt_nm1, 
                real *gb_tx_np1, real *gb_tx_n, real *gb_tx_nm1, 
                real *gb_ty_np1, real *gb_ty_n, real *gb_ty_nm1,
                real *gb_tz_np1, real *gb_tz_n, real *gb_tz_nm1,
                real *gb_xx_np1, real *gb_xx_n, real *gb_xx_nm1,
                real *gb_xy_np1, real *gb_xy_n, real *gb_xy_nm1,
                real *gb_xz_np1, real *gb_xz_n, real *gb_xz_nm1,
                real *gb_yy_np1, real *gb_yy_n, real *gb_yy_nm1,
                real *gb_yz_np1, real *gb_yz_n, real *gb_yz_nm1,
                real *gb_zz_np1, real *gb_zz_n, real *gb_zz_nm1,
                real *Hb_t_np1, real *Hb_t_n, real *Hb_t_nm1,
                real *Hb_x_np1, real *Hb_x_n, real *Hb_x_nm1,
                real *Hb_y_np1, real *Hb_y_n, real *Hb_y_nm1,
                real *Hb_z_np1, real *Hb_z_n, real *Hb_z_nm1,
                real *phi1_np1, real *phi1_n, real *phi1_nm1,
                real *kretsch_n,
                real *riemanncube_n,
                real *AdS_L, real *x, real *y, real *z, real *dt, real *chr, 
                real *ex, int *do_ex, int *Nx, int *Ny, int *Nz, int *axisym,
                real *ief_bh_r0,real *a_rot0,int *kerrads_background);

void smooth_ah_r_(real *AH_R,real *AH_w1,real *AH_eps, int *AH_Nchi, int *AH_Nphi);

void smooth_ah_r_b_(real *AH_R,real *AH_w1,real *AH_eps, int *AH_Nchi, int *AH_Nphi);

void reg_ah_r_(real *AH_R, int *AH_Nchi, int *AH_Nphi);

void adjust_ah_xc_(real *AH_R, real *AH_xc, int *AH_Nchi, int *AH_Nphi, real *dx, real *dy, real *dz, int *axisym);

void fill_ex_params_(real *AH_R, real *AH_xc, real *min_AH_R, real *max_AH_R, real *ex_r0, real *ex_xc, int *AH_Nchi, int *AH_Nphi, real *dx, real *dy, real *dz, int *axisym);



#endif

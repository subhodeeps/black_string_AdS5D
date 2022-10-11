c----------------------------------------------------------------------
c in cartesian coordinates t,x,y,z for x,y,z in [-1,1]
c
c An experimental evolution routine for the gb,phi1, computing
c the residual at the time just prior to updated
c
c choosing theta=pi/2
c
c L below is the AdS length scale
c----------------------------------------------------------------------
        subroutine g_evo_opt(gb_res,kg_res,cl_res,
     &                       gb_tt_np1,gb_tt_n,gb_tt_nm1,
     &                       gb_tx_np1,gb_tx_n,gb_tx_nm1,
     &                       gb_ty_np1,gb_ty_n,gb_ty_nm1,
     &                       gb_tz_np1,gb_tz_n,gb_tz_nm1,
     &                       gb_xx_np1,gb_xx_n,gb_xx_nm1,
     &                       gb_xy_np1,gb_xy_n,gb_xy_nm1,     
     &                       gb_xz_np1,gb_xz_n,gb_xz_nm1,
     &                       gb_yy_np1,gb_yy_n,gb_yy_nm1,
     &                       gb_yz_np1,gb_yz_n,gb_yz_nm1,
     &                       gb_zz_np1,gb_zz_n,gb_zz_nm1,
     &                       Hb_t_np1,Hb_t_n,Hb_t_nm1,
     &                       Hb_x_np1,Hb_x_n,Hb_x_nm1,
     &                       Hb_y_np1,Hb_y_n,Hb_y_nm1,
     &                       Hb_z_np1,Hb_z_n,Hb_z_nm1,
     &                       phi1_np1,phi1_n,phi1_nm1,
     &                       L,x,y,z,dt,chr,ex,
     &                       phys_bdy,ghost_width,Nx,Ny,Nz,
     &                       background,kappa_cd,rho_cd,
     &                       interptype,i_shift,regtype,
     &                       diss_kmax,tfunction,
     &                       ief_bh_r0,a_rot,kerrads_background)
        implicit none
        real*8 ief_bh_r0,a_rot
        integer kerrads_background
        logical calc_der,calc_adv_quant
        data calc_der/.true./
        data calc_adv_quant/.false./
        integer Nx,Ny,Nz
        integer phys_bdy(6),ghost_width(6)
        integer background
        integer interptype
        integer regtype
        integer i_shift
        integer diss_kmax
        integer max_ghost_width
        real*8 kappa_cd,rho_cd
        real*8 gb_res(Nx,Ny,Nz),kg_res(Nx,Ny,Nz),cl_res(Nx,Ny,Nz)
        real*8 gb_tt_np1(Nx,Ny,Nz),gb_tx_np1(Nx,Ny,Nz)
        real*8 gb_ty_np1(Nx,Ny,Nz)
        real*8 gb_tz_np1(Nx,Ny,Nz)
        real*8 gb_xx_np1(Nx,Ny,Nz),gb_xy_np1(Nx,Ny,Nz)
        real*8 gb_xz_np1(Nx,Ny,Nz)
        real*8 gb_yz_np1(Nx,Ny,Nz)
        real*8 gb_yy_np1(Nx,Ny,Nz),gb_zz_np1(Nx,Ny,Nz)
        real*8 gb_tt_n(Nx,Ny,Nz),gb_tx_n(Nx,Ny,Nz)
        real*8 gb_ty_n(Nx,Ny,Nz)
        real*8 gb_tz_n(Nx,Ny,Nz)
        real*8 gb_xx_n(Nx,Ny,Nz),gb_xy_n(Nx,Ny,Nz)
        real*8 gb_xz_n(Nx,Ny,Nz)
        real*8 gb_yz_n(Nx,Ny,Nz)
        real*8 gb_yy_n(Nx,Ny,Nz),gb_zz_n(Nx,Ny,Nz)
        real*8 gb_tt_nm1(Nx,Ny,Nz),gb_tx_nm1(Nx,Ny,Nz)
        real*8 gb_ty_nm1(Nx,Ny,Nz)
        real*8 gb_tz_nm1(Nx,Ny,Nz)
        real*8 gb_xx_nm1(Nx,Ny,Nz),gb_xy_nm1(Nx,Ny,Nz)
        real*8 gb_xz_nm1(Nx,Ny,Nz)
        real*8 gb_yz_nm1(Nx,Ny,Nz)
        real*8 gb_yy_nm1(Nx,Ny,Nz),gb_zz_nm1(Nx,Ny,Nz)
        real*8 Hb_t_n(Nx,Ny,Nz),Hb_x_n(Nx,Ny,Nz)
        real*8 Hb_y_n(Nx,Ny,Nz),Hb_z_n(Nx,Ny,Nz)
        real*8 Hb_t_np1(Nx,Ny,Nz),Hb_x_np1(Nx,Ny,Nz)
        real*8 Hb_y_np1(Nx,Ny,Nz),Hb_z_np1(Nx,Ny,Nz)
        real*8 Hb_t_nm1(Nx,Ny,Nz),Hb_x_nm1(Nx,Ny,Nz)
        real*8 Hb_y_nm1(Nx,Ny,Nz),Hb_z_nm1(Nx,Ny,Nz)
        real*8 phi1_np1(Nx,Ny,Nz),phi1_n(Nx,Ny,Nz),phi1_nm1(Nx,Ny,Nz)
        real*8 tfunction(Nx,Ny,Nz)
        real*8 L
        real*8 x(Nx),y(Ny),z(Nz),dt,chr(Nx,Ny,Nz),ex
        real*8 chr2(Nx,Ny,Nz)

        integer a,b,c,d,e
        integer rb,i,j,k,m
        integer is,ie,js,je,ks,ke,is_a_nan

        real*8 dx,dy,dz
        real*8 x0,y0,z0,rho0

        real*8 phi1_res,phi1_J

        real*8 PI
        parameter (PI=3.141592653589793d0)

        !--------------------------------------------------------------
        ! the following are first and second time derivatives of *n*
        ! level variables, and as these are the only derivatives we
        ! use we drop any _n identifier
        !--------------------------------------------------------------
        real*8 g0_tt_t, g0_tt_x, g0_tt_y, g0_tt_tt
        real*8 g0_tt_xx,g0_tt_yy,g0_tt_tx,g0_tt_ty
        real*8 g0_tt_xy
        real*8 g0_tt_z
        real*8 g0_tt_tz,g0_tt_xz,g0_tt_yz,g0_tt_zz

        real*8 g0_tx_t, g0_tx_x, g0_tx_y, g0_tx_tt
        real*8 g0_tx_xx,g0_tx_yy,g0_tx_tx,g0_tx_ty
        real*8 g0_tx_xy
        real*8 g0_tx_z
        real*8 g0_tx_tz,g0_tx_xz,g0_tx_yz,g0_tx_zz

        real*8 g0_ty_t, g0_ty_x, g0_ty_y, g0_ty_tt
        real*8 g0_ty_xx,g0_ty_yy,g0_ty_tx,g0_ty_ty
        real*8 g0_ty_xy
        real*8 g0_ty_z
        real*8 g0_ty_tz,g0_ty_xz,g0_ty_yz,g0_ty_zz

        real*8 g0_tz_t, g0_tz_x, g0_tz_y, g0_tz_tt
        real*8 g0_tz_xx,g0_tz_yy,g0_tz_tx,g0_tz_ty
        real*8 g0_tz_xy
        real*8 g0_tz_z
        real*8 g0_tz_tz,g0_tz_xz,g0_tz_yz,g0_tz_zz

        real*8 g0_xx_t, g0_xx_x, g0_xx_y, g0_xx_tt
        real*8 g0_xx_xx,g0_xx_yy,g0_xx_tx,g0_xx_ty
        real*8 g0_xx_xy
        real*8 g0_xx_z
        real*8 g0_xx_tz,g0_xx_xz,g0_xx_yz,g0_xx_zz

        real*8 g0_xy_t, g0_xy_x, g0_xy_y, g0_xy_tt
        real*8 g0_xy_xx,g0_xy_yy,g0_xy_tx,g0_xy_ty
        real*8 g0_xy_xy
        real*8 g0_xy_z
        real*8 g0_xy_tz,g0_xy_xz,g0_xy_yz,g0_xy_zz

        real*8 g0_xz_t, g0_xz_x, g0_xz_y, g0_xz_tt
        real*8 g0_xz_xx,g0_xz_yy,g0_xz_tx,g0_xz_ty
        real*8 g0_xz_xy
        real*8 g0_xz_z
        real*8 g0_xz_tz,g0_xz_xz,g0_xz_yz,g0_xz_zz

        real*8 g0_yy_t, g0_yy_x, g0_yy_y, g0_yy_tt
        real*8 g0_yy_xx,g0_yy_yy,g0_yy_tx,g0_yy_ty
        real*8 g0_yy_xy
        real*8 g0_yy_z
        real*8 g0_yy_tz,g0_yy_xz,g0_yy_yz,g0_yy_zz

        real*8 g0_yz_t, g0_yz_x, g0_yz_y, g0_yz_tt
        real*8 g0_yz_xx,g0_yz_yy,g0_yz_tx,g0_yz_ty
        real*8 g0_yz_xy
        real*8 g0_yz_z
        real*8 g0_yz_tz,g0_yz_xz,g0_yz_yz,g0_yz_zz

        real*8 g0_zz_t, g0_zz_x, g0_zz_y, g0_zz_tt
        real*8 g0_zz_xx,g0_zz_yy,g0_zz_tx,g0_zz_ty
        real*8 g0_zz_xy
        real*8 g0_zz_z
        real*8 g0_zz_tz,g0_zz_xz,g0_zz_yz,g0_zz_zz

        real*8 g0u_tt0,g0_tt0
        real*8 g0u_tx0,g0_tx0
        real*8 g0u_ty0,g0_ty0
        real*8 g0u_tz0,g0_tz0
        real*8 g0u_xx0,g0_xx0
        real*8 g0u_xy0,g0_xy0
        real*8 g0u_xz0,g0_xz0
        real*8 g0u_yy0,g0_yy0
        real*8 g0u_yz0,g0_yz0
        real*8 g0u_zz0,g0_zz0

        real*8 H0_t0,H0_x0,H0_y0
        real*8 H0_z0

        real*8 C_t,C_x,C_y
        real*8 C_t_tt_J,C_t_tx_J,C_t_ty_J,C_t_xx_J
        real*8 C_t_xy_J,C_t_yy_J,C_t_zz_J
        real*8 C_x_tt_J,C_x_tx_J,C_x_ty_J,C_x_xx_J
        real*8 C_x_xy_J,C_x_yy_J,C_x_zz_J
        real*8 C_y_tt_J,C_y_tx_J,C_y_ty_J,C_y_xx_J
        real*8 C_y_xy_J,C_y_yy_J,C_y_zz_J
        real*8 nu_t,nu_x,nu_y,nl_t,nl_x,nl_y
        real*8 d_gb_tt_res,d_gb_tx_res,d_gb_ty_res
        real*8 d_gb_xx_res,d_gb_xy_res,d_gb_yy_res
        real*8 d_zz_res
        real*8 d_gb_tt_J,d_gb_tx_J,d_gb_ty_J
        real*8 d_gb_tz_J
        real*8 d_gb_xx_J,d_gb_xy_J,d_gb_yy_J
        real*8 d_zz_J

        logical ltrace,is_nan,dump,first_nan
        logical first_evolved_pt
        logical extrap
        data extrap/.true./
        parameter (ltrace=.false.)
        data first_nan/.true./

        integer i2,j2

        !--------------------------------------------------------------
        ! using g0_ab = gads_ab + h0_ab
        !       g0^ab = gads^ab + h0^ab
        !       H0_a  = Hads_a + A_a
        ! where h0_ab, h0^ab are not inverses of each other
        !
        ! g0_ll(a,b)           = g0_ab           !use for diagnostics  !
        ! g0_uu(a,b)           = g0^ab           !eg: compare g0_ll    !
        ! g0_ll_x(a,b,c)       = g0_ab_,c        !    with gads_ll+h_ll!  
        ! g0_ll_xx(a,b,c,d)    = g0_ab_,cd
        ! g0_uu_x(a,b,c)       = g0^ab_,c
        !                      = -g0^ad g0^be g0_de_,c 
        !
        ! gads_ll(a,b)         = gads_ab
        ! gads_uu(a,b)         = gads^ab
        ! gads_ll_x(a,b,c)     = gads_ab_,c
        ! gads_ll_xx(a,b,c,d)  = gads_ab_,cd
        ! gads_uu_x(a,b,c)     = gads^ab_,c
        !                      = -gads^ad gads^be gads_de_,c 
        ! 
        ! h0_ll(a,b)           = gb_ab
        ! h0_uu(a,b)           = inverse(gads+gb)^ab - gads^ab
        ! h0_ll_x(a,b,c)       = (gb_ab)_,c
        ! h0_ll_xx(a,b,c,d)    = (gb_ab)_,cd
        ! h0_uu_x(a,b,c)       = g^ab_,c - gads^ab_,c
        !
        ! gammagg(a,b,c) = 0.5d0*gads_uu(a,d) 
        !                  *(gads_ll_x(c,d,b)-gads_ll_x(b,c,d)+gads_ll_x(d,b,c))
        ! gammahh(a,b,c) = 0.5d0*h0_uu(a,d) 
        !                  *(h0_ll_x(c,d,b)-h0_ll_x(b,c,d)+h0_ll_x(d,b,c))
        ! gammagh(a,b,c) = 0.5d0*gads_uu(a,d) 
        !                  *(h0_ll_x(c,d,b)-h0_ll_x(b,c,d)+h0_ll_x(d,b,c))
        ! gammahg(a,b,c) = 0.5d0*h0_uu(a,d) 
        !                  *(gads_ll_x(c,d,b)-gads_ll_x(b,c,d)+gads_ll_x(d,b,c))
        !
        ! cuuuu(a,b,c,d) = gads_uu(a,b)*gads_uu(c,d)
        !                   +h0_uu(a,b)*h0_uu(c,d) 
        !                   +gads_uu(a,b)*h0_uu(c,d) 
        !                   +h0_uu(a,b)*gads_uu(c,d) 
        !
        ! dlll(a,b,c)    = g0_ll_x(b,c,a)-g0_ll_x(a,b,c)+g0_ll_x(c,a,b)
        !
        ! A_l(a)     = Hb_a
        ! A_l_x(a,b) = Hb_a,b
        ! 
        ! phi10_x(a)= phi1_,a
        !
        ! grad_phi1_sq = g^cd*phi1_,c*phi1_,d
        !
        ! set_ab = 2*phi1_,a*phi1_,b - g_ab*grad_phi1_sq
        ! tr_set= g^cd*se_cd 
        !
        ! efe(a,b) = residual ... hardcoded expressions (see below)
        !
        ! t,x,y=1,2,3
        ! 
        ! NOTE: g0_ll_xx,gads_ll_xx,h0_ll_xx,efe,efe_J
        !       do *NOT* symmetric components filled in
        !
        !--------------------------------------------------------------
        real*8 efe(4,4),efe_J(4,4)
        real*8 term1(4,4),term2(4,4),term3(4,4),term4(4,4)
        real*8 term5(4,4),term6(4,4),term7(4,4),term8(4,4)
        real*8 gammagg(4,4,4),gammahh(4,4,4)
        real*8 gammagh(4,4,4),gammahg(4,4,4) 
        real*8 cuuuu(4,4,4,4),dlll(4,4,4)
        real*8 dphi1(4)

        real*8 ndotc,n_l(4),n_u(4),c_l(4),c_J_l(4)
        real*8 cd_ll(4,4),cd_J_ll(4,4)
 
        real*8 tr_set,grad_phi1_sq
        
        real*8 g0u_tt_ads0,g0u_xx_ads0,g0u_xy_ads0,g0u_yy_ads0
        real*8 g0u_zz_ads0

        real*8 H0_t_ads0,H0_x_ads0,H0_y_ads0

        real*8 dgb_J,ddgb_J,ddgb_J_tx,ddgb_J_ty
        real*8 ddgb_J_tz
        real*8 dc_J

        real*8 lambda4

        real*8 h0_tt0
        real*8 h0_tx0
        real*8 h0_ty0
        real*8 h0_xx0
        real*8 h0_xy0
        real*8 h0_yy0
        real*8 h0_zz0

        real*8 h0_tt_t, h0_tt_x, h0_tt_y, h0_tt_tt
        real*8 h0_tt_xx,h0_tt_yy,h0_tt_tx,h0_tt_ty
        real*8 h0_tt_xy

        real*8 h0_tx_t, h0_tx_x, h0_tx_y, h0_tx_tt
        real*8 h0_tx_xx,h0_tx_yy,h0_tx_tx,h0_tx_ty
        real*8 h0_tx_xy

        real*8 h0_ty_t, h0_ty_x, h0_ty_y, h0_ty_tt
        real*8 h0_ty_xx,h0_ty_yy,h0_ty_tx,h0_ty_ty
        real*8 h0_ty_xy

        real*8 h0_tz_t, h0_tz_x, h0_tz_y, h0_tz_tt
        real*8 h0_tz_xx,h0_tz_yy,h0_tz_tx,h0_tz_ty
        real*8 h0_tz_xy

        real*8 h0_xx_t, h0_xx_x, h0_xx_y, h0_xx_tt
        real*8 h0_xx_xx,h0_xx_yy,h0_xx_tx,h0_xx_ty
        real*8 h0_xx_xy

        real*8 h0_xy_t, h0_xy_x, h0_xy_y, h0_xy_tt
        real*8 h0_xy_xx,h0_xy_yy,h0_xy_tx,h0_xy_ty
        real*8 h0_xy_xy

        real*8 h0_xz_t, h0_xz_x, h0_xz_y, h0_xz_tt
        real*8 h0_xz_xx,h0_xz_yy,h0_xz_tx,h0_xz_ty
        real*8 h0_xz_xy

        real*8 h0_yy_t, h0_yy_x, h0_yy_y, h0_yy_tt
        real*8 h0_yy_xx,h0_yy_yy,h0_yy_tx,h0_yy_ty
        real*8 h0_yy_xy

        real*8 h0_yz_t, h0_yz_x, h0_yz_y, h0_yz_tt
        real*8 h0_yz_xx,h0_yz_yy,h0_yz_tx,h0_yz_ty
        real*8 h0_yz_xy

        real*8 h0_zz_t, h0_zz_x, h0_zz_y, h0_zz_tt
        real*8 h0_zz_xx,h0_zz_yy,h0_zz_tx,h0_zz_ty
        real*8 h0_zz_xy

        !--------------------------------------------------------------
        ! variables for tensor manipulations 
        !(indices are t,x,w,y,z)
        !--------------------------------------------------------------
        real*8 g0_ll(4,4),g0_uu(4,4)
        real*8 g0_ll_x(4,4,4),g0_uu_x(4,4,4),g0_ll_xx(4,4,4,4)
        real*8 gads_ll(4,4),gads_uu(4,4)
        real*8 gads_ll_x(4,4,4),gads_uu_x(4,4,4),gads_ll_xx(4,4,4,4)
        real*8 h0_ll(4,4),h0_uu(4,4)
        real*8 h0_ll_x(4,4,4),h0_uu_x(4,4,4),h0_ll_xx(4,4,4,4)
        real*8 gamma_ull(4,4,4),gamma_ull_x(4,4,4,4)
        real*8 riemann_ulll(4,4,4,4)
        real*8 ricci_ll(4,4),ricci_lu(4,4),ricci
        real*8 einstein_ll(4,4),set_ll(4,4)
        real*8 Hads_l(4),A_l(4),A_l_x(4,4)
        real*8 phi10_x(4),phi10_xx(4,4)

        !--------------------------------------------------------------
        ! initialize fixed-size variables 
        !--------------------------------------------------------------
        data lambda4/0.0/

        data g0u_tt_ads0,g0u_xx_ads0/0.0,0.0/
        data g0u_xy_ads0,g0u_yy_ads0/0.0,0.0/
        data g0u_zz_ads0/0.0/

        data H0_t_ads0,H0_x_ads0,H0_y_ads0/0.0,0.0,0.0/

        data dgb_J,ddgb_J,ddgb_J_tx,ddgb_J_ty/0.0,0.0,0.0,0.0/
        data ddgb_J_tz/0.0/
        data dc_J/0.0/

        data dlll/64*0.0/
        data cuuuu/256*0.0/

        data term1,term2/16*0.0,16*0.0/
        data term3,term4/16*0.0,16*0.0/
        data term5,term6/16*0.0,16*0.0/
        data term7,term8/16*0.0,16*0.0/

        data efe,efe_J/16*0.0,16*0.0/
        data cd_ll,cd_J_ll/16*0.0,16*0.0/

        data n_l,n_u,c_l,c_J_l/4*0.0,4*0.0,4*0.0,4*0.0/



        data rb,i,j,k/0,0,0,0/
        data i2,j2/0,0/
        data is,ie,js,je,ks,ke,is_a_nan/0,0,0,0,0,0,0/
        data a,b,c,d,e/0,0,0,0,0/

        data dx,dy,dz/0.0,0.0,0.0/

        data g0_tt_t, g0_tt_x, g0_tt_y, g0_tt_tt/0.0,0.0,0.0,0.0/
        data g0_tt_xx,g0_tt_yy,g0_tt_tx,g0_tt_ty/0.0,0.0,0.0,0.0/
        data g0_tt_xy/0.0/
        data g0_tt_z/0.0/
        data g0_tt_tz,g0_tt_xz,g0_tt_yz,g0_tt_zz/0.0,0.0,0.0,0.0/

        data g0_tx_t, g0_tx_x, g0_tx_y, g0_tx_tt/0.0,0.0,0.0,0.0/
        data g0_tx_xx,g0_tx_yy,g0_tx_tx,g0_tx_ty/0.0,0.0,0.0,0.0/
        data g0_tx_xy/0.0/
        data g0_tx_z/0.0/
        data g0_tx_tz,g0_tx_xz,g0_tx_yz,g0_tx_zz/0.0,0.0,0.0,0.0/

        data g0_ty_t, g0_ty_x, g0_ty_y, g0_ty_tt/0.0,0.0,0.0,0.0/
        data g0_ty_xx,g0_ty_yy,g0_ty_tx,g0_ty_ty/0.0,0.0,0.0,0.0/
        data g0_ty_xy/0.0/
        data g0_ty_z/0.0/
        data g0_ty_tz,g0_ty_xz,g0_ty_yz,g0_ty_zz/0.0,0.0,0.0,0.0/

        data g0_tz_t, g0_tz_x, g0_tz_y, g0_tz_tt/0.0,0.0,0.0,0.0/
        data g0_tz_xx,g0_tz_yy,g0_tz_tx,g0_tz_ty/0.0,0.0,0.0,0.0/
        data g0_tz_xy/0.0/
        data g0_tz_z/0.0/
        data g0_tz_tz,g0_tz_xz,g0_tz_yz,g0_tz_zz/0.0,0.0,0.0,0.0/

        data g0_xx_t, g0_xx_x, g0_xx_y, g0_xx_tt/0.0,0.0,0.0,0.0/
        data g0_xx_xx,g0_xx_yy,g0_xx_tx,g0_xx_ty/0.0,0.0,0.0,0.0/
        data g0_xx_xy/0.0/
        data g0_xx_z/0.0/
        data g0_xx_tz,g0_xx_xz,g0_xx_yz,g0_xx_zz/0.0,0.0,0.0,0.0/

        data g0_xy_t, g0_xy_x, g0_xy_y, g0_xy_tt/0.0,0.0,0.0,0.0/
        data g0_xy_xx,g0_xy_yy,g0_xy_tx,g0_xy_ty/0.0,0.0,0.0,0.0/
        data g0_xy_xy/0.0/
        data g0_xy_z/0.0/
        data g0_xy_tz,g0_xy_xz,g0_xy_yz,g0_xy_zz/0.0,0.0,0.0,0.0/

        data g0_xz_t, g0_xz_x, g0_xz_y, g0_xz_tt/0.0,0.0,0.0,0.0/
        data g0_xz_xx,g0_xz_yy,g0_xz_tx,g0_xz_ty/0.0,0.0,0.0,0.0/
        data g0_xz_xy/0.0/
        data g0_xz_z/0.0/
        data g0_xz_tz,g0_xz_xz,g0_xz_yz,g0_xz_zz/0.0,0.0,0.0,0.0/

        data g0_yy_t, g0_yy_x, g0_yy_y, g0_yy_tt/0.0,0.0,0.0,0.0/
        data g0_yy_xx,g0_yy_yy,g0_yy_tx,g0_yy_ty/0.0,0.0,0.0,0.0/
        data g0_yy_xy/0.0/
        data g0_yy_z/0.0/
        data g0_yy_tz,g0_yy_xz,g0_yy_yz,g0_yy_zz/0.0,0.0,0.0,0.0/

        data g0_yz_t, g0_yz_x, g0_yz_y, g0_yz_tt/0.0,0.0,0.0,0.0/
        data g0_yz_xx,g0_yz_yy,g0_yz_tx,g0_yz_ty/0.0,0.0,0.0,0.0/
        data g0_yz_xy/0.0/
        data g0_yz_z/0.0/
        data g0_yz_tz,g0_yz_xz,g0_yz_yz,g0_yz_zz/0.0,0.0,0.0,0.0/

        data g0_zz_t, g0_zz_x, g0_zz_y, g0_zz_tt/0.0,0.0,0.0,0.0/
        data g0_zz_xx,g0_zz_yy,g0_zz_tx,g0_zz_ty/0.0,0.0,0.0,0.0/
        data g0_zz_xy/0.0/
        data g0_zz_z/0.0/
        data g0_zz_tz,g0_zz_xz,g0_zz_yz,g0_zz_zz/0.0,0.0,0.0,0.0/

        data g0u_tt0,g0_tt0/0.0,0.0/
        data g0u_tx0,g0_tx0/0.0,0.0/
        data g0u_ty0,g0_ty0/0.0,0.0/
        data g0u_tz0,g0_tz0/0.0,0.0/
        data g0u_xx0,g0_xx0/0.0,0.0/
        data g0u_xy0,g0_xy0/0.0,0.0/
        data g0u_xz0,g0_xz0/0.0,0.0/
        data g0u_yy0,g0_yy0/0.0,0.0/
        data g0u_yz0,g0_yz0/0.0,0.0/
        data g0u_zz0,g0_zz0/0.0,0.0/

        data H0_t0,H0_x0,H0_y0/0.0,0.0,0.0/
        data H0_z0/0.0/

        data x0,y0,z0,rho0/0.0,0.0,0.0,0.0/

        data C_t,C_x,C_y/0.0,0.0,0.0/
        data C_t_tt_J,C_t_tx_J,C_t_ty_J,C_t_xx_J/0.0,0.0,0.0,0.0/
        data C_t_xy_J,C_t_yy_J,C_t_zz_J/0.0,0.0,0.0/
        data C_x_tt_J,C_x_tx_J,C_x_ty_J,C_x_xx_J/0.0,0.0,0.0,0.0/
        data C_x_xy_J,C_x_yy_J,C_x_zz_J/0.0,0.0,0.0/
        data C_y_tt_J,C_y_tx_J,C_y_ty_J,C_y_xx_J/0.0,0.0,0.0,0.0/
        data C_y_xy_J,C_y_yy_J,C_y_zz_J/0.0,0.0,0.0/
        data nu_t,nu_x,nu_y,nl_t,nl_x,nl_y/0.0,0.0,0.0,0.0,0.0,0.0/
        data d_gb_tt_res,d_gb_tx_res,d_gb_ty_res/0.0,0.0,0.0/
        data d_gb_xx_res,d_gb_xy_res,d_gb_yy_res/0.0,0.0,0.0/
        data d_zz_res/0.0/
        data d_gb_tt_J,d_gb_tx_J,d_gb_ty_J,d_gb_tz_J/0.0,0.0,0.0,0.0/
        data d_gb_xx_J,d_gb_xy_J,d_gb_yy_J/0.0,0.0,0.0/
        data d_zz_J/0.0/

        data grad_phi1_sq/1*0.0/
        data Hads_l,A_l,dphi1/4*0.0,4*0.0,4*0.0/
        data A_l_x/16*0.0/
        data g0_ll,g0_uu,gads_ll/16*0.0,16*0.0,16*0.0/
        data gads_uu,h0_ll,h0_uu/16*0.0,16*0.0,16*0.0/
        data gammagg,gammahh/64*0.0,64*0.0/
        data gammagh,gammahg/64*0.0,64*0.0/
        data g0_ll_x,g0_uu_x/64*0.0,64*0.0/
        data gads_ll_x,gads_uu_x/64*0.0,64*0.0/
        data h0_ll_x,h0_uu_x/64*0.0,64*0.0/
        data g0_ll_xx/256*0.0/
        data gads_ll_xx,h0_ll_xx/256*0.0,256*0.0/

        data gamma_ull/64*0.0/
        data gamma_ull_x/256*0.0/
        data riemann_ulll/256*0.0/
        data ricci/0.0/
        data ricci_ll,ricci_lu/16*0.0,16*0.0/
        data einstein_ll,set_ll/16*0.0,16*0.0/

        data phi10_x/4*0.0/
        data phi10_xx/16*0.0/

        !--------------------------------------------------------------
        if (ltrace) write(*,*) 'gb_zz_evo ... N=',Nx,Ny,Nz

        dx=x(2)-x(1)
        dy=y(2)-y(1)
        dz=z(2)-z(1)

        if (abs((y(2)-y(1))/dx-1).gt.1.0d-8) then
           write(*,*) 'error ... g_evo_opt not updated for dx!=dy!=dz'
           stop
        end if

        ! AdS4D cosmological constant
        !(lambda4=-(n-1)(n-2)/2/L^2) for n=4 dimensional AdS)
        lambda4=-3/L/L

        ! set index bounds for main loop
        is=2
        ie=Nx-1
        js=2
        je=Ny-1
        ks=2
        ke=Nz-1

        !(nearest-to-axis points are not evolved, according to regtype choice) 
!        if (regtype.eq.7 .or. regtype.eq.6) then
!          if (abs(y(1)).lt.dy/2) js=4
!        else if (regtype.eq.5 .or. regtype.eq.4 .or. regtype.eq.3) then
!          if (abs(y(1)).lt.dy/2) js=3
!        else
!          if (abs(y(1)).lt.dy/2) js=2
!        endif

        ! adjust index bounds to compensate for ghost_width
        if (ghost_width(1).gt.0) is=is+ghost_width(1)-1
        if (ghost_width(2).gt.0) ie=ie-(ghost_width(2)-1)
        if (ghost_width(3).gt.0) js=js+ghost_width(3)-1
        if (ghost_width(4).gt.0) je=je-(ghost_width(4)-1)
        if (ghost_width(5).gt.0) ks=ks+ghost_width(5)-1
        if (ghost_width(6).gt.0) ke=ke-(ghost_width(6)-1)

!        write(*,*) "Nx,Ny,Nz="
!     &             ,Nx,Ny,Nz
!        write(*,*) "ghost_width(1),ghost_width(2),ghost_width(3)
!     &             ,ghost_width(4),ghost_width(5),ghost_width(6)="
!     &             ,ghost_width(1),ghost_width(2),ghost_width(3)
!     &             ,ghost_width(4),ghost_width(5),ghost_width(6)

!       write(*,*) "is,ie,js,je,ks,ke,chr(is,js,ks-1)="
!     &             ,is,ie,js,je,ks,ke,chr(9,9,ks-1)

!        write(*,*) "2:gb_xx_n(is,js,ks)=",gb_xx_n(is,js,ks)
!        write(*,*) "3:gb_xx_n(is,js,ks-1)=",gb_xx_n(9,9,ks-1)

        ! check kmax value against ghost_width
        max_ghost_width=max(ghost_width(1),ghost_width(2),
     &                      ghost_width(3),ghost_width(4),
     &                      ghost_width(5),ghost_width(6))
        if (max_ghost_width.lt.2*diss_kmax) then
          write(*,*) 'WARNING ... ghost_width < 2*diss_kmax'
          write(*,*) 'max{ghost_width}=',max_ghost_width
          write(*,*) 'diss_kmax=',diss_kmax
          stop
        endif

        ! zero all outer boundary points
        if (phys_bdy(1).eq.1) then 
          do j=1,Ny
           do k=1,Nz
            gb_tt_np1(1,j,k) = 0
            gb_tx_np1(1,j,k) = 0
            gb_ty_np1(1,j,k) = 0
            gb_tz_np1(1,j,k) = 0
            gb_xx_np1(1,j,k) = 0
            gb_xy_np1(1,j,k) = 0
            gb_xz_np1(1,j,k) = 0
            gb_yy_np1(1,j,k) = 0
            gb_yz_np1(1,j,k) = 0
            gb_zz_np1(1,j,k) = 0
            phi1_np1(1,j,k) = 0
           end do
          end do
        end if
        if (phys_bdy(2).eq.1) then 
          do j=1,Ny
           do k=1,Nz
            gb_tt_np1(Nx,j,k) = 0
            gb_tx_np1(Nx,j,k) = 0
            gb_ty_np1(Nx,j,k) = 0
            gb_tz_np1(Nx,j,k) = 0
            gb_xx_np1(Nx,j,k) = 0
            gb_xy_np1(Nx,j,k) = 0
            gb_xz_np1(Nx,j,k) = 0
            gb_yy_np1(Nx,j,k) = 0
            gb_yz_np1(Nx,j,k) = 0
            gb_zz_np1(Nx,j,k) = 0
            phi1_np1(Nx,j,k) = 0
           end do
          end do
        end if
        if (phys_bdy(3).eq.1) then
          do i=1,Nx
           do k=1,Nz
            gb_tt_np1(i,1,k) = 0
            gb_tx_np1(i,1,k) = 0
            gb_ty_np1(i,1,k) = 0
            gb_tz_np1(i,1,k) = 0
            gb_xx_np1(i,1,k) = 0
            gb_xy_np1(i,1,k) = 0
            gb_xz_np1(i,1,k) = 0
            gb_yy_np1(i,1,k) = 0
            gb_yz_np1(i,1,k) = 0
            gb_zz_np1(i,1,k) = 0
            phi1_np1(i,1,k) = 0
           end do
          end do
        end if
        if (phys_bdy(4).eq.1) then 
          do i=1,Nx
           do k=1,Nz
            gb_tt_np1(i,Ny,k) = 0
            gb_tx_np1(i,Ny,k) = 0
            gb_ty_np1(i,Ny,k) = 0
            gb_tz_np1(i,Ny,k) = 0
            gb_xx_np1(i,Ny,k) = 0
            gb_xy_np1(i,Ny,k) = 0
            gb_xz_np1(i,Ny,k) = 0
            gb_yy_np1(i,Ny,k) = 0
            gb_yz_np1(i,Ny,k) = 0
            gb_zz_np1(i,Ny,k) = 0
            phi1_np1(i,Ny,k) = 0
           end do
          end do
        end if
        if (phys_bdy(5).eq.1) then
          do i=1,Nx
           do j=1,Ny
            gb_tt_np1(i,j,1) = 0
            gb_tx_np1(i,j,1) = 0
            gb_ty_np1(i,j,1) = 0
            gb_tz_np1(i,j,1) = 0
            gb_xx_np1(i,j,1) = 0
            gb_xy_np1(i,j,1) = 0
            gb_xz_np1(i,j,1) = 0
            gb_yy_np1(i,j,1) = 0
            gb_yz_np1(i,j,1) = 0
            gb_zz_np1(i,j,1) = 0
            phi1_np1(i,j,1) = 0
           end do
          end do
        end if
        if (phys_bdy(6).eq.1) then
          do i=1,Nx
           do j=1,Ny
            gb_tt_np1(i,j,Nz) = 0
            gb_tx_np1(i,j,Nz) = 0
            gb_ty_np1(i,j,Nz) = 0
            gb_tz_np1(i,j,Nz) = 0
            gb_xx_np1(i,j,Nz) = 0
            gb_xy_np1(i,j,Nz) = 0
            gb_xz_np1(i,j,Nz) = 0
            gb_yy_np1(i,j,Nz) = 0
            gb_yz_np1(i,j,Nz) = 0
            gb_zz_np1(i,j,Nz) = 0
            phi1_np1(i,j,Nz) = 0
           end do
          end do
        end if

        ! define chr2
        do i=is,ie
          do j=js,je
           do k=ks,ke
            if ((chr(i,j,k).ne.ex).and.
     &          (sqrt(x(i)**2+y(j)**2+z(k)**2).ge.(1.0d0-3*dx/2)).and.
     &          ((chr(i-1,j,k).eq.ex).or.(chr(i+1,j,k).eq.ex).or.
     &           (chr(i,j-1,k).eq.ex).or.(chr(i,j+1,k).eq.ex).or.
     &           (chr(i,j,k-1).eq.ex).or.(chr(i,j,k+1).eq.ex))) then
              chr2(i,j,k)=ex
!      write(*,*) 'chr2 excised point'
!      write(*,*) 'i,j,k,x(i),y(j),z(k),sqrt(x(i)**2+y(j)**2+z(k)**2)='
!     &           ,i,j,k,x(i),y(j),z(k),sqrt(x(i)**2+y(j)**2+z(k)**2)
            else 
              chr2(i,j,k)=ex-1
            end if
           end do
          end do
        end do

        !(MAIN LOOP) red-black loop through spacetime points x(i),y(j)  
        do rb=1,2
         do k=ks,ke
          do j=js,je
            do i=is+mod(j+k+rb,2),ie,2
              x0=x(i)
              y0=y(j)
              z0=z(k)
              rho0=sqrt(x0**2+y0**2+z0**2)
              dump=.false.

              if (ltrace) write(*,*) 'i,j,k:',i,j,k

              first_evolved_pt=.false.

              ! define first_evolved_pt
              if (chr(i,j,k).ne.ex.and.chr2(i,j,k).ne.ex.and.
     &             (chr2(i-1,j,k).eq.ex).or.(chr2(i+1,j,k).eq.ex).or.
     &             (chr2(i,j-1,k).eq.ex).or.(chr2(i,j+1,k).eq.ex).or.
     &             (chr2(i,j,k-1).eq.ex).or.(chr2(i,j,k+1).eq.ex)) then
                first_evolved_pt=.true.
              end if

              !(REGION) interior not one-point-away-from-ads-bdy points; evolve 
              if ((chr(i,j,k).ne.ex) 
     &            .and. (chr2(i,j,k).ne.ex)
     &           ) then

!!!!!!!!!!!!!!!!TO TEST efe(a,b)!!!!!!!!!!!!
!               do k=1,Nx
!                do m=1,Ny
!                 gb_tt_np1(k,m)=x(k)+y(m)**2
!                 gb_tt_n(k,m)=x(k)+y(m)**2
!                 gb_tt_nm1(k,m)=x(k)+y(m)**2
!                 gb_tx_np1(k,m)=x(k)+y(m)**2
!                 gb_tx_n(k,m)=x(k)+y(m)**2
!                 gb_tx_nm1(k,m)=x(k)+y(m)**2
!                 gb_ty_np1(k,m)=x(k)+y(m)**2
!                 gb_ty_n(k,m)=x(k)+y(m)**2
!                 gb_ty_nm1(k,m)=x(k)+y(m)**2
!                 gb_xx_np1(k,m)=x(k)+y(m)**2
!                 gb_xx_n(k,m)=x(k)+y(m)**2
!                 gb_xx_nm1(k,m)=x(k)+y(m)**2
!                 gb_xy_np1(k,m)=x(k)+y(m)**2
!                 gb_xy_n(k,m)=x(k)+y(m)**2
!                 gb_xy_nm1(k,m)=x(k)+y(m)**2
!                 gb_yy_np1(k,m)=x(k)+y(m)**2
!                 gb_yy_n(k,m)=x(k)+y(m)**2
!                 gb_yy_nm1(k,m)=x(k)+y(m)**2
!                 gb_zz_np1(k,m)=x(k)+y(m)**2
!                 gb_zz_n(k,m)=x(k)+y(m)**2
!                 gb_zz_nm1(k,m)=x(k)+y(m)**2
!
!                 Hb_t_np1(k,m)=x(k)**3+y(m)**4
!                 Hb_t_n(k,m)=x(k)**3+y(m)**4
!                 Hb_t_nm1(k,m)=x(k)**3+y(m)**4
!                 Hb_x_np1(k,m)=x(k)**3+y(m)**4
!                 Hb_x_n(k,m)=x(k)**3+y(m)**4
!                 Hb_x_nm1(k,m)=x(k)**3+y(m)**4
!                 Hb_y_np1(k,m)=x(k)**3+y(m)**4
!                 Hb_y_n(k,m)=x(k)**3+y(m)**4
!                 Hb_y_nm1(k,m)=x(k)**3+y(m)**4
!               end do
!              end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!      write (*,*) 'DEBUG from g_evo_opt'
!      write (*,*) 'L,x0,y0,z0,rho0,dx=',L,x0,y0,z0,rho0,dx
!      write (*,*) 'gb_tt_np1(i,j,k),gb_tt_n(i,j,k),gb_tt_nm1(i,j,k)='
!     &            ,gb_tt_np1(i,j,k),gb_tt_n(i,j,k),gb_tt_nm1(i,j,k)
!      write (*,*) 'gb_tx_np1(i,j,k),gb_tx_n(i,j,k),gb_tx_nm1(i,j,k)='
!     &            ,gb_tx_np1(i,j,k),gb_tx_n(i,j,k),gb_tx_nm1(i,j,k)
!      write (*,*) 'gb_ty_np1(i,j,k),gb_ty_n(i,j,k),gb_ty_nm1(i,j,k)='
!     &            ,gb_ty_np1(i,j,k),gb_ty_n(i,j,k),gb_ty_nm1(i,j,k)
!      write (*,*) 'gb_tz_np1(i,j,k),gb_tz_n(i,j,k),gb_tz_nm1(i,j,k)='
!     &            ,gb_tz_np1(i,j,k),gb_tz_n(i,j,k),gb_tz_nm1(i,j,k)
!      write (*,*) 'gb_xx_np1(i,j,k),gb_xx_n(i,j,k),gb_xx_nm1(i,j,k)='
!     &            ,gb_xx_np1(i,j,k),gb_xx_n(i,j,k),gb_xx_nm1(i,j,k)
!      write (*,*) 'gb_xy_np1(i,j,k),gb_xy_n(i,j,k),gb_xy_nm1(i,j,k)='
!     &            ,gb_xy_np1(i,j,k),gb_xy_n(i,j,k),gb_xy_nm1(i,j,k)
!      write (*,*) 'gb_xz_np1(i,j,k),gb_xz_n(i,j,k),gb_xz_nm1(i,j,k)='
!     &            ,gb_xz_np1(i,j,k),gb_xz_n(i,j,k),gb_xz_nm1(i,j,k)
!      write (*,*) 'gb_yy_np1(i,j,k),gb_yy_n(i,j,k),gb_yy_nm1(i,j,k)='
!     &            ,gb_yy_np1(i,j,k),gb_yy_n(i,j,k),gb_yy_nm1(i,j,k)
!      write (*,*) 'gb_yz_np1(i,j,k),gb_yz_n(i,j,k),gb_yz_nm1(i,j,k)='
!     &            ,gb_yz_np1(i,j,k),gb_yz_n(i,j,k),gb_yz_nm1(i,j,k)
!      write (*,*) 'gb_zz_np1(i,j,k),gb_zz_n(i,j,k),gb_zz_nm1(i,j,k)='
!     &            ,gb_zz_np1(i,j,k),gb_zz_n(i,j,k),gb_zz_nm1(i,j,k)
!      write (*,*) 'Hb_t_np1(i,j,k),Hb_t_n(i,j,k),Hb_t_nm1(i,j,k)='
!     &            ,Hb_t_np1(i,j,k),Hb_t_n(i,j,k),Hb_t_nm1(i,j,k)
!      write (*,*) 'Hb_x_np1(i,j,k),Hb_x_n(i,j,k),Hb_x_nm1(i,j,k)='
!     &            ,Hb_x_np1(i,j,k),Hb_x_n(i,j,k),Hb_x_nm1(i,j,k)
!      write (*,*) 'Hb_y_np1(i,j,k),Hb_y_n(i,j,k),Hb_y_nm1(i,j,k)='
!     &            ,Hb_y_np1(i,j,k),Hb_y_n(i,j,k),Hb_y_nm1(i,j,k)
!      write (*,*) 'Hb_z_np1(i,j,k),Hb_z_n(i,j,k),Hb_z_nm1(i,j,k)='
!     &            ,Hb_z_np1(i,j,k),Hb_z_n(i,j,k),Hb_z_nm1(i,j,k)

              ! computes tensors at point i,j
              call tensor_init(
     &                gb_tt_np1,gb_tt_n,gb_tt_nm1,
     &                gb_tx_np1,gb_tx_n,gb_tx_nm1,
     &                gb_ty_np1,gb_ty_n,gb_ty_nm1,
     &                gb_tz_np1,gb_tz_n,gb_tz_nm1,
     &                gb_xx_np1,gb_xx_n,gb_xx_nm1,
     &                gb_xy_np1,gb_xy_n,gb_xy_nm1,
     &                gb_xz_np1,gb_xz_n,gb_xz_nm1,
     &                gb_yy_np1,gb_yy_n,gb_yy_nm1,
     &                gb_yz_np1,gb_yz_n,gb_yz_nm1,
     &                gb_zz_np1,gb_zz_n,gb_zz_nm1,
     &                Hb_t_np1,Hb_t_n,Hb_t_nm1,
     &                Hb_x_np1,Hb_x_n,Hb_x_nm1,
     &                Hb_y_np1,Hb_y_n,Hb_y_nm1,
     &                Hb_z_np1,Hb_z_n,Hb_z_nm1,
     &                phi1_np1,phi1_n,phi1_nm1,
     &                g0_ll,g0_uu,g0_ll_x,g0_uu_x,g0_ll_xx,
     &                gads_ll,gads_uu,gads_ll_x,gads_uu_x,gads_ll_xx,
     &                h0_ll,h0_uu,h0_ll_x,h0_uu_x,h0_ll_xx,
     &                A_l,A_l_x,Hads_l,
     &                gamma_ull,gamma_ull_x,
     &                riemann_ulll,ricci_ll,ricci_lu,ricci,
     &                einstein_ll,set_ll,
     &                phi10_x,phi10_xx,
     &                x,y,z,dt,chr,L,ex,Nx,Ny,Nz,i,j,k,
     &                ief_bh_r0,a_rot,kerrads_background,
     &                calc_der,calc_adv_quant)


                do a=1,4
                  do b=1,4
                    do c=1,4
                      dlll(a,b,c)=
     &                    g0_ll_x(b,c,a)-g0_ll_x(a,b,c)+g0_ll_x(c,a,b)
                      gammagg(a,b,c)=0
                      gammahh(a,b,c)=0
                      gammagh(a,b,c)=0
                      gammahg(a,b,c)=0
                      do d=1,4
                        gammagg(a,b,c)=gammagg(a,b,c)
     &                                 +0.5d0*gads_uu(a,d)
     &                                      *(gads_ll_x(c,d,b)
     &                                       -gads_ll_x(b,c,d)
     &                                       +gads_ll_x(d,b,c))
                        gammahh(a,b,c)=gammahh(a,b,c)
     &                                 +0.5d0*h0_uu(a,d)
     &                                      *(h0_ll_x(c,d,b)
     &                                       -h0_ll_x(b,c,d)
     &                                       +h0_ll_x(d,b,c))
                        gammagh(a,b,c)=gammagh(a,b,c)
     &                                 +0.5d0*gads_uu(a,d)
     &                                      *(h0_ll_x(c,d,b)
     &                                       -h0_ll_x(b,c,d)
     &                                       +h0_ll_x(d,b,c))
                        gammahg(a,b,c)=gammahg(a,b,c)
     &                                 +0.5d0*h0_uu(a,d)
     &                                      *(gads_ll_x(c,d,b)
     &                                       -gads_ll_x(b,c,d)
     &                                       +gads_ll_x(d,b,c))
                        cuuuu(a,b,c,d)=gads_uu(a,b)*gads_uu(c,d)+
     &                                 h0_uu(a,b)*h0_uu(c,d)+
     &                                 gads_uu(a,b)*h0_uu(c,d)+
     &                                 h0_uu(a,b)*gads_uu(c,d)
                      end do
                    end do
                  end do
                end do

                !---------------------------------------------------------------- 
                ! Analytically remove the pure AdS terms in the C_l,
                ! denoted ->, so 
                ! C_a = (Hads_a+A_a) - boxx_a
                ! for
                ! 
                ! Hads_a - boxx_a =  (1/2 gads^cd gads_cd,a - gads^cd
                ! gads_ca,d )
                !                   -(1/2 g^cd g_cd,a - g^cd g_ca,d )
                !
                !                 ->-(1/2 h^cd h_cd,a - h^cd h_ca,d )
                !                   -(1/2 gads^cd h_cd,a - gads^cd
                !                   h_ca,d )
                !                   -(1/2 h^cd gads_cd,a - h^cd
                !                   gads_ca,d )
                !----------------------------------------------------------------
                do a=1,4
                  c_l(a)=A_l(a)
     &                   -( 0.5d0*( h0_uu(1,1)*h0_ll_x(1,1,a)+
     &                              h0_uu(2,2)*h0_ll_x(2,2,a)+
     &                              h0_uu(3,3)*h0_ll_x(3,3,a)+
     &                              h0_uu(4,4)*h0_ll_x(4,4,a)+
     &                           2*(h0_uu(1,2)*h0_ll_x(1,2,a)+
     &                              h0_uu(1,3)*h0_ll_x(1,3,a)+
     &                              h0_uu(1,4)*h0_ll_x(1,4,a)+
     &                              h0_uu(2,3)*h0_ll_x(2,3,a)+
     &                              h0_uu(2,4)*h0_ll_x(2,4,a)+
     &                              h0_uu(3,4)*h0_ll_x(3,4,a)) )
     &                     -1.0d0*( h0_uu(1,1)*h0_ll_x(1,a,1)+
     &                              h0_uu(2,2)*h0_ll_x(2,a,2)+
     &                              h0_uu(3,3)*h0_ll_x(3,a,3)+
     &                              h0_uu(4,4)*h0_ll_x(4,a,4)+
     &                             (h0_uu(1,2)*h0_ll_x(1,a,2)+
     &                              h0_uu(2,1)*h0_ll_x(2,a,1)+
     &                              h0_uu(1,3)*h0_ll_x(1,a,3)+
     &                              h0_uu(3,1)*h0_ll_x(3,a,1)+
     &                              h0_uu(1,4)*h0_ll_x(1,a,4)+
     &                              h0_uu(4,1)*h0_ll_x(4,a,1)+
     &                              h0_uu(2,3)*h0_ll_x(2,a,3)+
     &                              h0_uu(3,2)*h0_ll_x(3,a,2)+
     &                              h0_uu(2,4)*h0_ll_x(2,a,4)+
     &                              h0_uu(4,2)*h0_ll_x(4,a,2)+
     &                              h0_uu(3,4)*h0_ll_x(3,a,4)+
     &                              h0_uu(4,3)*h0_ll_x(4,a,3)) ) )
     &                   -( 0.5d0*( gads_uu(1,1)*h0_ll_x(1,1,a)+
     &                              gads_uu(2,2)*h0_ll_x(2,2,a)+
     &                              gads_uu(3,3)*h0_ll_x(3,3,a)+
     &                              gads_uu(4,4)*h0_ll_x(4,4,a)+
     &                           2*(gads_uu(1,2)*h0_ll_x(1,2,a)+
     &                              gads_uu(1,3)*h0_ll_x(1,3,a)+
     &                              gads_uu(1,4)*h0_ll_x(1,4,a)+
     &                              gads_uu(2,3)*h0_ll_x(2,3,a)+
     &                              gads_uu(2,4)*h0_ll_x(2,4,a)+
     &                              gads_uu(3,4)*h0_ll_x(3,4,a)) )
     &                     -1.0d0*( gads_uu(1,1)*h0_ll_x(1,a,1)+
     &                              gads_uu(2,2)*h0_ll_x(2,a,2)+
     &                              gads_uu(3,3)*h0_ll_x(3,a,3)+
     &                              gads_uu(4,4)*h0_ll_x(4,a,4)+
     &                             (gads_uu(1,2)*h0_ll_x(1,a,2)+
     &                              gads_uu(2,1)*h0_ll_x(2,a,1)+
     &                              gads_uu(1,3)*h0_ll_x(1,a,3)+
     &                              gads_uu(3,1)*h0_ll_x(3,a,1)+
     &                              gads_uu(1,4)*h0_ll_x(1,a,4)+
     &                              gads_uu(4,1)*h0_ll_x(4,a,1)+
     &                              gads_uu(2,3)*h0_ll_x(2,a,3)+
     &                              gads_uu(3,2)*h0_ll_x(3,a,2)+
     &                              gads_uu(2,4)*h0_ll_x(2,a,4)+
     &                              gads_uu(4,2)*h0_ll_x(4,a,2)+
     &                              gads_uu(3,4)*h0_ll_x(3,a,4)+
     &                              gads_uu(4,3)*h0_ll_x(4,a,3)) ) )
     &                   -( 0.5d0*( h0_uu(1,1)*gads_ll_x(1,1,a)+
     &                              h0_uu(2,2)*gads_ll_x(2,2,a)+
     &                              h0_uu(3,3)*gads_ll_x(3,3,a)+
     &                              h0_uu(4,4)*gads_ll_x(4,4,a)+
     &                           2*(h0_uu(1,2)*gads_ll_x(1,2,a)+
     &                              h0_uu(1,3)*gads_ll_x(1,3,a)+
     &                              h0_uu(1,4)*gads_ll_x(1,4,a)+
     &                              h0_uu(2,3)*gads_ll_x(2,3,a)+
     &                              h0_uu(2,4)*gads_ll_x(2,4,a)+
     &                              h0_uu(3,4)*gads_ll_x(3,4,a)) )
     &                     -1.0d0*( h0_uu(1,1)*gads_ll_x(1,a,1)+
     &                              h0_uu(2,2)*gads_ll_x(2,a,2)+
     &                              h0_uu(3,3)*gads_ll_x(3,a,3)+
     &                              h0_uu(4,4)*gads_ll_x(4,a,4)+
     &                             (h0_uu(1,2)*gads_ll_x(1,a,2)+
     &                              h0_uu(2,1)*gads_ll_x(2,a,1)+
     &                              h0_uu(1,3)*gads_ll_x(1,a,3)+
     &                              h0_uu(3,1)*gads_ll_x(3,a,1)+
     &                              h0_uu(1,4)*gads_ll_x(1,a,4)+
     &                              h0_uu(4,1)*gads_ll_x(4,a,1)+
     &                              h0_uu(2,3)*gads_ll_x(2,a,3)+
     &                              h0_uu(3,2)*gads_ll_x(3,a,2)+
     &                              h0_uu(2,4)*gads_ll_x(2,a,4)+
     &                              h0_uu(4,2)*gads_ll_x(4,a,2)+
     &                              h0_uu(3,4)*gads_ll_x(3,a,4)+
     &                              h0_uu(4,3)*gads_ll_x(4,a,3)) ) )
                end do

                n_l(1)=-1/sqrt(-g0_uu(1,1))
                do a=1,4
                  n_u(a)=n_l(1)*g0_uu(a,1)+
     &                   n_l(2)*g0_uu(a,2)+
     &                   n_l(3)*g0_uu(a,3)+
     &                   n_l(4)*g0_uu(a,4)
                end do

                ndotc  =n_u(1)*c_l(1)+
     &                  n_u(2)*c_l(2)+
     &                  n_u(3)*c_l(3)+
     &                  n_u(4)*c_l(4)
 
                grad_phi1_sq=phi10_x(1)*phi10_x(1)*g0_uu(1,1)+
     &                       phi10_x(2)*phi10_x(2)*g0_uu(2,2)+
     &                       phi10_x(3)*phi10_x(3)*g0_uu(3,3)+
     &                       phi10_x(4)*phi10_x(4)*g0_uu(4,4)+
     &                    2*(phi10_x(1)*phi10_x(2)*g0_uu(1,2)+
     &                       phi10_x(1)*phi10_x(3)*g0_uu(1,3)+
     &                       phi10_x(1)*phi10_x(4)*g0_uu(1,4)+
     &                       phi10_x(2)*phi10_x(3)*g0_uu(2,3)+
     &                       phi10_x(2)*phi10_x(4)*g0_uu(2,4)+
     &                       phi10_x(3)*phi10_x(4)*g0_uu(3,4))

!      write (*,*) 'DEBUG from g_evo_opt'
!      write (*,*) 'L,x0,y0,z0,rho0,dx=',L,x0,y0,z0,rho0,dx
!      write (*,*) 'phi1_n(i,j,k)=',phi1_n(i,j,k)
!      write (*,*) 'grad_phi1_sq=',grad_phi1_sq
                do a=1,4
                  do b=1,4
                    set_ll(a,b)=phi10_x(a)*phi10_x(b)
     &                         -g0_ll(a,b)*(grad_phi1_sq/2)
                  end do
                end do

                tr_set =set_ll(1,1)*g0_uu(1,1)+
     &                  set_ll(2,2)*g0_uu(2,2)+
     &                  set_ll(3,3)*g0_uu(3,3)+
     &                  set_ll(4,4)*g0_uu(4,4)+
     &               2*(set_ll(1,2)*g0_uu(1,2)+
     &                  set_ll(1,3)*g0_uu(1,3)+
     &                  set_ll(1,4)*g0_uu(1,4)+
     &                  set_ll(2,3)*g0_uu(2,3)+
     &                  set_ll(2,4)*g0_uu(2,4)+
     &                  set_ll(3,4)*g0_uu(3,4))

!      write (*,*) 'tr_set=',tr_set

                !---------------------------------------------------------------- 
                ! Analytically remove the pure AdS terms in the EFEs, denoted ->, so 
                ! efe_ab =   term1_ab + term2_ab + term3_ab + term4_ab 
                !          + term5_ab + term6_ab + term7_ab + term8_ab 
                !          - 8*PI*(se_ab-1/2*tr(se)*g_ab)
                ! for
                ! 
                ! term1_ab = -1/2 g^cd g_ab,cd  
                !          ->-1/2 (h^cd h_ab,cd + gads^cd h_ab,cd + h^cd gads_ab,cd) 
                ! term2_ab = -1/2 g^cd,a g_bc,d
                !          ->-1/2 (h^cd,a h_bc,d + gads^cd,a h_bc,d + h^cd,a gads_bc,d)
                ! term3_ab = -1/2 g^cd,b g_ac,d
                !          ->-1/2 (h^cd,b h_ac,d + gads^cd,b h_ac,d + h^cd,b gads_ac,d)
                ! term4_ab = -1/2 H_a,b
                !          ->-1/2 Hb_a,b
                ! term5_ab = -1/2 H_b,a
                !          ->-1/2 Hb_b,a
                ! term6_ab = H_c G^c_ab
                !          ->Hb_c Ghh^c_ab + Hb_c Ggh^c_ab + Hb_c Ghg^c_ab  
                ! term7_ab = -G^c_db G^d_ca
                !          ->-(  Ghh^c_db Ghh^d_ca + Ggg^c_db Ghh^d_ca
                !              + Ggg^c_db Ggh^d_ca + Ggg^c_db Ghg^d_ca
                !              + Ghh^c_db Ggg^d_ca + Ghh^c_db Ggh^d_ca
                !              + Ghh^c_db Ghg^d_ca + Ggh^c_db Ggg^d_ca
                !              + Ggh^c_db Ghh^d_ca + Ggh^c_db Ggh^d_ca
                !              + Ggh^c_db Ghg^d_ca + Ghg^c_db Ggg^d_ca
                !              + Ghg^c_db Ghh^d_ca + Ghg^c_db Ggh^d_ca
                !              + Ghg^c_db Ghg^d_ca  )
                ! term8_ab = - lambda4 g_ab
                !          ->- lambda4 h_ab
                !
                ! where G   = guu(g_ll_x-g_ll_x+g_ll_x)
                ! where Ggg = gadsuu(gads_ll_x-gads_ll_x+gads_ll_x)
                ! where Ghh = huu(h_ll_x-h_ll_x+h_ll_x)
                ! where Ggh = gadsuu(h_ll_x-h_ll_x+h_ll_x)
                ! where Ghg = huu(gads_ll_x-gads_ll_x+gads_ll_x)
                !
                !---------------------------------------------------------------- 
                do a=1,4
                  do b=a,4
                    term1(a,b)=-0.5d0*(                             
     &                            h0_uu(1,1)*h0_ll_xx(a,b,1,1)+
     &                            h0_uu(2,2)*h0_ll_xx(a,b,2,2)+
     &                            h0_uu(3,3)*h0_ll_xx(a,b,3,3)+
     &                            h0_uu(4,4)*h0_ll_xx(a,b,4,4)+
     &                         2*(h0_uu(1,2)*h0_ll_xx(a,b,1,2)+
     &                            h0_uu(1,3)*h0_ll_xx(a,b,1,3)+
     &                            h0_uu(1,4)*h0_ll_xx(a,b,1,4)+
     &                            h0_uu(2,3)*h0_ll_xx(a,b,2,3)+
     &                            h0_uu(2,4)*h0_ll_xx(a,b,2,4)+
     &                            h0_uu(3,4)*h0_ll_xx(a,b,3,4))
     &                         +
     &                            gads_uu(1,1)*h0_ll_xx(a,b,1,1)+
     &                            gads_uu(2,2)*h0_ll_xx(a,b,2,2)+
     &                            gads_uu(3,3)*h0_ll_xx(a,b,3,3)+
     &                            gads_uu(4,4)*h0_ll_xx(a,b,4,4)+
     &                         2*(gads_uu(1,2)*h0_ll_xx(a,b,1,2)+
     &                            gads_uu(1,3)*h0_ll_xx(a,b,1,3)+
     &                            gads_uu(1,4)*h0_ll_xx(a,b,1,4)+
     &                            gads_uu(2,3)*h0_ll_xx(a,b,2,3)+
     &                            gads_uu(2,4)*h0_ll_xx(a,b,2,4)+
     &                            gads_uu(3,4)*h0_ll_xx(a,b,3,4))
     &                         +
     &                            h0_uu(1,1)*gads_ll_xx(a,b,1,1)+
     &                            h0_uu(2,2)*gads_ll_xx(a,b,2,2)+
     &                            h0_uu(3,3)*gads_ll_xx(a,b,3,3)+
     &                            h0_uu(4,4)*gads_ll_xx(a,b,4,4)+
     &                         2*(h0_uu(1,2)*gads_ll_xx(a,b,1,2)+
     &                            h0_uu(1,3)*gads_ll_xx(a,b,1,3)+
     &                            h0_uu(1,4)*gads_ll_xx(a,b,1,4)+
     &                            h0_uu(2,3)*gads_ll_xx(a,b,2,3)+
     &                            h0_uu(2,4)*gads_ll_xx(a,b,2,4)+
     &                            h0_uu(3,4)*gads_ll_xx(a,b,3,4))
     &                                )
     &
                    term2(a,b)=-0.5d0*(                             
     &                            h0_uu_x(1,1,a)* h0_ll_x(b,1,1) +
     &                            h0_uu_x(1,2,a)*(h0_ll_x(b,1,2) +
     &                                            h0_ll_x(b,2,1))+
     &                            h0_uu_x(1,3,a)*(h0_ll_x(b,1,3) +
     &                                            h0_ll_x(b,3,1))+
     &                            h0_uu_x(1,4,a)*(h0_ll_x(b,1,4) +
     &                                            h0_ll_x(b,4,1))+
     &                            h0_uu_x(2,2,a)* h0_ll_x(b,2,2) +
     &                            h0_uu_x(2,3,a)*(h0_ll_x(b,2,3) +
     &                                            h0_ll_x(b,3,2))+
     &                            h0_uu_x(2,4,a)*(h0_ll_x(b,2,4) +
     &                                            h0_ll_x(b,4,2))+
     &                            h0_uu_x(3,3,a)* h0_ll_x(b,3,3) +
     &                            h0_uu_x(3,4,a)*(h0_ll_x(b,3,4) +
     &                                            h0_ll_x(b,4,3))+
     &                            h0_uu_x(4,4,a)* h0_ll_x(b,4,4) 
     &                         +
     &                            gads_uu_x(1,1,a)* h0_ll_x(b,1,1) +
     &                            gads_uu_x(1,2,a)*(h0_ll_x(b,1,2) +
     &                                              h0_ll_x(b,2,1))+
     &                            gads_uu_x(1,3,a)*(h0_ll_x(b,1,3) +
     &                                              h0_ll_x(b,3,1))+
     &                            gads_uu_x(1,4,a)*(h0_ll_x(b,1,4) +
     &                                              h0_ll_x(b,4,1))+
     &                            gads_uu_x(2,2,a)* h0_ll_x(b,2,2) +
     &                            gads_uu_x(2,3,a)*(h0_ll_x(b,2,3) +
     &                                              h0_ll_x(b,3,2))+
     &                            gads_uu_x(2,4,a)*(h0_ll_x(b,2,4) +
     &                                              h0_ll_x(b,4,2))+
     &                            gads_uu_x(3,3,a)* h0_ll_x(b,3,3) +
     &                            gads_uu_x(3,4,a)*(h0_ll_x(b,3,4) +
     &                                              h0_ll_x(b,4,3))+
     &                            gads_uu_x(4,4,a)* h0_ll_x(b,4,4) 
     &                         +
     &                            h0_uu_x(1,1,a)* gads_ll_x(b,1,1) +
     &                            h0_uu_x(1,2,a)*(gads_ll_x(b,1,2) + 
     &                                            gads_ll_x(b,2,1))+ 
     &                            h0_uu_x(1,3,a)*(gads_ll_x(b,1,3) + 
     &                                            gads_ll_x(b,3,1))+ 
     &                            h0_uu_x(1,4,a)*(gads_ll_x(b,1,4) +
     &                                            gads_ll_x(b,4,1))+
     &                            h0_uu_x(2,2,a)* gads_ll_x(b,2,2) +
     &                            h0_uu_x(2,3,a)*(gads_ll_x(b,2,3) +
     &                                            gads_ll_x(b,3,2))+
     &                            h0_uu_x(2,4,a)*(gads_ll_x(b,2,4) +
     &                                            gads_ll_x(b,4,2))+
     &                            h0_uu_x(3,3,a)* gads_ll_x(b,3,3) +
     &                            h0_uu_x(3,4,a)*(gads_ll_x(b,3,4) +
     &                                            gads_ll_x(b,4,3))+
     &                            h0_uu_x(4,4,a)* gads_ll_x(b,4,4)
     &                              )
     &
                    term3(a,b)=-0.5d0*(                            
     &                            h0_uu_x(1,1,b)* h0_ll_x(a,1,1) +
     &                            h0_uu_x(1,2,b)*(h0_ll_x(a,1,2) +
     &                                            h0_ll_x(a,2,1))+
     &                            h0_uu_x(1,3,b)*(h0_ll_x(a,1,3) +
     &                                            h0_ll_x(a,3,1))+
     &                            h0_uu_x(1,4,b)*(h0_ll_x(a,1,4) +
     &                                            h0_ll_x(a,4,1))+
     &                            h0_uu_x(2,2,b)* h0_ll_x(a,2,2) +
     &                            h0_uu_x(2,3,b)*(h0_ll_x(a,2,3) +
     &                                            h0_ll_x(a,3,2))+
     &                            h0_uu_x(2,4,b)*(h0_ll_x(a,2,4) +
     &                                            h0_ll_x(a,4,2))+
     &                            h0_uu_x(3,3,b)* h0_ll_x(a,3,3) +
     &                            h0_uu_x(3,4,b)*(h0_ll_x(a,3,4) +
     &                                            h0_ll_x(a,4,3))+
     &                            h0_uu_x(4,4,b)* h0_ll_x(a,4,4) 
     &                         +
     &                            gads_uu_x(1,1,b)* h0_ll_x(a,1,1) +
     &                            gads_uu_x(1,2,b)*(h0_ll_x(a,1,2) +
     &                                              h0_ll_x(a,2,1))+
     &                            gads_uu_x(1,3,b)*(h0_ll_x(a,1,3) +
     &                                              h0_ll_x(a,3,1))+
     &                            gads_uu_x(1,4,b)*(h0_ll_x(a,1,4) +
     &                                              h0_ll_x(a,4,1))+
     &                            gads_uu_x(2,2,b)* h0_ll_x(a,2,2) +
     &                            gads_uu_x(2,3,b)*(h0_ll_x(a,2,3) +
     &                                              h0_ll_x(a,3,2))+
     &                            gads_uu_x(2,4,b)*(h0_ll_x(a,2,4) +
     &                                              h0_ll_x(a,4,2))+
     &                            gads_uu_x(3,3,b)* h0_ll_x(a,3,3) +
     &                            gads_uu_x(3,4,b)*(h0_ll_x(a,3,4) +
     &                                              h0_ll_x(a,4,3))+
     &                            gads_uu_x(4,4,b)* h0_ll_x(a,4,4) 
     &                         +
     &                            h0_uu_x(1,1,b)* gads_ll_x(a,1,1) +
     &                            h0_uu_x(1,2,b)*(gads_ll_x(a,1,2) +  
     &                                            gads_ll_x(a,2,1))+  
     &                            h0_uu_x(1,3,b)*(gads_ll_x(a,1,3) +  
     &                                            gads_ll_x(a,3,1))+   
     &                            h0_uu_x(1,4,b)*(gads_ll_x(a,1,4) +
     &                                            gads_ll_x(a,4,1))+
     &                            h0_uu_x(2,2,b)* gads_ll_x(a,2,2) +
     &                            h0_uu_x(2,3,b)*(gads_ll_x(a,2,3) +
     &                                            gads_ll_x(a,3,2))+
     &                            h0_uu_x(2,4,b)*(gads_ll_x(a,2,4) +
     &                                            gads_ll_x(a,4,2))+
     &                            h0_uu_x(3,3,b)* gads_ll_x(a,3,3) +
     &                            h0_uu_x(3,4,b)*(gads_ll_x(a,3,4) +
     &                                            gads_ll_x(a,4,3))+
     &                            h0_uu_x(4,4,b)* gads_ll_x(a,4,4) 
     &                              )

                    term4(a,b)=-0.5d0*A_l_x(a,b)                  

                    term5(a,b)=-0.5d0*A_l_x(b,a)           

                    term6(a,b)=     (           
     &                            Hads_l(1)*gammahh(1,a,b)+      
     &                            Hads_l(2)*gammahh(2,a,b)+
     &                            Hads_l(3)*gammahh(3,a,b)+
     &                            Hads_l(4)*gammahh(4,a,b)
     &                         +
     &                            Hads_l(1)*gammagh(1,a,b)+
     &                            Hads_l(2)*gammagh(2,a,b)+
     &                            Hads_l(3)*gammagh(3,a,b)+
     &                            Hads_l(4)*gammagh(4,a,b)
     &                         +
     &                            Hads_l(1)*gammahg(1,a,b)+
     &                            Hads_l(2)*gammahg(2,a,b)+
     &                            Hads_l(3)*gammahg(3,a,b)+
     &                            Hads_l(4)*gammahg(4,a,b)
     &                         +                  
     &                            A_l(1)*gammagg(1,a,b)  +
     &                            A_l(2)*gammagg(2,a,b)  +
     &                            A_l(3)*gammagg(3,a,b)  +
     &                            A_l(4)*gammagg(4,a,b)  
     &                         +
     &                            A_l(1)*gammahh(1,a,b)  +
     &                            A_l(2)*gammahh(2,a,b)  +
     &                            A_l(3)*gammahh(3,a,b)  +
     &                            A_l(4)*gammahh(4,a,b)  
     &                         +
     &                            A_l(1)*gammagh(1,a,b)  +
     &                            A_l(2)*gammagh(2,a,b)  +
     &                            A_l(3)*gammagh(3,a,b)  +
     &                            A_l(4)*gammagh(4,a,b)  
     &                         +
     &                            A_l(1)*gammahg(1,a,b)  +
     &                            A_l(2)*gammahg(2,a,b)  +
     &                            A_l(3)*gammahg(3,a,b)  +
     &                            A_l(4)*gammahg(4,a,b)  
     &                              ) 
                             
                    term7(a,b)=    -(
     &                            gammahh(1,1,b)*gammahh(1,1,a)+
     &                            gammahh(1,2,b)*gammahh(2,1,a)+
     &                            gammahh(1,3,b)*gammahh(3,1,a)+
     &                            gammahh(1,4,b)*gammahh(4,1,a)+
     &                            gammahh(2,1,b)*gammahh(1,2,a)+
     &                            gammahh(2,2,b)*gammahh(2,2,a)+
     &                            gammahh(2,3,b)*gammahh(3,2,a)+
     &                            gammahh(2,4,b)*gammahh(4,2,a)+
     &                            gammahh(3,1,b)*gammahh(1,3,a)+
     &                            gammahh(3,2,b)*gammahh(2,3,a)+
     &                            gammahh(3,3,b)*gammahh(3,3,a)+
     &                            gammahh(3,4,b)*gammahh(4,3,a)+
     &                            gammahh(4,1,b)*gammahh(1,4,a)+
     &                            gammahh(4,2,b)*gammahh(2,4,a)+
     &                            gammahh(4,3,b)*gammahh(3,4,a)+
     &                            gammahh(4,4,b)*gammahh(4,4,a)
     &                         +
     &                            gammagg(1,1,b)*gammahh(1,1,a)+
     &                            gammagg(1,2,b)*gammahh(2,1,a)+
     &                            gammagg(1,3,b)*gammahh(3,1,a)+
     &                            gammagg(1,4,b)*gammahh(4,1,a)+
     &                            gammagg(2,1,b)*gammahh(1,2,a)+
     &                            gammagg(2,2,b)*gammahh(2,2,a)+
     &                            gammagg(2,3,b)*gammahh(3,2,a)+
     &                            gammagg(2,4,b)*gammahh(4,2,a)+
     &                            gammagg(3,1,b)*gammahh(1,3,a)+
     &                            gammagg(3,2,b)*gammahh(2,3,a)+
     &                            gammagg(3,3,b)*gammahh(3,3,a)+
     &                            gammagg(3,4,b)*gammahh(4,3,a)+
     &                            gammagg(4,1,b)*gammahh(1,4,a)+
     &                            gammagg(4,2,b)*gammahh(2,4,a)+
     &                            gammagg(4,3,b)*gammahh(3,4,a)+
     &                            gammagg(4,4,b)*gammahh(4,4,a)
     &                         +
     &                            gammagg(1,1,b)*gammagh(1,1,a)+
     &                            gammagg(1,2,b)*gammagh(2,1,a)+
     &                            gammagg(1,3,b)*gammagh(3,1,a)+
     &                            gammagg(1,4,b)*gammagh(4,1,a)+
     &                            gammagg(2,1,b)*gammagh(1,2,a)+
     &                            gammagg(2,2,b)*gammagh(2,2,a)+
     &                            gammagg(2,3,b)*gammagh(3,2,a)+
     &                            gammagg(2,4,b)*gammagh(4,2,a)+
     &                            gammagg(3,1,b)*gammagh(1,3,a)+
     &                            gammagg(3,2,b)*gammagh(2,3,a)+
     &                            gammagg(3,3,b)*gammagh(3,3,a)+
     &                            gammagg(3,4,b)*gammagh(4,3,a)+
     &                            gammagg(4,1,b)*gammagh(1,4,a)+
     &                            gammagg(4,2,b)*gammagh(2,4,a)+
     &                            gammagg(4,3,b)*gammagh(3,4,a)+
     &                            gammagg(4,4,b)*gammagh(4,4,a)
     &                         +
     &                            gammagg(1,1,b)*gammahg(1,1,a)+
     &                            gammagg(1,2,b)*gammahg(2,1,a)+
     &                            gammagg(1,3,b)*gammahg(3,1,a)+
     &                            gammagg(1,4,b)*gammahg(4,1,a)+
     &                            gammagg(2,1,b)*gammahg(1,2,a)+
     &                            gammagg(2,2,b)*gammahg(2,2,a)+
     &                            gammagg(2,3,b)*gammahg(3,2,a)+
     &                            gammagg(2,4,b)*gammahg(4,2,a)+
     &                            gammagg(3,1,b)*gammahg(1,3,a)+
     &                            gammagg(3,2,b)*gammahg(2,3,a)+
     &                            gammagg(3,3,b)*gammahg(3,3,a)+
     &                            gammagg(3,4,b)*gammahg(4,3,a)+
     &                            gammagg(4,1,b)*gammahg(1,4,a)+
     &                            gammagg(4,2,b)*gammahg(2,4,a)+
     &                            gammagg(4,3,b)*gammahg(3,4,a)+
     &                            gammagg(4,4,b)*gammahg(4,4,a)
     &                         +
     &                            gammahh(1,1,b)*gammagg(1,1,a)+
     &                            gammahh(1,2,b)*gammagg(2,1,a)+
     &                            gammahh(1,3,b)*gammagg(3,1,a)+
     &                            gammahh(1,4,b)*gammagg(4,1,a)+
     &                            gammahh(2,1,b)*gammagg(1,2,a)+
     &                            gammahh(2,2,b)*gammagg(2,2,a)+
     &                            gammahh(2,3,b)*gammagg(3,2,a)+
     &                            gammahh(2,4,b)*gammagg(4,2,a)+
     &                            gammahh(3,1,b)*gammagg(1,3,a)+
     &                            gammahh(3,2,b)*gammagg(2,3,a)+
     &                            gammahh(3,3,b)*gammagg(3,3,a)+
     &                            gammahh(3,4,b)*gammagg(4,3,a)+
     &                            gammahh(4,1,b)*gammagg(1,4,a)+
     &                            gammahh(4,2,b)*gammagg(2,4,a)+
     &                            gammahh(4,3,b)*gammagg(3,4,a)+
     &                            gammahh(4,4,b)*gammagg(4,4,a)
     &                         +
     &                            gammahh(1,1,b)*gammagh(1,1,a)+
     &                            gammahh(1,2,b)*gammagh(2,1,a)+
     &                            gammahh(1,3,b)*gammagh(3,1,a)+
     &                            gammahh(1,4,b)*gammagh(4,1,a)+
     &                            gammahh(2,1,b)*gammagh(1,2,a)+
     &                            gammahh(2,2,b)*gammagh(2,2,a)+
     &                            gammahh(2,3,b)*gammagh(3,2,a)+
     &                            gammahh(2,4,b)*gammagh(4,2,a)+
     &                            gammahh(3,1,b)*gammagh(1,3,a)+
     &                            gammahh(3,2,b)*gammagh(2,3,a)+
     &                            gammahh(3,3,b)*gammagh(3,3,a)+
     &                            gammahh(3,4,b)*gammagh(4,3,a)+
     &                            gammahh(4,1,b)*gammagh(1,4,a)+
     &                            gammahh(4,2,b)*gammagh(2,4,a)+
     &                            gammahh(4,3,b)*gammagh(3,4,a)+
     &                            gammahh(4,4,b)*gammagh(4,4,a)
     &                         +
     &                            gammahh(1,1,b)*gammahg(1,1,a)+
     &                            gammahh(1,2,b)*gammahg(2,1,a)+
     &                            gammahh(1,3,b)*gammahg(3,1,a)+
     &                            gammahh(1,4,b)*gammahg(4,1,a)+
     &                            gammahh(2,1,b)*gammahg(1,2,a)+
     &                            gammahh(2,2,b)*gammahg(2,2,a)+
     &                            gammahh(2,3,b)*gammahg(3,2,a)+
     &                            gammahh(2,4,b)*gammahg(4,2,a)+
     &                            gammahh(3,1,b)*gammahg(1,3,a)+
     &                            gammahh(3,2,b)*gammahg(2,3,a)+
     &                            gammahh(3,3,b)*gammahg(3,3,a)+
     &                            gammahh(3,4,b)*gammahg(4,3,a)+
     &                            gammahh(4,1,b)*gammahg(1,4,a)+
     &                            gammahh(4,2,b)*gammahg(2,4,a)+
     &                            gammahh(4,3,b)*gammahg(3,4,a)+
     &                            gammahh(4,4,b)*gammahg(4,4,a)
     &                         +
     &                            gammagh(1,1,b)*gammagg(1,1,a)+
     &                            gammagh(1,2,b)*gammagg(2,1,a)+
     &                            gammagh(1,3,b)*gammagg(3,1,a)+
     &                            gammagh(1,4,b)*gammagg(4,1,a)+
     &                            gammagh(2,1,b)*gammagg(1,2,a)+
     &                            gammagh(2,2,b)*gammagg(2,2,a)+
     &                            gammagh(2,3,b)*gammagg(3,2,a)+
     &                            gammagh(2,4,b)*gammagg(4,2,a)+
     &                            gammagh(3,1,b)*gammagg(1,3,a)+
     &                            gammagh(3,2,b)*gammagg(2,3,a)+
     &                            gammagh(3,3,b)*gammagg(3,3,a)+
     &                            gammagh(3,4,b)*gammagg(4,3,a)+
     &                            gammagh(4,1,b)*gammagg(1,4,a)+
     &                            gammagh(4,2,b)*gammagg(2,4,a)+
     &                            gammagh(4,3,b)*gammagg(3,4,a)+
     &                            gammagh(4,4,b)*gammagg(4,4,a)
     &                         +
     &                            gammagh(1,1,b)*gammahh(1,1,a)+
     &                            gammagh(1,2,b)*gammahh(2,1,a)+
     &                            gammagh(1,3,b)*gammahh(3,1,a)+
     &                            gammagh(1,4,b)*gammahh(4,1,a)+
     &                            gammagh(2,1,b)*gammahh(1,2,a)+
     &                            gammagh(2,2,b)*gammahh(2,2,a)+
     &                            gammagh(2,3,b)*gammahh(3,2,a)+
     &                            gammagh(2,4,b)*gammahh(4,2,a)+
     &                            gammagh(3,1,b)*gammahh(1,3,a)+
     &                            gammagh(3,2,b)*gammahh(2,3,a)+
     &                            gammagh(3,3,b)*gammahh(3,3,a)+
     &                            gammagh(3,4,b)*gammahh(4,3,a)+
     &                            gammagh(4,1,b)*gammahh(1,4,a)+
     &                            gammagh(4,2,b)*gammahh(2,4,a)+
     &                            gammagh(4,3,b)*gammahh(3,4,a)+
     &                            gammagh(4,4,b)*gammahh(4,4,a)
     &                         +
     &                            gammagh(1,1,b)*gammagh(1,1,a)+
     &                            gammagh(1,2,b)*gammagh(2,1,a)+
     &                            gammagh(1,3,b)*gammagh(3,1,a)+
     &                            gammagh(1,4,b)*gammagh(4,1,a)+
     &                            gammagh(2,1,b)*gammagh(1,2,a)+
     &                            gammagh(2,2,b)*gammagh(2,2,a)+
     &                            gammagh(2,3,b)*gammagh(3,2,a)+
     &                            gammagh(2,4,b)*gammagh(4,2,a)+
     &                            gammagh(3,1,b)*gammagh(1,3,a)+
     &                            gammagh(3,2,b)*gammagh(2,3,a)+
     &                            gammagh(3,3,b)*gammagh(3,3,a)+
     &                            gammagh(3,4,b)*gammagh(4,3,a)+
     &                            gammagh(4,1,b)*gammagh(1,4,a)+
     &                            gammagh(4,2,b)*gammagh(2,4,a)+
     &                            gammagh(4,3,b)*gammagh(3,4,a)+
     &                            gammagh(4,4,b)*gammagh(4,4,a)
     &                         +
     &                            gammagh(1,1,b)*gammahg(1,1,a)+
     &                            gammagh(1,2,b)*gammahg(2,1,a)+
     &                            gammagh(1,3,b)*gammahg(3,1,a)+
     &                            gammagh(1,4,b)*gammahg(4,1,a)+
     &                            gammagh(2,1,b)*gammahg(1,2,a)+
     &                            gammagh(2,2,b)*gammahg(2,2,a)+
     &                            gammagh(2,3,b)*gammahg(3,2,a)+
     &                            gammagh(2,4,b)*gammahg(4,2,a)+
     &                            gammagh(3,1,b)*gammahg(1,3,a)+
     &                            gammagh(3,2,b)*gammahg(2,3,a)+
     &                            gammagh(3,3,b)*gammahg(3,3,a)+
     &                            gammagh(3,4,b)*gammahg(4,3,a)+
     &                            gammagh(4,1,b)*gammahg(1,4,a)+
     &                            gammagh(4,2,b)*gammahg(2,4,a)+
     &                            gammagh(4,3,b)*gammahg(3,4,a)+
     &                            gammagh(4,4,b)*gammahg(4,4,a)
     &                         +
     &                            gammahg(1,1,b)*gammagg(1,1,a)+
     &                            gammahg(1,2,b)*gammagg(2,1,a)+
     &                            gammahg(1,3,b)*gammagg(3,1,a)+
     &                            gammahg(1,4,b)*gammagg(4,1,a)+
     &                            gammahg(2,1,b)*gammagg(1,2,a)+
     &                            gammahg(2,2,b)*gammagg(2,2,a)+
     &                            gammahg(2,3,b)*gammagg(3,2,a)+
     &                            gammahg(2,4,b)*gammagg(4,2,a)+
     &                            gammahg(3,1,b)*gammagg(1,3,a)+
     &                            gammahg(3,2,b)*gammagg(2,3,a)+
     &                            gammahg(3,3,b)*gammagg(3,3,a)+
     &                            gammahg(3,4,b)*gammagg(4,3,a)+
     &                            gammahg(4,1,b)*gammagg(1,4,a)+
     &                            gammahg(4,2,b)*gammagg(2,4,a)+
     &                            gammahg(4,3,b)*gammagg(3,4,a)+
     &                            gammahg(4,4,b)*gammagg(4,4,a)
     &                         +
     &                            gammahg(1,1,b)*gammahh(1,1,a)+
     &                            gammahg(1,2,b)*gammahh(2,1,a)+
     &                            gammahg(1,3,b)*gammahh(3,1,a)+
     &                            gammahg(1,4,b)*gammahh(4,1,a)+
     &                            gammahg(2,1,b)*gammahh(1,2,a)+
     &                            gammahg(2,2,b)*gammahh(2,2,a)+
     &                            gammahg(2,3,b)*gammahh(3,2,a)+
     &                            gammahg(2,4,b)*gammahh(4,2,a)+
     &                            gammahg(3,1,b)*gammahh(1,3,a)+
     &                            gammahg(3,2,b)*gammahh(2,3,a)+
     &                            gammahg(3,3,b)*gammahh(3,3,a)+
     &                            gammahg(3,4,b)*gammahh(4,3,a)+
     &                            gammahg(4,1,b)*gammahh(1,4,a)+
     &                            gammahg(4,2,b)*gammahh(2,4,a)+
     &                            gammahg(4,3,b)*gammahh(3,4,a)+
     &                            gammahg(4,4,b)*gammahh(4,4,a)
     &                         +
     &                            gammahg(1,1,b)*gammagh(1,1,a)+
     &                            gammahg(1,2,b)*gammagh(2,1,a)+
     &                            gammahg(1,3,b)*gammagh(3,1,a)+
     &                            gammahg(1,4,b)*gammagh(4,1,a)+
     &                            gammahg(2,1,b)*gammagh(1,2,a)+
     &                            gammahg(2,2,b)*gammagh(2,2,a)+
     &                            gammahg(2,3,b)*gammagh(3,2,a)+
     &                            gammahg(2,4,b)*gammagh(4,2,a)+
     &                            gammahg(3,1,b)*gammagh(1,3,a)+
     &                            gammahg(3,2,b)*gammagh(2,3,a)+
     &                            gammahg(3,3,b)*gammagh(3,3,a)+
     &                            gammahg(3,4,b)*gammagh(4,3,a)+
     &                            gammahg(4,1,b)*gammagh(1,4,a)+
     &                            gammahg(4,2,b)*gammagh(2,4,a)+
     &                            gammahg(4,3,b)*gammagh(3,4,a)+
     &                            gammahg(4,4,b)*gammagh(4,4,a)
     &                         +
     &                            gammahg(1,1,b)*gammahg(1,1,a)+
     &                            gammahg(1,2,b)*gammahg(2,1,a)+
     &                            gammahg(1,3,b)*gammahg(3,1,a)+
     &                            gammahg(1,4,b)*gammahg(4,1,a)+
     &                            gammahg(2,1,b)*gammahg(1,2,a)+
     &                            gammahg(2,2,b)*gammahg(2,2,a)+
     &                            gammahg(2,3,b)*gammahg(3,2,a)+
     &                            gammahg(2,4,b)*gammahg(4,2,a)+
     &                            gammahg(3,1,b)*gammahg(1,3,a)+
     &                            gammahg(3,2,b)*gammahg(2,3,a)+
     &                            gammahg(3,3,b)*gammahg(3,3,a)+
     &                            gammahg(3,4,b)*gammahg(4,3,a)+
     &                            gammahg(4,1,b)*gammahg(1,4,a)+
     &                            gammahg(4,2,b)*gammahg(2,4,a)+
     &                            gammahg(4,3,b)*gammahg(3,4,a)+
     &                            gammahg(4,4,b)*gammahg(4,4,a)
     &                              )
                    term8(a,b)=-lambda4*h0_ll(a,b)

                    efe(a,b)=term1(a,b)+term2(a,b)+term3(a,b)+term4(a,b)
     &                      +term5(a,b)+term6(a,b)+term7(a,b)+term8(a,b)
     &                      -8*PI*(set_ll(a,b)-tr_set*g0_ll(a,b)/2)

                  end do
                end do


!!!!!!!!!!!!!!!!TO TEST efe(a,b)!!!!!!!!!!!!
!                do a=1,4
!                 do b=a,4
!                  write(*,*)'DEBUG from g_evo_opt'
!                  write(*,*)'i,j,k,x(i),y(j),z(k)=',i,j,k,x(i),y(j),z(k)
!                  write(*,*)'a,b,efe(a,b)=',a,b,efe(a,b)
!                  write(*,*)'term1(a,b)=',term1(a,b)
!                  write(*,*)'term2(a,b)=',term2(a,b)
!                  write(*,*)'term3(a,b)=',term3(a,b)
!                  write(*,*)'term4(a,b)=',term4(a,b)
!                  write(*,*)'term5(a,b)=',term5(a,b)
!                  write(*,*)'term6(a,b)=',term6(a,b)
!                  write(*,*)'term7(a,b)=',term7(a,b)
!                  write(*,*)'term8(a,b)=',term8(a,b)
!                  write(*,*)'a,b,set_ll(a,b),tr_set='
!     &                      ,a,b,set_ll(a,b),tr_set
!
!                 end do
!                end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                !--------------------------------------------------------------------------
                ! phi1_res = phi1,ab g^ab + phi1,b g^ab,a + phi1,c g^cb gamma^a_ab
                !         (= g^ab phi1,ab - g^ab gamma^c_ab phi1,c) 
                !--------------------------------------------------------------------------
                phi1_res= phi10_xx(1,1)*g0_uu(1,1)+
     &                    phi10_xx(2,2)*g0_uu(2,2)+
     &                    phi10_xx(3,3)*g0_uu(3,3)+
     &                    phi10_xx(4,4)*g0_uu(4,4)+
     &                 2*(phi10_xx(1,2)*g0_uu(1,2)+
     &                    phi10_xx(1,3)*g0_uu(1,3)+
     &                    phi10_xx(1,4)*g0_uu(1,4)+
     &                    phi10_xx(2,3)*g0_uu(2,3)+
     &                    phi10_xx(2,4)*g0_uu(2,4)+
     &                    phi10_xx(3,4)*g0_uu(3,4))
     &                +
     &                    phi10_x(1)*g0_uu_x(1,1,1)+
     &                    phi10_x(1)*g0_uu_x(2,1,2)+
     &                    phi10_x(1)*g0_uu_x(3,1,3)+
     &                    phi10_x(1)*g0_uu_x(4,1,4)+
     &                    phi10_x(2)*g0_uu_x(1,2,1)+
     &                    phi10_x(2)*g0_uu_x(2,2,2)+
     &                    phi10_x(2)*g0_uu_x(3,2,3)+
     &                    phi10_x(2)*g0_uu_x(4,2,4)+
     &                    phi10_x(3)*g0_uu_x(1,3,1)+
     &                    phi10_x(3)*g0_uu_x(2,3,2)+
     &                    phi10_x(3)*g0_uu_x(3,3,3)+
     &                    phi10_x(3)*g0_uu_x(4,3,4)+
     &                    phi10_x(4)*g0_uu_x(1,4,1)+
     &                    phi10_x(4)*g0_uu_x(2,4,2)+
     &                    phi10_x(4)*g0_uu_x(3,4,3)+
     &                    phi10_x(4)*g0_uu_x(4,4,4)
     &                +
     &                    phi10_x(1)*g0_uu(1,1)*gamma_ull(1,1,1)+
     &                    phi10_x(1)*g0_uu(1,1)*gamma_ull(2,2,1)+ 
     &                    phi10_x(1)*g0_uu(1,1)*gamma_ull(3,3,1)+
     &                    phi10_x(1)*g0_uu(1,1)*gamma_ull(4,4,1)+
     &                    phi10_x(1)*g0_uu(1,2)*gamma_ull(1,1,2)+
     &                    phi10_x(1)*g0_uu(1,2)*gamma_ull(2,2,2)+
     &                    phi10_x(1)*g0_uu(1,2)*gamma_ull(3,3,2)+
     &                    phi10_x(1)*g0_uu(1,2)*gamma_ull(4,4,2)+
     &                    phi10_x(1)*g0_uu(1,3)*gamma_ull(1,1,3)+
     &                    phi10_x(1)*g0_uu(1,3)*gamma_ull(2,2,3)+
     &                    phi10_x(1)*g0_uu(1,3)*gamma_ull(3,3,3)+
     &                    phi10_x(1)*g0_uu(1,3)*gamma_ull(4,4,3)+
     &                    phi10_x(1)*g0_uu(1,4)*gamma_ull(1,1,4)+
     &                    phi10_x(1)*g0_uu(1,4)*gamma_ull(2,2,4)+
     &                    phi10_x(1)*g0_uu(1,4)*gamma_ull(3,3,4)+
     &                    phi10_x(1)*g0_uu(1,4)*gamma_ull(4,4,4)+
     &                    phi10_x(2)*g0_uu(2,1)*gamma_ull(1,1,1)+
     &                    phi10_x(2)*g0_uu(2,1)*gamma_ull(2,2,1)+ 
     &                    phi10_x(2)*g0_uu(2,1)*gamma_ull(3,3,1)+
     &                    phi10_x(2)*g0_uu(2,1)*gamma_ull(4,4,1)+
     &                    phi10_x(2)*g0_uu(2,2)*gamma_ull(1,1,2)+
     &                    phi10_x(2)*g0_uu(2,2)*gamma_ull(2,2,2)+
     &                    phi10_x(2)*g0_uu(2,2)*gamma_ull(3,3,2)+
     &                    phi10_x(2)*g0_uu(2,2)*gamma_ull(4,4,2)+
     &                    phi10_x(2)*g0_uu(2,3)*gamma_ull(1,1,3)+
     &                    phi10_x(2)*g0_uu(2,3)*gamma_ull(2,2,3)+
     &                    phi10_x(2)*g0_uu(2,3)*gamma_ull(3,3,3)+
     &                    phi10_x(2)*g0_uu(2,3)*gamma_ull(4,4,3)+
     &                    phi10_x(2)*g0_uu(2,4)*gamma_ull(1,1,4)+
     &                    phi10_x(2)*g0_uu(2,4)*gamma_ull(2,2,4)+
     &                    phi10_x(2)*g0_uu(2,4)*gamma_ull(3,3,4)+
     &                    phi10_x(2)*g0_uu(2,4)*gamma_ull(4,4,4)+
     &                    phi10_x(3)*g0_uu(3,1)*gamma_ull(1,1,1)+
     &                    phi10_x(3)*g0_uu(3,1)*gamma_ull(2,2,1)+ 
     &                    phi10_x(3)*g0_uu(3,1)*gamma_ull(3,3,1)+
     &                    phi10_x(3)*g0_uu(3,1)*gamma_ull(4,4,1)+
     &                    phi10_x(3)*g0_uu(3,2)*gamma_ull(1,1,2)+
     &                    phi10_x(3)*g0_uu(3,2)*gamma_ull(2,2,2)+
     &                    phi10_x(3)*g0_uu(3,2)*gamma_ull(3,3,2)+
     &                    phi10_x(3)*g0_uu(3,2)*gamma_ull(4,4,2)+
     &                    phi10_x(3)*g0_uu(3,3)*gamma_ull(1,1,3)+
     &                    phi10_x(3)*g0_uu(3,3)*gamma_ull(2,2,3)+
     &                    phi10_x(3)*g0_uu(3,3)*gamma_ull(3,3,3)+
     &                    phi10_x(3)*g0_uu(3,3)*gamma_ull(4,4,3)+
     &                    phi10_x(3)*g0_uu(3,4)*gamma_ull(1,1,4)+
     &                    phi10_x(3)*g0_uu(3,4)*gamma_ull(2,2,4)+
     &                    phi10_x(3)*g0_uu(3,4)*gamma_ull(3,3,4)+
     &                    phi10_x(3)*g0_uu(3,4)*gamma_ull(4,4,4)+
     &                    phi10_x(4)*g0_uu(4,1)*gamma_ull(1,1,1)+
     &                    phi10_x(4)*g0_uu(4,1)*gamma_ull(2,2,1)+ 
     &                    phi10_x(4)*g0_uu(4,1)*gamma_ull(3,3,1)+
     &                    phi10_x(4)*g0_uu(4,1)*gamma_ull(4,4,1)+
     &                    phi10_x(4)*g0_uu(4,2)*gamma_ull(1,1,2)+
     &                    phi10_x(4)*g0_uu(4,2)*gamma_ull(2,2,2)+
     &                    phi10_x(4)*g0_uu(4,2)*gamma_ull(3,3,2)+
     &                    phi10_x(4)*g0_uu(4,2)*gamma_ull(4,4,2)+
     &                    phi10_x(4)*g0_uu(4,3)*gamma_ull(1,1,3)+
     &                    phi10_x(4)*g0_uu(4,3)*gamma_ull(2,2,3)+
     &                    phi10_x(4)*g0_uu(4,3)*gamma_ull(3,3,3)+
     &                    phi10_x(4)*g0_uu(4,3)*gamma_ull(4,4,3)+
     &                    phi10_x(4)*g0_uu(4,4)*gamma_ull(1,1,4)+
     &                    phi10_x(4)*g0_uu(4,4)*gamma_ull(2,2,4)+
     &                    phi10_x(4)*g0_uu(4,4)*gamma_ull(3,3,4)+
     &                    phi10_x(4)*g0_uu(4,4)*gamma_ull(4,4,4)

!       write (*,*) 'phi1_res=',phi1_res

                !---------------------------------------------------------------- 
                ! computes diag. Jacobian of g_np1->L.g_np1 transformation
                ! by differentiating L.g wrt. g(a,b)_ij_np1 diag. entries
                ! 
                ! ddgb_J_tx,ddgb_J_ty differ from ddgb_J due to forward/backward stencils
                ! at excision surfaces that affect the cross-derivatives tx,ty 
                ! (these are the only contributions, since the diag Jacobian is diff wrt. g_ij_np1)
                !---------------------------------------------------------------- 
                dgb_J=1.0d0/2.0d0/dt
                ddgb_J=1.0d0/dt/dt

!!!!!!!!MY VERSION!!!!!!!

!i
        if (i.eq.1) then
               if ((.not.extrap)
     &            .and.(chr(i+1,j,k).ne.ex)
     &            .and.(chr(i+2,j,k).ne.ex)
     &            .and.(chr(i+3,j,k).ne.ex)) then
                   ddgb_J_tx=-1/dt/dx
               else if ((chr(i+1,j,k).ne.ex
     &                 .and.chr(i+2,j,k).ne.ex)) then
                   ddgb_J_tx=-3/4/dt/dx
               else if (chr(i+1,j,k).ne.ex) then
                   ddgb_J_tx=-1/2/dt/dx
!                  write(*,*) 'g_evo_opt: warning ... i=1 first order'
!                  write(*,*) '    i,j,Nx,Ny,dx=',i,j,Nx,Ny,dx
               else
                   ddgb_J_tx=0
                   return
               end if
        else if (i.eq.2) then
         if ((chr(i-1,j,k).ne.ex).and.(chr(i+1,j,k).ne.ex)) then
                   ddgb_J_tx=0
         else if (chr(i-1,j,k).eq.ex) then
               if ((.not.extrap)
     &            .and.(chr(i+1,j,k).ne.ex)
     &            .and.(chr(i+2,j,k).ne.ex)
     &            .and.(chr(i+3,j,k).ne.ex)) then
                   ddgb_J_tx=-1/dt/dx
               else if ((chr(i+1,j,k).ne.ex)
     &                 .and.(chr(i+2,j,k).ne.ex)) then
                   ddgb_J_tx=-3/4/dt/dx
               else if (chr(i+1,j,k).ne.ex) then
                   ddgb_J_tx=-1/2/dt/dx
!                  write(*,*) 'g_evo_opt: warning ... i=1 first order'
!                  write(*,*) '    i,j,Nx,Ny,dx=',i,j,Nx,Ny,dx
               else
                   ddgb_J_tx=0
                   return
               end if
         else   !this is the case where (i-1,j,k) is not excised and (i+1,j,k) is excised 
                   ddgb_J_tx=1/2/dt/dx
!                  write(*,*) 'g_evo_opt: warning ... i=1 first order'
!                  write(*,*) '    i,j,Nx,Ny,dx=',i,j,Nx,Ny,dx
         end if
        else if (i.eq.3) then
         if ((chr(i-1,j,k).ne.ex).and.(chr(i+1,j,k).ne.ex)) then
                   ddgb_J_tx=0
         else if (chr(i-1,j,k).eq.ex) then
               if ((.not.extrap)
     &            .and.(chr(i+1,j,k).ne.ex)
     &            .and.(chr(i+2,j,k).ne.ex)
     &            .and.(chr(i+3,j,k).ne.ex)) then
                   ddgb_J_tx=-1/dt/dx
               else if ((chr(i+1,j,k).ne.ex)
     &                 .and.(chr(i+2,j,k).ne.ex)) then
                   ddgb_J_tx=-3/4/dt/dx
               else if (chr(i+1,j,k).ne.ex) then
                   ddgb_J_tx=-1/2/dt/dx
!                  write(*,*) 'g_evo_opt: warning ... i=1 first order'
!                  write(*,*) '    i,j,Nx,Ny,dx=',i,j,Nx,Ny,dx
               else
                   ddgb_J_tx=0
                   return
               end if
         else 
               if (chr(i-2,j,k).ne.ex) then
                   ddgb_J_tx=3/4/dt/dx
               else
                   ddgb_J_tx=1/2/dt/dx
!                  write(*,*) 'g_evo_opt: warning ... i=1 first order'
!                  write(*,*) '    i,j,Nx,Ny,dx=',i,j,Nx,Ny,dx
               end if
         end if
        else if ((i.ge.4).and.(i.le.(Nx-3))) then
         if ((chr(i-1,j,k).ne.ex).and.(chr(i+1,j,k).ne.ex)) then
                   ddgb_J_tx=0
         else if (chr(i-1,j,k).eq.ex) then
               if ((.not.extrap)
     &            .and.(chr(i+1,j,k).ne.ex)
     &            .and.(chr(i+2,j,k).ne.ex)
     &            .and.(chr(i+3,j,k).ne.ex)) then
                   ddgb_J_tx=-1/dt/dx
               else if ((chr(i+1,j,k).ne.ex)
     &                 .and.(chr(i+2,j,k).ne.ex)) then
                   ddgb_J_tx=-3/4/dt/dx
               else if (chr(i+1,j,k).ne.ex) then
                   ddgb_J_tx=-1/2/dt/dx
!                  write(*,*) 'g_evo_opt: warning ... i=1 first order'
!                  write(*,*) '    i,j,Nx,Ny,dx=',i,j,Nx,Ny,dx
               else
                   ddgb_J_tx=0
                   return
               end if
         else
               if ((.not.extrap)
     &            .and.(chr(i-3,j,k).ne.ex)
     &            .and.(chr(i-2,j,k).ne.ex)) then
                   ddgb_J_tx=1/dt/dx
               else if (chr(i-2,j,k).ne.ex) then
                   ddgb_J_tx=3/4/dt/dx
               else
                   ddgb_J_tx=1/2/dt/dx
!                  write(*,*) 'g_evo_opt: warning ... i=1 first order'
!                  write(*,*) '    i,j,Nx,Ny,dx=',i,j,Nx,Ny,dx
               end if
         end if
        else if (i.eq.(Nx-2)) then
         if ((chr(i+1,j,k).ne.ex).and.(chr(i-1,j,k).ne.ex)) then
                   ddgb_J_tx=0
         else if (chr(i+1,j,k).eq.ex) then
               if ((.not.extrap)
     &            .and.(chr(i-1,j,k).ne.ex)
     &            .and.(chr(i-2,j,k).ne.ex)
     &            .and.(chr(i-3,j,k).ne.ex)) then
                   ddgb_J_tx=1/dt/dx
               else if ((chr(i-1,j,k).ne.ex)
     &                 .and.(chr(i-2,j,k).ne.ex)) then
                   ddgb_J_tx=3/4/dt/dx
               else if (chr(i-1,j,k).ne.ex) then
                   ddgb_J_tx=1/2/dt/dx
!                  write(*,*) 'g_evo_opt: warning ... i=1 first order'
!                  write(*,*) '    i,j,Nx,Ny,dx=',i,j,Nx,Ny,dx
               else
                   ddgb_J_tx=0
                   return
               end if
         else 
               if (chr(i+2,j,k).ne.ex) then
                   ddgb_J_tx=-3/4/dt/dx
               else
                   ddgb_J_tx=-1/2/dt/dx
!                  write(*,*) 'g_evo_opt: warning ... i=1 first order'
!                  write(*,*) '    i,j,Nx,Ny,dx=',i,j,Nx,Ny,dx
               end if
         end if
        else if (i.eq.(Nx-1)) then
         if ((chr(i+1,j,k).ne.ex).and.(chr(i-1,j,k).ne.ex)) then
                   ddgb_J_tx=0
         else if (chr(i+1,j,k).eq.ex) then
               if ((.not.extrap)
     &            .and.(chr(i-1,j,k).ne.ex)
     &            .and.(chr(i-2,j,k).ne.ex)
     &            .and.(chr(i-3,j,k).ne.ex)) then
                   ddgb_J_tx=1/dt/dx
               else if ((chr(i-1,j,k).ne.ex)
     &                 .and.(chr(i-2,j,k).ne.ex)) then
                   ddgb_J_tx=3/4/dt/dx
               else if (chr(i-1,j,k).ne.ex) then
                   ddgb_J_tx=1/2/dt/dx
!                  write(*,*) 'g_evo_opt: warning ... i=1 first order'
!                  write(*,*) '    i,j,Nx,Ny,dx=',i,j,Nx,Ny,dx
               else
                   ddgb_J_tx=0
                   return
               end if
         else
                   ddgb_J_tx=-1/2/dt/dx
!                  write(*,*) 'g_evo_opt: warning ... i=1 first order'
!                  write(*,*) '    i,j,Nx,Ny,dx=',i,j,Nx,Ny,dx
         end if
        else if (i.eq.Nx) then
               if ((.not.extrap)
     &            .and.(chr(i-1,j,k).ne.ex)
     &            .and.(chr(i-2,j,k).ne.ex)
     &            .and.(chr(i-3,j,k).ne.ex)) then
                   ddgb_J_tx=1/dt/dx
               else if ((chr(i-1,j,k).ne.ex)
     &                 .and.(chr(i-2,j,k).ne.ex)) then
                   ddgb_J_tx=3/4/dt/dx
               else if (chr(i-1,j,k).ne.ex) then
                   ddgb_J_tx=1/2/dt/dx
!                  write(*,*) 'g_evo_opt: warning ... i=1 first order'
!                  write(*,*) '    i,j,Nx,Ny,dx=',i,j,Nx,Ny,dx
               else
                   ddgb_J_tx=0
                   return
               end if
        end if



!j
        if (j.eq.1) then
               if ((.not.extrap)
     &            .and.(chr(i,j+1,k).ne.ex)
     &            .and.(chr(i,j+2,k).ne.ex)
     &            .and.(chr(i,j+3,k).ne.ex)) then
                   ddgb_J_ty=-1/dt/dy
               else if ((chr(i,j+1,k).ne.ex
     &                 .and.chr(i,j+2,k).ne.ex)) then
                   ddgb_J_ty=-3/4/dt/dy
               else if (chr(i,j+1,k).ne.ex) then
                   ddgb_J_ty=-1/2/dt/dy
!                  write(*,*) 'g_evo_opt: warning ... j=1 first order'
!                  write(*,*) '    i,j,Nx,Ny,dy=',i,j,Nx,Ny,dy
               else
                   ddgb_J_ty=0
                   return
               end if
        else if (j.eq.2) then
         if ((chr(i,j-1,k).ne.ex).and.(chr(i,j+1,k).ne.ex)) then
                   ddgb_J_ty=0
         else if (chr(i,j-1,k).eq.ex) then
               if ((.not.extrap)
     &            .and.(chr(i,j+1,k).ne.ex)
     &            .and.(chr(i,j+2,k).ne.ex)
     &            .and.(chr(i,j+3,k).ne.ex)) then
                   ddgb_J_ty=-1/dt/dy
               else if ((chr(i,j+1,k).ne.ex)
     &                 .and.(chr(i,j+2,k).ne.ex)) then
                   ddgb_J_ty=-3/4/dt/dy
               else if (chr(i,j+1,k).ne.ex) then
                   ddgb_J_ty=-1/2/dt/dy
!                  write(*,*) 'g_evo_opt: warning ... j=1 first order'
!                  write(*,*) '    i,j,Nx,Ny,dy=',i,j,Nx,Ny,dy
               else
                   ddgb_J_ty=0
                   return
               end if
         else   !this is the case where (i,j-1,k) is not excised and (i,j+1,k) is excised 
                   ddgb_J_ty=1/2/dt/dy
!                  write(*,*) 'g_evo_opt: warning ... j=1 first order'
!                  write(*,*) '    i,j,Nx,Ny,dy=',i,j,Nx,Ny,dy
         end if
        else if (j.eq.3) then
         if ((chr(i,j-1,k).ne.ex).and.(chr(i,j+1,k).ne.ex)) then
                   ddgb_J_ty=0
         else if (chr(i,j-1,k).eq.ex) then
               if ((.not.extrap)
     &            .and.(chr(i,j+1,k).ne.ex)
     &            .and.(chr(i,j+2,k).ne.ex)
     &            .and.(chr(i,j+3,k).ne.ex)) then
                   ddgb_J_ty=-1/dt/dy
               else if ((chr(i,j+1,k).ne.ex)
     &                 .and.(chr(i,j+2,k).ne.ex)) then
                   ddgb_J_ty=-3/4/dt/dy
               else if (chr(i,j+1,k).ne.ex) then
                   ddgb_J_ty=-1/2/dt/dy
!                  write(*,*) 'g_evo_opt: warning ... j=1 first order'
!                  write(*,*) '    i,j,Nx,Ny,dy=',i,j,Nx,Ny,dy
               else
                   ddgb_J_ty=0
                   return
               end if
         else 
               if (chr(i,j-2,k).ne.ex) then
                   ddgb_J_ty=3/4/dt/dy
               else
                   ddgb_J_ty=1/2/dt/dy
!                  write(*,*) 'g_evo_opt: warning ... j=1 first order'
!                  write(*,*) '    i,j,Nx,Ny,dy=',i,j,Nx,Ny,dy
               end if
         end if
        else if ((j.ge.4).and.(j.le.(Ny-3))) then
         if ((chr(i,j-1,k).ne.ex).and.(chr(i,j+1,k).ne.ex)) then
                   ddgb_J_ty=0
         else if (chr(i,j-1,k).eq.ex) then
               if ((.not.extrap)
     &            .and.(chr(i,j+1,k).ne.ex)
     &            .and.(chr(i,j+2,k).ne.ex)
     &            .and.(chr(i,j+3,k).ne.ex)) then
                   ddgb_J_ty=-1/dt/dy
               else if ((chr(i,j+1,k).ne.ex)
     &                 .and.(chr(i,j+2,k).ne.ex)) then
                   ddgb_J_ty=-3/4/dt/dy
               else if (chr(i,j+1,k).ne.ex) then
                   ddgb_J_ty=-1/2/dt/dy
!                  write(*,*) 'g_evo_opt: warning ... j=1 first order'
!                  write(*,*) '    i,j,Nx,Ny,dy=',i,j,Nx,Ny,dy
               else
                   ddgb_J_ty=0
                   return
               end if
         else
               if ((.not.extrap)
     &            .and.(chr(i,j-3,k).ne.ex)
     &            .and.(chr(i,j-2,k).ne.ex)) then
                   ddgb_J_ty=1/dt/dy
               else if (chr(i,j-2,k).ne.ex) then
                   ddgb_J_ty=3/4/dt/dy
               else
                   ddgb_J_ty=1/2/dt/dy
!                  write(*,*) 'g_evo_opt: warning ... j=1 first order'
!                  write(*,*) '    i,j,Nx,Ny,dy=',i,j,Nx,Ny,dy
               end if
         end if
        else if (j.eq.(Ny-2)) then
         if ((chr(i,j+1,k).ne.ex).and.(chr(i,j-1,k).ne.ex)) then
                   ddgb_J_ty=0
         else if (chr(i,j+1,k).eq.ex) then
               if ((.not.extrap)
     &            .and.(chr(i,j-1,k).ne.ex)
     &            .and.(chr(i,j-2,k).ne.ex)
     &            .and.(chr(i,j-3,k).ne.ex)) then
                   ddgb_J_ty=1/dt/dy
               else if ((chr(i,j-1,k).ne.ex)
     &                 .and.(chr(i,j-2,k).ne.ex)) then
                   ddgb_J_ty=3/4/dt/dy
               else if (chr(i,j-1,k).ne.ex) then
                   ddgb_J_ty=1/2/dt/dy
!                  write(*,*) 'g_evo_opt: warning ... j=1 first order'
!                  write(*,*) '    i,j,Nx,Ny,dy=',i,j,Nx,Ny,dy
               else
                   ddgb_J_ty=0
                   return
               end if
         else 
               if (chr(i,j+2,k).ne.ex) then
                   ddgb_J_ty=-3/4/dt/dy
               else
                   ddgb_J_ty=-1/2/dt/dy
!                  write(*,*) 'g_evo_opt: warning ... j=1 first order'
!                  write(*,*) '    i,j,Nx,Ny,dy=',i,j,Nx,Ny,dy
               end if
         end if
        else if (j.eq.(Ny-1)) then
         if ((chr(i,j+1,k).ne.ex).and.(chr(i,j-1,k).ne.ex)) then
                   ddgb_J_ty=0
         else if (chr(i,j+1,k).eq.ex) then
               if ((.not.extrap)
     &            .and.(chr(i,j-1,k).ne.ex)
     &            .and.(chr(i,j-2,k).ne.ex)
     &            .and.(chr(i,j-3,k).ne.ex)) then
                   ddgb_J_ty=1/dt/dy
               else if ((chr(i,j-1,k).ne.ex)
     &                 .and.(chr(i,j-2,k).ne.ex)) then
                   ddgb_J_ty=3/4/dt/dy
               else if (chr(i,j-1,k).ne.ex) then
                   ddgb_J_ty=1/2/dt/dy
!                  write(*,*) 'g_evo_opt: warning ... j=1 first order'
!                  write(*,*) '    i,j,Nx,Ny,dy=',i,j,Nx,Ny,dy
               else
                   ddgb_J_ty=0
                   return
               end if
         else
                   ddgb_J_ty=-1/2/dt/dy
!                  write(*,*) 'g_evo_opt: warning ... j=1 first order'
!                  write(*,*) '    i,j,Nx,Ny,dy=',i,j,Nx,Ny,dy
         end if
        else if (j.eq.Ny) then
               if ((.not.extrap)
     &            .and.(chr(i,j-1,k).ne.ex)
     &            .and.(chr(i,j-2,k).ne.ex)
     &            .and.(chr(i,j-3,k).ne.ex)) then
                   ddgb_J_ty=1/dt/dy
               else if ((chr(i,j-1,k).ne.ex)
     &                 .and.(chr(i,j-2,k).ne.ex)) then
                   ddgb_J_ty=3/4/dt/dy
               else if (chr(i,j-1,k).ne.ex) then
                   ddgb_J_ty=1/2/dt/dy
!                  write(*,*) 'g_evo_opt: warning ... j=1 first order'
!                  write(*,*) '    i,j,Nx,Ny,dy=',i,j,Nx,Ny,dy
               else
                   ddgb_J_ty=0
                   return
               end if
        end if

!k
!        ddgb_J_tz=0

        if (k.eq.1) then
               if ((.not.extrap)
     &            .and.(chr(i,j,k+1).ne.ex)
     &            .and.(chr(i,j,k+2).ne.ex)
     &            .and.(chr(i,j,k+3).ne.ex)) then
                   ddgb_J_tz=-1/dt/dz
               else if ((chr(i,j,k+1).ne.ex
     &                 .and.chr(i,j,k+2).ne.ex)) then
                   ddgb_J_tz=-3/4/dt/dz
               else if (chr(i,j,k+1).ne.ex) then
                   ddgb_J_tz=-1/2/dt/dz
!                  write(*,*) 'g_evo_opt: warning ... k=1 first order'
!                  write(*,*) '    i,j,k,Nx,Ny,Nz,dz=',i,j,k,Nx,Ny,Nz,dz
               else
                   ddgb_J_tz=0
                   return
               end if
        else if (k.eq.2) then
         if ((chr(i,j,k-1).ne.ex).and.(chr(i,j,k+1).ne.ex)) then
                   ddgb_J_tz=0
         else if (chr(i,j,k-1).eq.ex) then
               if ((.not.extrap)
     &            .and.(chr(i,j,k+1).ne.ex)
     &            .and.(chr(i,j,k+2).ne.ex)
     &            .and.(chr(i,j,k+3).ne.ex)) then
                   ddgb_J_tz=-1/dt/dz
               else if ((chr(i,j,k+1).ne.ex)
     &                 .and.(chr(i,j,k+2).ne.ex)) then
                   ddgb_J_tz=-3/4/dt/dz
               else if (chr(i,j,k+1).ne.ex) then
                   ddgb_J_tz=-1/2/dt/dz
!                  write(*,*) 'g_evo_opt: warning ... k=1 first order'
!                  write(*,*) '    i,j,k,Nx,Ny,Nz,dz=',i,j,k,Nx,Ny,Nz,dz
               else
                   ddgb_J_tz=0
                   return
               end if
         else   !this is the case where (i,j,k-1) is not excised and (i,j,k+1) is excised 
                   ddgb_J_tz=1/2/dt/dz
!                  write(*,*) 'g_evo_opt: warning ... k=1 first order'
!                  write(*,*) '    i,j,k,Nx,Ny,Nz,dz=',i,j,k,Nx,Ny,Nz,dz
         end if
        else if (k.eq.3) then
         if ((chr(i,j,k-1).ne.ex).and.(chr(i,j,k+1).ne.ex)) then
                   ddgb_J_tz=0
         else if (chr(i,j,k-1).eq.ex) then
               if ((.not.extrap)
     &            .and.(chr(i,j,k+1).ne.ex)
     &            .and.(chr(i,j,k+2).ne.ex)
     &            .and.(chr(i,j,k+3).ne.ex)) then
                   ddgb_J_tz=-1/dt/dz
               else if ((chr(i,j,k+1).ne.ex)
     &                 .and.(chr(i,j,k+2).ne.ex)) then
                   ddgb_J_tz=-3/4/dt/dz
               else if (chr(i,j,k+1).ne.ex) then
                   ddgb_J_tz=-1/2/dt/dz
!                  write(*,*) 'g_evo_opt: warning ... k=1 first order'
!                  write(*,*) '    i,j,k,Nx,Ny,Nz,dz=',i,j,k,Nx,Ny,Nz,dz
               else
                   ddgb_J_tz=0
                   return
               end if
         else
               if (chr(i,j,k-2).ne.ex) then
                   ddgb_J_tz=3/4/dt/dz
               else
                   ddgb_J_tz=1/2/dt/dz
!                  write(*,*) 'g_evo_opt: warning ... k=1 first order'
!                  write(*,*) '    i,j,k,Nx,Ny,Nz,dz=',i,j,k,Nx,Ny,Nz,dz
               end if
         end if
        else if ((k.ge.4).and.(k.le.(Nz-3))) then
         if ((chr(i,j,k-1).ne.ex).and.(chr(i,j,k+1).ne.ex)) then
                   ddgb_J_tz=0
         else if (chr(i,j,k-1).eq.ex) then
               if ((.not.extrap)
     &            .and.(chr(i,j,k+1).ne.ex)
     &            .and.(chr(i,j,k+2).ne.ex)
     &            .and.(chr(i,j,k+3).ne.ex)) then
                   ddgb_J_tz=-1/dt/dz
               else if ((chr(i,j,k+1).ne.ex)
     &                 .and.(chr(i,j,k+2).ne.ex)) then
                   ddgb_J_tz=-3/4/dt/dz
               else if (chr(i,j,k+1).ne.ex) then
                   ddgb_J_tz=-1/2/dt/dz
!                  write(*,*) 'g_evo_opt: warning ... k=1 first order'
!                  write(*,*) '    i,j,k,Nx,Ny,Nz,dz=',i,j,k,Nx,Ny,Nz,dz
               else
                   ddgb_J_tz=0
                   return
               end if
         else
               if ((.not.extrap)
     &            .and.(chr(i,j,k-3).ne.ex)
     &            .and.(chr(i,j,k-2).ne.ex)) then
                   ddgb_J_tz=1/dt/dz
               else if (chr(i,j,k-2).ne.ex) then
                   ddgb_J_tz=3/4/dt/dz
               else
                   ddgb_J_tz=1/2/dt/dz
!                  write(*,*) 'g_evo_opt: warning ... k=1 first order'
!                  write(*,*) '    i,j,k,Nx,Ny,Nz,dz=',i,j,k,Nx,Ny,Nz,dz
               end if
         end if
        else if (k.eq.(Nz-2)) then
         if ((chr(i,j,k+1).ne.ex).and.(chr(i,j,k-1).ne.ex)) then
                   ddgb_J_tz=0
         else if (chr(i,j,k+1).eq.ex) then
               if ((.not.extrap)
     &            .and.(chr(i,j,k-1).ne.ex)
     &            .and.(chr(i,j,k-2).ne.ex)
     &            .and.(chr(i,j,k-3).ne.ex)) then
                   ddgb_J_tz=1/dt/dz
               else if ((chr(i,j,k-1).ne.ex)
     &                 .and.(chr(i,j,k-2).ne.ex)) then
                   ddgb_J_tz=3/4/dt/dz
               else if (chr(i,j,k-1).ne.ex) then
                   ddgb_J_tz=1/2/dt/dz
!                  write(*,*) 'g_evo_opt: warning ... k=1 first order'
!                  write(*,*) '    i,j,k,Nx,Ny,Nz,dz=',i,j,k,Nx,Ny,Nz,dz
               else
                   ddgb_J_tz=0
                   return
               end if
         else
               if (chr(i,j,k+2).ne.ex) then
                   ddgb_J_tz=-3/4/dt/dz
               else
                   ddgb_J_tz=-1/2/dt/dz
!                  write(*,*) 'g_evo_opt: warning ... k=1 first order'
!                  write(*,*) '    i,j,k,Nx,Ny,Nz,dz=',i,j,k,Nx,Ny,Nz,dz
               end if
         end if
        else if (k.eq.(Nz-1)) then
         if ((chr(i,j,k+1).ne.ex).and.(chr(i,j,k-1).ne.ex)) then
                   ddgb_J_tz=0
         else if (chr(i,j,k+1).eq.ex) then
               if ((.not.extrap)
     &            .and.(chr(i,j,k-1).ne.ex)
     &            .and.(chr(i,j,k-2).ne.ex)
     &            .and.(chr(i,j,k-3).ne.ex)) then
                   ddgb_J_tz=1/dt/dz
               else if ((chr(i,j,k-1).ne.ex)
     &                 .and.(chr(i,j,k-2).ne.ex)) then
                   ddgb_J_tz=3/4/dt/dz
               else if (chr(i,j,k-1).ne.ex) then
                   ddgb_J_tz=1/2/dt/dz
!                  write(*,*) 'g_evo_opt: warning ... k=1 first order'
!                  write(*,*) '    i,j,k,Nx,Ny,Nz,dz=',i,j,k,Nx,Ny,Nz,dz
               else
                   ddgb_J_tz=0
                   return
               end if
         else
                   ddgb_J_tz=-1/2/dt/dz
!                  write(*,*) 'g_evo_opt: warning ... k=1 first order'
!                  write(*,*) '    i,j,k,Nx,Ny,Nz,dz=',i,j,k,Nx,Ny,Nz,dz
         end if
        else if (k.eq.Nz) then
               if ((.not.extrap)
     &            .and.(chr(i,j,k-1).ne.ex)
     &            .and.(chr(i,j,k-2).ne.ex)
     &            .and.(chr(i,j,k-3).ne.ex)) then
                   ddgb_J_tz=1/dt/dz
               else if ((chr(i,j,k-1).ne.ex)
     &                 .and.(chr(i,j,k-2).ne.ex)) then
                   ddgb_J_tz=3/4/dt/dz
               else if (chr(i,j,k-1).ne.ex) then
                   ddgb_J_tz=1/2/dt/dz
!                  write(*,*) 'g_evo_opt: warning ... k=1 first order'
!                  write(*,*) '    i,j,k,Nx,Ny,Nz,dz=',i,j,k,Nx,Ny,Nz,dz
               else
                   ddgb_J_tz=0
                   return
               end if
        end if




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





!!!!!OLD VERSION
!                if (i.eq.1.or.(chr(i-1,j,k).eq.ex)) then
!                   if ((.not.extrap)
!     &                 .and.(i.le.(Nx-3))
!     &                 .and.((chr(i+1,j,k).ne.ex
!     &                 .and.chr(i+2,j,k).ne.ex
!     &                 .and.chr(i+3,j,k).ne.ex))) then
!                      ddgb_J_tx=-1/dt/dx
!                   else if (i.le.(Nx-2)
!     &                      .and.((chr(i+1,j,k).ne.ex
!     &                      .and.chr(i+2,j,k).ne.ex))) then
!                      ddgb_J_tx=-3/4/dt/dx
!                   else if (i.le.(Nx-1).and.chr(i+1,j,k).ne.ex) then
!                      ddgb_J_tx=-1/2/dt/dx
!                   else
!                      write(*,*) 'g_evo_opt: error in chr stencil (A)'
!                      write(*,*) '    i,j,Nx,Ny,dx=',i,j,Nx,Ny,dx
!                      write(*,*) '    (first error only)'
!                      ddgb_J_tx=0
!                   end if
!                else if (i.eq.Nx.or.(chr(i+1,j,k).eq.ex)) then
!                   if ((.not.extrap)
!     &                 .and.(i.ge.4)
!     &                 .and.((chr(i-1,j,k).ne.ex
!     &                 .and.chr(i-2,j,k).ne.ex
!     &                 .and.chr(i-3,j,k).ne.ex))) then
!                      ddgb_J_tx=1/dt/dx
!                   else if (i.ge.3
!     &                      .and.((chr(i-1,j,k).ne.ex
!     &                      .and.chr(i-2,j,k).ne.ex))) then
!                      ddgb_J_tx=3/4/dt/dx
!                   else if (i.ge.2.and.chr(i-1,j,k).ne.ex) then
!                      ddgb_J_tx=1/2/dt/dx
!                   else
!                      write(*,*) 'g_evo_opt: error in chr stencil (B)'
!                      write(*,*) '    i,j,Nx,Ny,dx=',i,j,Nx,Ny,dx
!                      write(*,*) '    (first error only)'
!                      ddgb_J_tx=0
!                   end if
!                else
!                   if ((chr(i+1,j,k).ne.ex.and.chr(i-1,j,k).ne.ex)) then
!                      ddgb_J_tx=0
!                   else
!                      write(*,*) 'g_evo_opt: error in chr stencil (C)'
!                      write(*,*) '    i,j,Nx,Ny,dx=',i,j,Nx,Ny,dx
!                      write(*,*) '    (first error only)'
!                      ddgb_J_tx=0
!                   end if
!                end if
!
!                if ((j.eq.1).or.(chr(i,j-1,k).eq.ex)) then
!                   if ((.not.extrap)
!     &                 .and.(j.le.(Ny-3))
!     &                 .and.((chr(i,j+1,k).ne.ex
!     &                 .and.chr(i,j+2,k).ne.ex
!     &                 .and.chr(i,j+3,k).ne.ex))) then
!                      ddgb_J_ty=-1/dt/dy
!                   else if (j.le.(Ny-2).and.((chr(i,j+1,k).ne.ex
!     &                      .and.chr(i,j+2,k).ne.ex))) then
!                      ddgb_J_ty=-3/4/dt/dy              
!                   else if (j.le.(Ny-1).and.chr(i,j+1,k).ne.ex) then
!                      ddgb_J_ty=-1/2/dt/dy
!                   else
!                      write(*,*) 'g_evo_opt: error in chr stencil (D)'
!                      write(*,*) '    i,j,Nx,Ny,dy=',i,j,Nx,Ny,dy
!                      write(*,*) '    (first error only)'
!                      ddgb_J_ty=0
!                   end if
!                else if ((j.eq.Ny).or.(chr(i,j+1,k).eq.ex)) then
!                   if ((.not.extrap)
!     &                 .and.(j.ge.4)
!     &                 .and.((chr(i,j-1,k).ne.ex
!     &                 .and.chr(i,j-2,k).ne.ex
!     &                 .and.chr(i,j-3,k).ne.ex))) then
!                      ddgb_J_ty=1/dt/dy
!                   else if (j.ge.3.and.((chr(i,j-1,k).ne.ex
!     &                      .and.chr(i,j-2,k).ne.ex))) then
!                      ddgb_J_ty=3/4/dt/dy
!                   else if (j.ge.2.and.chr(i,j-1,k).ne.ex) then
!                      ddgb_J_ty=1/2/dt/dy
!                   else
!                      write(*,*) 'g_evo_opt: error in chr stencil (E)'
!                      write(*,*) '    i,j,Nx,Ny,dy=',i,j,Nx,Ny,dy
!                      write(*,*) '    (first error only)'
!                      ddgb_J_ty=0
!                   end if
!                else
!                   if ((chr(i,j+1,k).ne.ex.and.chr(i,j-1,k).ne.ex)) then
!                      ddgb_J_ty=0
!                   else
!                      write(*,*) 'g_evo_opt: error in chr stencil (F)'
!                      write(*,*) '    i,j,Nx,Ny,dy=',i,j,Nx,Ny,dy
!                      write(*,*) '    (first error only)'
!                      ddgb_J_ty=0
!                   end if
!                end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                efe_J(1,1)=    -0.5d0*(
     &                            g0_uu(1,1)*ddgb_J
     &                            +2*g0_uu(1,2)*ddgb_J_tx
     &                            +2*g0_uu(1,3)*ddgb_J_ty
     &                            +2*g0_uu(1,4)*ddgb_J_tz
     &                                )
     &                      
     &                         -0.5d0*(
     &                            -dgb_J*
     &                            (g0_uu(1,1)*g0_uu(1,1)*g0_ll_x(1,1,1)+
     &                             g0_uu(1,1)*g0_uu(2,1)*g0_ll_x(1,1,2)+
     &                             g0_uu(1,1)*g0_uu(3,1)*g0_ll_x(1,1,3)+
     &                             g0_uu(1,1)*g0_uu(4,1)*g0_ll_x(1,1,4)+
     &                             g0_uu(2,1)*g0_uu(1,1)*g0_ll_x(1,2,1)+
     &                             g0_uu(2,1)*g0_uu(2,1)*g0_ll_x(1,2,2)+
     &                             g0_uu(2,1)*g0_uu(3,1)*g0_ll_x(1,2,3)+
     &                             g0_uu(2,1)*g0_uu(4,1)*g0_ll_x(1,2,4)+
     &                             g0_uu(3,1)*g0_uu(1,1)*g0_ll_x(1,3,1)+
     &                             g0_uu(3,1)*g0_uu(2,1)*g0_ll_x(1,3,2)+
     &                             g0_uu(3,1)*g0_uu(3,1)*g0_ll_x(1,3,3)+
     &                             g0_uu(3,1)*g0_uu(4,1)*g0_ll_x(1,3,4)+
     &                             g0_uu(4,1)*g0_uu(1,1)*g0_ll_x(1,4,1)+
     &                             g0_uu(4,1)*g0_uu(2,1)*g0_ll_x(1,4,2)+
     &                             g0_uu(4,1)*g0_uu(3,1)*g0_ll_x(1,4,3)+
     &                             g0_uu(4,1)*g0_uu(4,1)*g0_ll_x(1,4,4))
     &                            +dgb_J*
     &                            (g0_uu_x(1,1,1))
     &                                )*2
     &
     &                         +      (
     &                            0.5d0*dgb_J*
     &                            ((Hads_l(1)+A_l(1))*g0_uu(1,1)+
     &                             (Hads_l(2)+A_l(2))*g0_uu(2,1)+
     &                             (Hads_l(3)+A_l(3))*g0_uu(3,1)+
     &                             (Hads_l(4)+A_l(4))*g0_uu(4,1))
     &                                )  
     &                      
     &                         -      (
     &                            0.25d0*dgb_J*
     &                            (cuuuu(1,1,1,1)*dlll(1,1,1)+
     &                             cuuuu(1,1,1,2)*dlll(1,1,2)+
     &                             cuuuu(1,1,1,3)*dlll(1,1,3)+
     &                             cuuuu(1,1,1,4)*dlll(1,1,4)+
     &                             cuuuu(2,1,1,1)*dlll(2,1,1)+
     &                             cuuuu(2,1,1,2)*dlll(2,1,2)+
     &                             cuuuu(2,1,1,3)*dlll(2,1,3)+
     &                             cuuuu(2,1,1,4)*dlll(2,1,4)+
     &                             cuuuu(3,1,1,1)*dlll(3,1,1)+
     &                             cuuuu(3,1,1,2)*dlll(3,1,2)+
     &                             cuuuu(3,1,1,3)*dlll(3,1,3)+
     &                             cuuuu(3,1,1,4)*dlll(3,1,4)+
     &                             cuuuu(4,1,1,1)*dlll(4,1,1)+
     &                             cuuuu(4,1,1,2)*dlll(4,1,2)+
     &                             cuuuu(4,1,1,3)*dlll(4,1,3)+
     &                             cuuuu(4,1,1,4)*dlll(4,1,4)
     &                           +
     &                             cuuuu(1,1,1,1)*dlll(1,1,1)+
     &                             cuuuu(1,1,2,1)*dlll(2,1,1)+
     &                             cuuuu(1,1,3,1)*dlll(3,1,1)+
     &                             cuuuu(1,1,4,1)*dlll(4,1,1)+
     &                             cuuuu(1,2,1,1)*dlll(1,1,2)+
     &                             cuuuu(1,2,2,1)*dlll(2,1,2)+
     &                             cuuuu(1,2,3,1)*dlll(3,1,2)+
     &                             cuuuu(1,2,4,1)*dlll(4,1,2)+
     &                             cuuuu(1,3,1,1)*dlll(1,1,3)+
     &                             cuuuu(1,3,2,1)*dlll(2,1,3)+
     &                             cuuuu(1,3,3,1)*dlll(3,1,3)+
     &                             cuuuu(1,3,4,1)*dlll(4,1,3)+
     &                             cuuuu(1,4,1,1)*dlll(1,1,4)+
     &                             cuuuu(1,4,2,1)*dlll(2,1,4)+
     &                             cuuuu(1,4,3,1)*dlll(3,1,4)+
     &                             cuuuu(1,4,4,1)*dlll(4,1,4))
     &                                )

                do b=2,3
                  efe_J(1,b)=  -0.5d0*(
     &                            g0_uu(1,1)*ddgb_J
     &                            +2*g0_uu(1,2)*ddgb_J_tx
     &                            +2*g0_uu(1,3)*ddgb_J_ty
     &                            +2*g0_uu(1,4)*ddgb_J_tz
     &                                )
     &                         
     &                         -0.5d0*(
     &                            -dgb_J*
     &                            (g0_uu(1,1)*g0_uu(1,b)*g0_ll_x(b,1,1)+
     &                             g0_uu(1,1)*g0_uu(2,b)*g0_ll_x(b,1,2)+
     &                             g0_uu(1,1)*g0_uu(3,b)*g0_ll_x(b,1,3)+
     &                             g0_uu(1,1)*g0_uu(4,b)*g0_ll_x(b,1,4)+
     &                             g0_uu(2,1)*g0_uu(1,b)*g0_ll_x(b,2,1)+
     &                             g0_uu(2,1)*g0_uu(2,b)*g0_ll_x(b,2,2)+
     &                             g0_uu(2,1)*g0_uu(3,b)*g0_ll_x(b,2,3)+
     &                             g0_uu(2,1)*g0_uu(4,b)*g0_ll_x(b,2,4)+
     &                             g0_uu(3,1)*g0_uu(1,b)*g0_ll_x(b,3,1)+
     &                             g0_uu(3,1)*g0_uu(2,b)*g0_ll_x(b,3,2)+
     &                             g0_uu(3,1)*g0_uu(3,b)*g0_ll_x(b,3,3)+
     &                             g0_uu(3,1)*g0_uu(4,b)*g0_ll_x(b,3,4)+
     &                             g0_uu(4,1)*g0_uu(1,b)*g0_ll_x(b,4,1)+
     &                             g0_uu(4,1)*g0_uu(2,b)*g0_ll_x(b,4,2)+
     &                             g0_uu(4,1)*g0_uu(3,b)*g0_ll_x(b,4,3)+
     &                             g0_uu(4,1)*g0_uu(4,b)*g0_ll_x(b,4,4)+
     &
     &                             g0_uu(1,b)*g0_uu(1,1)*g0_ll_x(b,1,1)+
     &                             g0_uu(1,b)*g0_uu(2,1)*g0_ll_x(b,1,2)+
     &                             g0_uu(1,b)*g0_uu(3,1)*g0_ll_x(b,1,3)+
     &                             g0_uu(1,b)*g0_uu(4,1)*g0_ll_x(b,1,4)+
     &                             g0_uu(2,b)*g0_uu(1,1)*g0_ll_x(b,2,1)+
     &                             g0_uu(2,b)*g0_uu(2,1)*g0_ll_x(b,2,2)+
     &                             g0_uu(2,b)*g0_uu(3,1)*g0_ll_x(b,2,3)+
     &                             g0_uu(2,b)*g0_uu(4,1)*g0_ll_x(b,2,4)+
     &                             g0_uu(3,b)*g0_uu(1,1)*g0_ll_x(b,3,1)+
     &                             g0_uu(3,b)*g0_uu(2,1)*g0_ll_x(b,3,2)+
     &                             g0_uu(3,b)*g0_uu(3,1)*g0_ll_x(b,3,3)+
     &                             g0_uu(3,b)*g0_uu(4,1)*g0_ll_x(b,3,4)+
     &                             g0_uu(4,b)*g0_uu(1,1)*g0_ll_x(b,4,1)+
     &                             g0_uu(4,b)*g0_uu(2,1)*g0_ll_x(b,4,2)+
     &                             g0_uu(4,b)*g0_uu(3,1)*g0_ll_x(b,4,3)+
     &                             g0_uu(4,b)*g0_uu(4,1)*g0_ll_x(b,4,4))
     &                            +dgb_J*
     &                            (g0_uu_x(1,1,1))
     &                                )
     &                      
     &                         -0.5d0*(
     &                            dgb_J*
     &                            (g0_uu_x(b,1,b))
     &                                )
     &
     &                         -      (
     &                            0.5d0*dgb_J*
     &                            (cuuuu(1,1,1,b)*dlll(1,b,1)+
     &                             cuuuu(1,1,2,b)*dlll(2,b,1)+
     &                             cuuuu(1,1,3,b)*dlll(3,b,1)+
     &                             cuuuu(1,1,4,b)*dlll(4,b,1)+
     &                             cuuuu(1,2,1,b)*dlll(1,b,2)+
     &                             cuuuu(1,2,2,b)*dlll(2,b,2)+
     &                             cuuuu(1,2,3,b)*dlll(3,b,2)+
     &                             cuuuu(1,2,4,b)*dlll(4,b,2)+
     &                             cuuuu(1,3,1,b)*dlll(1,b,3)+
     &                             cuuuu(1,3,2,b)*dlll(2,b,3)+
     &                             cuuuu(1,3,3,b)*dlll(3,b,3)+
     &                             cuuuu(1,3,4,b)*dlll(4,b,3)+
     &                             cuuuu(1,4,1,b)*dlll(1,b,4)+
     &                             cuuuu(1,4,2,b)*dlll(2,b,4)+
     &                             cuuuu(1,4,3,b)*dlll(3,b,4)+
     &                             cuuuu(1,4,4,b)*dlll(4,b,4))
     &                                )
                end do

                   efe_J(1,4)=  -0.5d0*(
     &                            g0_uu(1,1)*ddgb_J
     &                            +2*g0_uu(1,2)*ddgb_J_tx
     &                            +2*g0_uu(1,3)*ddgb_J_ty
     &                            +2*g0_uu(1,4)*ddgb_J_tz
     &                                )
     &
     &                         -0.5d0*(
     &                            -dgb_J*
     &                            (g0_uu(1,1)*g0_uu(1,4)*g0_ll_x(4,1,1)+
     &                             g0_uu(1,1)*g0_uu(2,4)*g0_ll_x(4,1,2)+
     &                             g0_uu(1,1)*g0_uu(3,4)*g0_ll_x(4,1,3)+
     &                             g0_uu(1,1)*g0_uu(4,4)*g0_ll_x(4,1,4)+
     &                             g0_uu(2,1)*g0_uu(1,4)*g0_ll_x(4,2,1)+
     &                             g0_uu(2,1)*g0_uu(2,4)*g0_ll_x(4,2,2)+
     &                             g0_uu(2,1)*g0_uu(3,4)*g0_ll_x(4,2,3)+
     &                             g0_uu(2,1)*g0_uu(4,4)*g0_ll_x(4,2,4)+
     &                             g0_uu(3,1)*g0_uu(1,4)*g0_ll_x(4,3,1)+
     &                             g0_uu(3,1)*g0_uu(2,4)*g0_ll_x(4,3,2)+
     &                             g0_uu(3,1)*g0_uu(3,4)*g0_ll_x(4,3,3)+
     &                             g0_uu(3,1)*g0_uu(4,4)*g0_ll_x(4,3,4)+
     &                             g0_uu(4,1)*g0_uu(1,4)*g0_ll_x(4,4,1)+
     &                             g0_uu(4,1)*g0_uu(2,4)*g0_ll_x(4,4,2)+
     &                             g0_uu(4,1)*g0_uu(3,4)*g0_ll_x(4,4,3)+
     &                             g0_uu(4,1)*g0_uu(4,4)*g0_ll_x(4,4,4)+
     &
     &                             g0_uu(1,4)*g0_uu(1,1)*g0_ll_x(4,1,1)+
     &                             g0_uu(1,4)*g0_uu(2,1)*g0_ll_x(4,1,2)+
     &                             g0_uu(1,4)*g0_uu(3,1)*g0_ll_x(4,1,3)+
     &                             g0_uu(1,4)*g0_uu(4,1)*g0_ll_x(4,1,4)+
     &                             g0_uu(2,4)*g0_uu(1,1)*g0_ll_x(4,2,1)+
     &                             g0_uu(2,4)*g0_uu(2,1)*g0_ll_x(4,2,2)+
     &                             g0_uu(2,4)*g0_uu(3,1)*g0_ll_x(4,2,3)+
     &                             g0_uu(2,4)*g0_uu(4,1)*g0_ll_x(4,2,4)+
     &                             g0_uu(3,4)*g0_uu(1,1)*g0_ll_x(4,3,1)+
     &                             g0_uu(3,4)*g0_uu(2,1)*g0_ll_x(4,3,2)+
     &                             g0_uu(3,4)*g0_uu(3,1)*g0_ll_x(4,3,3)+
     &                             g0_uu(3,4)*g0_uu(4,1)*g0_ll_x(4,3,4)+
     &                             g0_uu(4,4)*g0_uu(1,1)*g0_ll_x(4,4,1)+
     &                             g0_uu(4,4)*g0_uu(2,1)*g0_ll_x(4,4,2)+
     &                             g0_uu(4,4)*g0_uu(3,1)*g0_ll_x(4,4,3)+
     &                             g0_uu(4,4)*g0_uu(4,1)*g0_ll_x(4,4,4))
     &                            +dgb_J*
     &                            (g0_uu_x(1,1,1))
     &                                )
     &
     &                         -0.5d0*(
     &                            dgb_J*
     &                            (g0_uu_x(4,1,4))
     &                                )
     &
     &                         -      (
     &                            0.5d0*dgb_J*
     &                            (cuuuu(1,1,1,4)*dlll(1,4,1)+
     &                             cuuuu(1,1,2,4)*dlll(2,4,1)+
     &                             cuuuu(1,1,3,4)*dlll(3,4,1)+
     &                             cuuuu(1,1,4,4)*dlll(4,4,1)+
     &                             cuuuu(1,2,1,4)*dlll(1,4,2)+
     &                             cuuuu(1,2,2,4)*dlll(2,4,2)+
     &                             cuuuu(1,2,3,4)*dlll(3,4,2)+
     &                             cuuuu(1,2,4,4)*dlll(4,4,2)+
     &                             cuuuu(1,3,1,4)*dlll(1,4,3)+
     &                             cuuuu(1,3,2,4)*dlll(2,4,3)+
     &                             cuuuu(1,3,3,4)*dlll(3,4,3)+
     &                             cuuuu(1,3,4,4)*dlll(4,4,3)+
     &                             cuuuu(1,4,1,4)*dlll(1,4,4)+
     &                             cuuuu(1,4,2,4)*dlll(2,4,4)+
     &                             cuuuu(1,4,3,4)*dlll(3,4,4)+
     &                             cuuuu(1,4,4,4)*dlll(4,4,4))
     &                                )


                do a=2,3
                  do b=a,3
                    efe_J(a,b)=-0.5d0*(
     &                            g0_uu(1,1)*ddgb_J
     &                            +2*g0_uu(1,2)*ddgb_J_tx
     &                            +2*g0_uu(1,3)*ddgb_J_ty
     &                            +2*g0_uu(1,4)*ddgb_J_tz
     &                                )
     &                      
     &                         -0.5d0*(
     &                            +dgb_J*
     &                            (g0_uu_x(a,1,a))
     &                                )
     &                      
     &                         -0.5d0*(
     &                            +dgb_J*
     &                            (g0_uu_x(b,1,b))
     &                                )
     &
     &                         +      (
     &                            -0.5d0*dgb_J*
     &                            ((Hads_l(1)+A_l(1))*g0_uu(1,1)+
     &                             (Hads_l(2)+A_l(2))*g0_uu(2,1)+
     &                             (Hads_l(3)+A_l(3))*g0_uu(3,1)+
     &                             (Hads_l(4)+A_l(4))*g0_uu(4,1))
     &                                )  
     &                      
     &                         -      (
     &                            0.25d0*dgb_J*
     &                            (cuuuu(1,a,1,1)*dlll(1,a,1)+
     &                             cuuuu(1,a,1,2)*dlll(1,a,2)+
     &                             cuuuu(1,a,1,3)*dlll(1,a,3)+
     &                             cuuuu(1,a,1,4)*dlll(1,a,4)+
     &                             cuuuu(2,a,1,1)*dlll(2,a,1)+
     &                             cuuuu(2,a,1,2)*dlll(2,a,2)+
     &                             cuuuu(2,a,1,3)*dlll(2,a,3)+
     &                             cuuuu(2,a,1,4)*dlll(2,a,4)+
     &                             cuuuu(3,a,1,1)*dlll(3,a,1)+
     &                             cuuuu(3,a,1,2)*dlll(3,a,2)+
     &                             cuuuu(3,a,1,3)*dlll(3,a,3)+
     &                             cuuuu(3,a,1,4)*dlll(3,a,4)+
     &                             cuuuu(4,a,1,1)*dlll(4,a,1)+
     &                             cuuuu(4,a,1,2)*dlll(4,a,2)+
     &                             cuuuu(4,a,1,3)*dlll(4,a,3)+
     &                             cuuuu(4,a,1,4)*dlll(4,a,4)-
     &                             cuuuu(1,1,a,1)*dlll(1,a,1)-
     &                             cuuuu(1,1,a,2)*dlll(1,a,2)-
     &                             cuuuu(1,1,a,3)*dlll(1,a,3)-
     &                             cuuuu(1,1,a,4)*dlll(1,a,4)-
     &                             cuuuu(2,1,a,1)*dlll(2,a,1)-
     &                             cuuuu(2,1,a,2)*dlll(2,a,2)-
     &                             cuuuu(2,1,a,3)*dlll(2,a,3)-
     &                             cuuuu(2,1,a,4)*dlll(2,a,4)-
     &                             cuuuu(3,1,a,1)*dlll(3,a,1)-
     &                             cuuuu(3,1,a,2)*dlll(3,a,2)-
     &                             cuuuu(3,1,a,3)*dlll(3,a,3)-
     &                             cuuuu(3,1,a,4)*dlll(3,a,4)-
     &                             cuuuu(4,1,a,1)*dlll(4,a,1)-
     &                             cuuuu(4,1,a,2)*dlll(4,a,2)-
     &                             cuuuu(4,1,a,3)*dlll(4,a,3)-
     &                             cuuuu(4,1,a,4)*dlll(4,a,4)
     &                           +
     &                             cuuuu(1,1,1,b)*dlll(1,b,1)+
     &                             cuuuu(1,1,2,b)*dlll(2,b,1)+
     &                             cuuuu(1,1,3,b)*dlll(3,b,1)+
     &                             cuuuu(1,1,4,b)*dlll(4,b,1)+
     &                             cuuuu(1,2,1,b)*dlll(1,b,2)+
     &                             cuuuu(1,2,2,b)*dlll(2,b,2)+
     &                             cuuuu(1,2,3,b)*dlll(3,b,2)+
     &                             cuuuu(1,2,4,b)*dlll(4,b,2)+
     &                             cuuuu(1,3,1,b)*dlll(1,b,3)+
     &                             cuuuu(1,3,2,b)*dlll(2,b,3)+
     &                             cuuuu(1,3,3,b)*dlll(3,b,3)+
     &                             cuuuu(1,3,4,b)*dlll(4,b,3)+
     &                             cuuuu(1,4,1,b)*dlll(1,b,4)+
     &                             cuuuu(1,4,2,b)*dlll(2,b,4)+
     &                             cuuuu(1,4,3,b)*dlll(3,b,4)+
     &                             cuuuu(1,4,4,b)*dlll(4,b,4)-
     &                             cuuuu(b,1,1,1)*dlll(1,b,1)-
     &                             cuuuu(b,1,2,1)*dlll(2,b,1)-
     &                             cuuuu(b,1,3,1)*dlll(3,b,1)-
     &                             cuuuu(b,1,4,1)*dlll(4,b,1)-
     &                             cuuuu(b,2,1,1)*dlll(1,b,2)-
     &                             cuuuu(b,2,2,1)*dlll(2,b,2)-
     &                             cuuuu(b,2,3,1)*dlll(3,b,2)-
     &                             cuuuu(b,2,4,1)*dlll(4,b,2)-
     &                             cuuuu(b,3,1,1)*dlll(1,b,3)-
     &                             cuuuu(b,3,2,1)*dlll(2,b,3)-
     &                             cuuuu(b,3,3,1)*dlll(3,b,3)-
     &                             cuuuu(b,3,4,1)*dlll(4,b,3)-
     &                             cuuuu(b,4,1,1)*dlll(1,b,4)-
     &                             cuuuu(b,4,2,1)*dlll(2,b,4)-
     &                             cuuuu(b,4,3,1)*dlll(3,b,4)-
     &                             cuuuu(b,4,4,1)*dlll(4,b,4))
     &                                )
                  end do
                end do

                efe_J(2,4)=-0.5d0*(
     &                            g0_uu(1,1)*ddgb_J
     &                            +2*g0_uu(1,2)*ddgb_J_tx
     &                            +2*g0_uu(1,3)*ddgb_J_ty
     &                            +2*g0_uu(1,4)*ddgb_J_tz
     &                                )
     &
     &                         -0.5d0*(
     &                            +dgb_J*
     &                            (g0_uu_x(2,1,2))
     &                                )
     &
     &                         -0.5d0*(
     &                            +dgb_J*
     &                            (g0_uu_x(4,1,4))
     &                                )
     &
     &                         +      (
     &                            -0.5d0*dgb_J*
     &                            ((Hads_l(1)+A_l(1))*g0_uu(1,1)+
     &                             (Hads_l(2)+A_l(2))*g0_uu(2,1)+
     &                             (Hads_l(3)+A_l(3))*g0_uu(3,1)+
     &                             (Hads_l(4)+A_l(4))*g0_uu(4,1))
     &                                )
     &
     &                         -      (
     &                            0.25d0*dgb_J*
     &                            (cuuuu(1,2,1,1)*dlll(1,2,1)+
     &                             cuuuu(1,2,1,2)*dlll(1,2,2)+
     &                             cuuuu(1,2,1,3)*dlll(1,2,3)+
     &                             cuuuu(1,2,1,4)*dlll(1,2,4)+
     &                             cuuuu(2,2,1,1)*dlll(2,2,1)+
     &                             cuuuu(2,2,1,2)*dlll(2,2,2)+
     &                             cuuuu(2,2,1,3)*dlll(2,2,3)+
     &                             cuuuu(2,2,1,4)*dlll(2,2,4)+
     &                             cuuuu(3,2,1,1)*dlll(3,2,1)+
     &                             cuuuu(3,2,1,2)*dlll(3,2,2)+
     &                             cuuuu(3,2,1,3)*dlll(3,2,3)+
     &                             cuuuu(3,2,1,4)*dlll(3,2,4)+
     &                             cuuuu(4,2,1,1)*dlll(4,2,1)+
     &                             cuuuu(4,2,1,2)*dlll(4,2,2)+
     &                             cuuuu(4,2,1,3)*dlll(4,2,3)+
     &                             cuuuu(4,2,1,4)*dlll(4,2,4)-
     &                             cuuuu(1,1,2,1)*dlll(1,2,1)-
     &                             cuuuu(1,1,2,2)*dlll(1,2,2)-
     &                             cuuuu(1,1,2,3)*dlll(1,2,3)-
     &                             cuuuu(1,1,2,4)*dlll(1,2,4)-
     &                             cuuuu(2,1,2,1)*dlll(2,2,1)-
     &                             cuuuu(2,1,2,2)*dlll(2,2,2)-
     &                             cuuuu(2,1,2,3)*dlll(2,2,3)-
     &                             cuuuu(2,1,2,4)*dlll(2,2,4)-
     &                             cuuuu(3,1,2,1)*dlll(3,2,1)-
     &                             cuuuu(3,1,2,2)*dlll(3,2,2)-
     &                             cuuuu(3,1,2,3)*dlll(3,2,3)-
     &                             cuuuu(3,1,2,4)*dlll(3,2,4)-
     &                             cuuuu(4,1,2,1)*dlll(4,2,1)-
     &                             cuuuu(4,1,2,2)*dlll(4,2,2)-
     &                             cuuuu(4,1,2,3)*dlll(4,2,3)-
     &                             cuuuu(4,1,2,4)*dlll(4,2,4)
     &                           +
     &                             cuuuu(1,1,1,4)*dlll(1,4,1)+
     &                             cuuuu(1,1,2,4)*dlll(2,4,1)+
     &                             cuuuu(1,1,3,4)*dlll(3,4,1)+
     &                             cuuuu(1,1,4,4)*dlll(4,4,1)+
     &                             cuuuu(1,2,1,4)*dlll(1,4,2)+
     &                             cuuuu(1,2,2,4)*dlll(2,4,2)+
     &                             cuuuu(1,2,3,4)*dlll(3,4,2)+
     &                             cuuuu(1,2,4,4)*dlll(4,4,2)+
     &                             cuuuu(1,3,1,4)*dlll(1,4,3)+
     &                             cuuuu(1,3,2,4)*dlll(2,4,3)+
     &                             cuuuu(1,3,3,4)*dlll(3,4,3)+
     &                             cuuuu(1,3,4,4)*dlll(4,4,3)+
     &                             cuuuu(1,4,1,4)*dlll(1,4,4)+
     &                             cuuuu(1,4,2,4)*dlll(2,4,4)+
     &                             cuuuu(1,4,3,4)*dlll(3,4,4)+
     &                             cuuuu(1,4,4,4)*dlll(4,4,4)-
     &                             cuuuu(4,1,1,1)*dlll(1,4,1)-
     &                             cuuuu(4,1,2,1)*dlll(2,4,1)-
     &                             cuuuu(4,1,3,1)*dlll(3,4,1)-
     &                             cuuuu(4,1,4,1)*dlll(4,4,1)-
     &                             cuuuu(4,2,1,1)*dlll(1,4,2)-
     &                             cuuuu(4,2,2,1)*dlll(2,4,2)-
     &                             cuuuu(4,2,3,1)*dlll(3,4,2)-
     &                             cuuuu(4,2,4,1)*dlll(4,4,2)-
     &                             cuuuu(4,3,1,1)*dlll(1,4,3)-
     &                             cuuuu(4,3,2,1)*dlll(2,4,3)-
     &                             cuuuu(4,3,3,1)*dlll(3,4,3)-
     &                             cuuuu(4,3,4,1)*dlll(4,4,3)-
     &                             cuuuu(4,4,1,1)*dlll(1,4,4)-
     &                             cuuuu(4,4,2,1)*dlll(2,4,4)-
     &                             cuuuu(4,4,3,1)*dlll(3,4,4)-
     &                             cuuuu(4,4,4,1)*dlll(4,4,4))
     &                                )

                efe_J(3,4)=-0.5d0*(
     &                            g0_uu(1,1)*ddgb_J
     &                            +2*g0_uu(1,2)*ddgb_J_tx
     &                            +2*g0_uu(1,3)*ddgb_J_ty
     &                            +2*g0_uu(1,4)*ddgb_J_tz
     &                                )
     &
     &                         -0.5d0*(
     &                            +dgb_J*
     &                            (g0_uu_x(3,1,3))
     &                                )
     &
     &                         -0.5d0*(
     &                            +dgb_J*
     &                            (g0_uu_x(4,1,4))
     &                                )
     &
     &                         +      (
     &                            -0.5d0*dgb_J*
     &                            ((Hads_l(1)+A_l(1))*g0_uu(1,1)+
     &                             (Hads_l(2)+A_l(2))*g0_uu(2,1)+
     &                             (Hads_l(3)+A_l(3))*g0_uu(3,1)+
     &                             (Hads_l(4)+A_l(4))*g0_uu(4,1))
     &                                )
     &
     &                         -      (
     &                            0.25d0*dgb_J*
     &                            (cuuuu(1,3,1,1)*dlll(1,3,1)+
     &                             cuuuu(1,3,1,2)*dlll(1,3,2)+
     &                             cuuuu(1,3,1,3)*dlll(1,3,3)+
     &                             cuuuu(1,3,1,4)*dlll(1,3,4)+
     &                             cuuuu(2,3,1,1)*dlll(2,3,1)+
     &                             cuuuu(2,3,1,2)*dlll(2,3,2)+
     &                             cuuuu(2,3,1,3)*dlll(2,3,3)+
     &                             cuuuu(2,3,1,4)*dlll(2,3,4)+
     &                             cuuuu(3,3,1,1)*dlll(3,3,1)+
     &                             cuuuu(3,3,1,2)*dlll(3,3,2)+
     &                             cuuuu(3,3,1,3)*dlll(3,3,3)+
     &                             cuuuu(3,3,1,4)*dlll(3,3,4)+
     &                             cuuuu(4,3,1,1)*dlll(4,3,1)+
     &                             cuuuu(4,3,1,2)*dlll(4,3,2)+
     &                             cuuuu(4,3,1,3)*dlll(4,3,3)+
     &                             cuuuu(4,3,1,4)*dlll(4,3,4)-
     &                             cuuuu(1,1,3,1)*dlll(1,3,1)-
     &                             cuuuu(1,1,3,2)*dlll(1,3,2)-
     &                             cuuuu(1,1,3,3)*dlll(1,3,3)-
     &                             cuuuu(1,1,3,4)*dlll(1,3,4)-
     &                             cuuuu(2,1,3,1)*dlll(2,3,1)-
     &                             cuuuu(2,1,3,2)*dlll(2,3,2)-
     &                             cuuuu(2,1,3,3)*dlll(2,3,3)-
     &                             cuuuu(2,1,3,4)*dlll(2,3,4)-
     &                             cuuuu(3,1,3,1)*dlll(3,3,1)-
     &                             cuuuu(3,1,3,2)*dlll(3,3,2)-
     &                             cuuuu(3,1,3,3)*dlll(3,3,3)-
     &                             cuuuu(3,1,3,4)*dlll(3,3,4)-
     &                             cuuuu(4,1,3,1)*dlll(4,3,1)-
     &                             cuuuu(4,1,3,2)*dlll(4,3,2)-
     &                             cuuuu(4,1,3,3)*dlll(4,3,3)-
     &                             cuuuu(4,1,3,4)*dlll(4,3,4)
     &                           +
     &                             cuuuu(1,1,1,4)*dlll(1,4,1)+
     &                             cuuuu(1,1,2,4)*dlll(2,4,1)+
     &                             cuuuu(1,1,3,4)*dlll(3,4,1)+
     &                             cuuuu(1,1,4,4)*dlll(4,4,1)+
     &                             cuuuu(1,2,1,4)*dlll(1,4,2)+
     &                             cuuuu(1,2,2,4)*dlll(2,4,2)+
     &                             cuuuu(1,2,3,4)*dlll(3,4,2)+
     &                             cuuuu(1,2,4,4)*dlll(4,4,2)+
     &                             cuuuu(1,3,1,4)*dlll(1,4,3)+
     &                             cuuuu(1,3,2,4)*dlll(2,4,3)+
     &                             cuuuu(1,3,3,4)*dlll(3,4,3)+
     &                             cuuuu(1,3,4,4)*dlll(4,4,3)+
     &                             cuuuu(1,4,1,4)*dlll(1,4,4)+
     &                             cuuuu(1,4,2,4)*dlll(2,4,4)+
     &                             cuuuu(1,4,3,4)*dlll(3,4,4)+
     &                             cuuuu(1,4,4,4)*dlll(4,4,4)-
     &                             cuuuu(4,1,1,1)*dlll(1,4,1)-
     &                             cuuuu(4,1,2,1)*dlll(2,4,1)-
     &                             cuuuu(4,1,3,1)*dlll(3,4,1)-
     &                             cuuuu(4,1,4,1)*dlll(4,4,1)-
     &                             cuuuu(4,2,1,1)*dlll(1,4,2)-
     &                             cuuuu(4,2,2,1)*dlll(2,4,2)-
     &                             cuuuu(4,2,3,1)*dlll(3,4,2)-
     &                             cuuuu(4,2,4,1)*dlll(4,4,2)-
     &                             cuuuu(4,3,1,1)*dlll(1,4,3)-
     &                             cuuuu(4,3,2,1)*dlll(2,4,3)-
     &                             cuuuu(4,3,3,1)*dlll(3,4,3)-
     &                             cuuuu(4,3,4,1)*dlll(4,4,3)-
     &                             cuuuu(4,4,1,1)*dlll(1,4,4)-
     &                             cuuuu(4,4,2,1)*dlll(2,4,4)-
     &                             cuuuu(4,4,3,1)*dlll(3,4,4)-
     &                             cuuuu(4,4,4,1)*dlll(4,4,4))
     &                                )

!!!!!!!!!!2+1 version!!!!!
!                efe_J(4,4)=    -0.5d0*(
!     &                            y0**2*g0_uu(1,1)*ddgb_J
!     &                            +4*y0*g0_uu(1,3)
!     &                                                   *dgb_J
!     &                            +2*y0**2*g0_uu(1,2)
!     &                                                *ddgb_J_tx
!     &                            +2*y0**2*g0_uu(1,3)
!     &                                                *ddgb_J_ty
!     &                                )
!     &                      
!     &                         -0.5d0*(
!     &                            +y0**2*dgb_J*
!     &                            (g0_uu_x(4,1,4))
!     &                                )
!     &                      
!     &                         -0.5d0*(
!     &                            +y0**2*dgb_J*
!     &                            (g0_uu_x(4,1,4))
!     &                                )
!     &
!     &                         +      (
!     &                            -0.5d0*y0**2*dgb_J*
!     &                            ((Hads_l(1)+A_l(1))*g0_uu(1,1)+
!     &                             (Hads_l(2)+A_l(2))*g0_uu(2,1)+
!     &                             (Hads_l(3)+A_l(3))*g0_uu(3,1)+
!     &                             (Hads_l(4)+A_l(4))*g0_uu(4,1))
!     &                                )  
!     &                      
!     &                         -      (
!     &                            0.25d0*y0**2*dgb_J*
!     &                            (cuuuu(1,4,1,1)*dlll(1,4,1)+
!     &                             cuuuu(1,4,1,2)*dlll(1,4,2)+
!     &                             cuuuu(1,4,1,3)*dlll(1,4,3)+
!     &                             cuuuu(1,4,1,4)*dlll(1,4,4)+
!     &                             cuuuu(2,4,1,1)*dlll(2,4,1)+
!     &                             cuuuu(2,4,1,2)*dlll(2,4,2)+
!     &                             cuuuu(2,4,1,3)*dlll(2,4,3)+
!     &                             cuuuu(2,4,1,4)*dlll(2,4,4)+
!     &                             cuuuu(3,4,1,1)*dlll(3,4,1)+
!     &                             cuuuu(3,4,1,2)*dlll(3,4,2)+
!     &                             cuuuu(3,4,1,3)*dlll(3,4,3)+
!     &                             cuuuu(3,4,1,4)*dlll(3,4,4)+
!     &                             cuuuu(4,4,1,1)*dlll(4,4,1)+
!     &                             cuuuu(4,4,1,2)*dlll(4,4,2)+
!     &                             cuuuu(4,4,1,3)*dlll(4,4,3)+
!     &                             cuuuu(4,4,1,4)*dlll(4,4,4)-
!     &                             cuuuu(1,1,4,1)*dlll(1,4,1)-
!     &                             cuuuu(1,1,4,2)*dlll(1,4,2)-
!     &                             cuuuu(1,1,4,3)*dlll(1,4,3)-
!     &                             cuuuu(1,1,4,4)*dlll(1,4,4)-
!     &                             cuuuu(2,1,4,1)*dlll(2,4,1)-
!     &                             cuuuu(2,1,4,2)*dlll(2,4,2)-
!     &                             cuuuu(2,1,4,3)*dlll(2,4,3)-
!     &                             cuuuu(2,1,4,4)*dlll(2,4,4)-
!     &                             cuuuu(3,1,4,1)*dlll(3,4,1)-
!     &                             cuuuu(3,1,4,2)*dlll(3,4,2)-
!     &                             cuuuu(3,1,4,3)*dlll(3,4,3)-
!     &                             cuuuu(3,1,4,4)*dlll(3,4,4)-
!     &                             cuuuu(4,1,4,1)*dlll(4,4,1)-
!     &                             cuuuu(4,1,4,2)*dlll(4,4,2)-
!     &                             cuuuu(4,1,4,3)*dlll(4,4,3)-
!     &                             cuuuu(4,1,4,4)*dlll(4,4,4)
!     &                           +
!     &                             cuuuu(1,1,1,4)*dlll(1,4,1)+
!     &                             cuuuu(1,1,2,4)*dlll(2,4,1)+
!     &                             cuuuu(1,1,3,4)*dlll(3,4,1)+
!     &                             cuuuu(1,1,4,4)*dlll(4,4,1)+
!     &                             cuuuu(1,2,1,4)*dlll(1,4,2)+
!     &                             cuuuu(1,2,2,4)*dlll(2,4,2)+
!     &                             cuuuu(1,2,3,4)*dlll(3,4,2)+
!     &                             cuuuu(1,2,4,4)*dlll(4,4,2)+
!     &                             cuuuu(1,3,1,4)*dlll(1,4,3)+
!     &                             cuuuu(1,3,2,4)*dlll(2,4,3)+
!     &                             cuuuu(1,3,3,4)*dlll(3,4,3)+
!     &                             cuuuu(1,3,4,4)*dlll(4,4,3)+
!     &                             cuuuu(1,4,1,4)*dlll(1,4,4)+
!     &                             cuuuu(1,4,2,4)*dlll(2,4,4)+
!     &                             cuuuu(1,4,3,4)*dlll(3,4,4)+
!     &                             cuuuu(1,4,4,4)*dlll(4,4,4)-
!     &                             cuuuu(4,1,1,1)*dlll(1,4,1)-
!     &                             cuuuu(4,1,2,1)*dlll(2,4,1)-
!     &                             cuuuu(4,1,3,1)*dlll(3,4,1)-
!     &                             cuuuu(4,1,4,1)*dlll(4,4,1)-
!     &                             cuuuu(4,2,1,1)*dlll(1,4,2)-
!     &                             cuuuu(4,2,2,1)*dlll(2,4,2)-
!     &                             cuuuu(4,2,3,1)*dlll(3,4,2)-
!     &                             cuuuu(4,2,4,1)*dlll(4,4,2)-
!     &                             cuuuu(4,3,1,1)*dlll(1,4,3)-
!     &                             cuuuu(4,3,2,1)*dlll(2,4,3)-
!     &                             cuuuu(4,3,3,1)*dlll(3,4,3)-
!     &                             cuuuu(4,3,4,1)*dlll(4,4,3)-
!     &                             cuuuu(4,4,1,1)*dlll(1,4,4)-
!     &                             cuuuu(4,4,2,1)*dlll(2,4,4)-
!     &                             cuuuu(4,4,3,1)*dlll(3,4,4)-
!     &                             cuuuu(4,4,4,1)*dlll(4,4,4))
!     &                                )
!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    efe_J(4,4)=-0.5d0*(
     &                            g0_uu(1,1)*ddgb_J
     &                            +2*g0_uu(1,2)*ddgb_J_tx
     &                            +2*g0_uu(1,3)*ddgb_J_ty
     &                            +2*g0_uu(1,4)*ddgb_J_tz
     &                                )
     &
     &                         -0.5d0*(
     &                            +dgb_J*
     &                            (g0_uu_x(4,1,4))
     &                                )
     &
     &                         -0.5d0*(
     &                            +dgb_J*
     &                            (g0_uu_x(4,1,4))
     &                                )
     &
     &                         +      (
     &                            -0.5d0*dgb_J*
     &                            ((Hads_l(1)+A_l(1))*g0_uu(1,1)+
     &                             (Hads_l(2)+A_l(2))*g0_uu(2,1)+
     &                             (Hads_l(3)+A_l(3))*g0_uu(3,1)+
     &                             (Hads_l(4)+A_l(4))*g0_uu(4,1))
     &                                )
     &
     &                         -      (
     &                            0.25d0*dgb_J*
     &                            (cuuuu(1,4,1,1)*dlll(1,4,1)+
     &                             cuuuu(1,4,1,2)*dlll(1,4,2)+
     &                             cuuuu(1,4,1,3)*dlll(1,4,3)+
     &                             cuuuu(1,4,1,4)*dlll(1,4,4)+
     &                             cuuuu(2,4,1,1)*dlll(2,4,1)+
     &                             cuuuu(2,4,1,2)*dlll(2,4,2)+
     &                             cuuuu(2,4,1,3)*dlll(2,4,3)+
     &                             cuuuu(2,4,1,4)*dlll(2,4,4)+
     &                             cuuuu(3,4,1,1)*dlll(3,4,1)+
     &                             cuuuu(3,4,1,2)*dlll(3,4,2)+
     &                             cuuuu(3,4,1,3)*dlll(3,4,3)+
     &                             cuuuu(3,4,1,4)*dlll(3,4,4)+
     &                             cuuuu(4,4,1,1)*dlll(4,4,1)+
     &                             cuuuu(4,4,1,2)*dlll(4,4,2)+
     &                             cuuuu(4,4,1,3)*dlll(4,4,3)+
     &                             cuuuu(4,4,1,4)*dlll(4,4,4)-
     &                             cuuuu(1,1,4,1)*dlll(1,4,1)-
     &                             cuuuu(1,1,4,2)*dlll(1,4,2)-
     &                             cuuuu(1,1,4,3)*dlll(1,4,3)-
     &                             cuuuu(1,1,4,4)*dlll(1,4,4)-
     &                             cuuuu(2,1,4,1)*dlll(2,4,1)-
     &                             cuuuu(2,1,4,2)*dlll(2,4,2)-
     &                             cuuuu(2,1,4,3)*dlll(2,4,3)-
     &                             cuuuu(2,1,4,4)*dlll(2,4,4)-
     &                             cuuuu(3,1,4,1)*dlll(3,4,1)-
     &                             cuuuu(3,1,4,2)*dlll(3,4,2)-
     &                             cuuuu(3,1,4,3)*dlll(3,4,3)-
     &                             cuuuu(3,1,4,4)*dlll(3,4,4)-
     &                             cuuuu(4,1,4,1)*dlll(4,4,1)-
     &                             cuuuu(4,1,4,2)*dlll(4,4,2)-
     &                             cuuuu(4,1,4,3)*dlll(4,4,3)-
     &                             cuuuu(4,1,4,4)*dlll(4,4,4)
     &                           +
     &                             cuuuu(1,1,1,4)*dlll(1,4,1)+
     &                             cuuuu(1,1,2,4)*dlll(2,4,1)+
     &                             cuuuu(1,1,3,4)*dlll(3,4,1)+
     &                             cuuuu(1,1,4,4)*dlll(4,4,1)+
     &                             cuuuu(1,2,1,4)*dlll(1,4,2)+
     &                             cuuuu(1,2,2,4)*dlll(2,4,2)+
     &                             cuuuu(1,2,3,4)*dlll(3,4,2)+
     &                             cuuuu(1,2,4,4)*dlll(4,4,2)+
     &                             cuuuu(1,3,1,4)*dlll(1,4,3)+
     &                             cuuuu(1,3,2,4)*dlll(2,4,3)+
     &                             cuuuu(1,3,3,4)*dlll(3,4,3)+
     &                             cuuuu(1,3,4,4)*dlll(4,4,3)+
     &                             cuuuu(1,4,1,4)*dlll(1,4,4)+
     &                             cuuuu(1,4,2,4)*dlll(2,4,4)+
     &                             cuuuu(1,4,3,4)*dlll(3,4,4)+
     &                             cuuuu(1,4,4,4)*dlll(4,4,4)-
     &                             cuuuu(4,1,1,1)*dlll(1,4,1)-
     &                             cuuuu(4,1,2,1)*dlll(2,4,1)-
     &                             cuuuu(4,1,3,1)*dlll(3,4,1)-
     &                             cuuuu(4,1,4,1)*dlll(4,4,1)-
     &                             cuuuu(4,2,1,1)*dlll(1,4,2)-
     &                             cuuuu(4,2,2,1)*dlll(2,4,2)-
     &                             cuuuu(4,2,3,1)*dlll(3,4,2)-
     &                             cuuuu(4,2,4,1)*dlll(4,4,2)-
     &                             cuuuu(4,3,1,1)*dlll(1,4,3)-
     &                             cuuuu(4,3,2,1)*dlll(2,4,3)-
     &                             cuuuu(4,3,3,1)*dlll(3,4,3)-
     &                             cuuuu(4,3,4,1)*dlll(4,4,3)-
     &                             cuuuu(4,4,1,1)*dlll(1,4,4)-
     &                             cuuuu(4,4,2,1)*dlll(2,4,4)-
     &                             cuuuu(4,4,3,1)*dlll(3,4,4)-
     &                             cuuuu(4,4,4,1)*dlll(4,4,4))
     &                                )

                !----------------------------------------------------------------
                ! computes diag. Jacobian of phi1_np1->L.phi1_np1 transformation
                ! by differentiating L.phi1=box.phi1-dV/dphi1 wrt. phi1_np1
                ! and remember: phi10=phi1*(1-rho0**2)**2
                ! 
                ! (re-use the dgb_J,ddgb_J,ddgb_J_tx,ddgb_J_ty defined above)
                !----------------------------------------------------------------
                phi1_J=    (
     &                 g0_uu(1,1)*ddgb_J
     &                           *(1-rho0**2)**2
     &                 -4*(3)*x0*g0_uu(1,2)*dgb_J
     &                           *(1-rho0**2)**2
     &                 +2*g0_uu(1,2)*ddgb_J_tx
     &                           *(1-rho0**2)**2
     &                 -4*(3)*y0*g0_uu(1,3)*dgb_J
     &                           *(1-rho0**2)**2
     &                 +2*g0_uu(1,3)*ddgb_J_ty
     &                           *(1-rho0**2)**2
     &                 +2*g0_uu(1,4)*ddgb_J_tz
     &                           *(1-rho0**2)**2
     &                     )
     &                +
     &                  dgb_J*(1-rho0**2)**2
     &                  *(
     &                    g0_uu_x(1,1,1)+
     &                    g0_uu_x(2,1,2)+
     &                    g0_uu_x(3,1,3)+
     &                    g0_uu_x(4,1,4)
     &                   )
     &                +
     &                  dgb_J*(1-rho0**2)**2
     &                  *(
     &                    g0_uu(1,1)*gamma_ull(1,1,1)+
     &                    g0_uu(1,1)*gamma_ull(2,2,1)+
     &                    g0_uu(1,1)*gamma_ull(3,3,1)+
     &                    g0_uu(1,1)*gamma_ull(4,4,1)+
     &                    g0_uu(1,2)*gamma_ull(1,1,2)+
     &                    g0_uu(1,2)*gamma_ull(2,2,2)+
     &                    g0_uu(1,2)*gamma_ull(3,3,2)+
     &                    g0_uu(1,2)*gamma_ull(4,4,2)+
     &                    g0_uu(1,3)*gamma_ull(1,1,3)+
     &                    g0_uu(1,3)*gamma_ull(2,2,3)+
     &                    g0_uu(1,3)*gamma_ull(3,3,3)+
     &                    g0_uu(1,3)*gamma_ull(4,4,3)+
     &                    g0_uu(1,4)*gamma_ull(1,1,4)+
     &                    g0_uu(1,4)*gamma_ull(2,2,4)+
     &                    g0_uu(1,4)*gamma_ull(3,3,4)+
     &                    g0_uu(1,4)*gamma_ull(4,4,4)
     &                   )

                ! constraint damping terms added to efe,efe_J
                do a=1,4
                  do b=1,4
                    cd_ll(a,b)=-kappa_cd*
     &                 (
     &                  n_l(a)*c_l(b)+n_l(b)*c_l(a)
     &                 -(1+rho_cd)*g0_ll(a,b)*ndotc
     &                 )
                  end do
                end do

                dc_J=1.0d0/2.0d0/dt

                cd_J_ll(1,1)=-kappa_cd*
     &              (
     &               n_l(1)*g0_uu(1,1)*dc_J
     &              -(1+rho_cd)*g0_ll(1,1)*n_u(1)*
     &                    (0.5d0*g0_uu(1,1)*dc_J)
     &              )
                cd_J_ll(1,2)=-kappa_cd*
     &              (
     &               n_l(1)*g0_uu(1,1)*dc_J
     &              -(1+rho_cd)*g0_ll(1,2)*n_u(2)*
     &                    (g0_uu(1,1)*dc_J)
     &              )
                cd_J_ll(1,3)=-kappa_cd*
     &              (
     &               n_l(1)*g0_uu(1,1)*dc_J
     &              -(1+rho_cd)*g0_ll(1,3)*n_u(3)*
     &                    (g0_uu(1,1)*dc_J)
     &              )
                cd_J_ll(1,4)=-kappa_cd*
     &              (
     &               n_l(1)*g0_uu(1,1)*dc_J
     &              -(1+rho_cd)*g0_ll(1,4)*n_u(4)*
     &                    (g0_uu(1,1)*dc_J)
     &              )

                cd_J_ll(2,2)=-kappa_cd*
     &              (
     &               2*n_l(2)*g0_uu(1,2)*dc_J
     &              -(1+rho_cd)*g0_ll(2,2)*
     &                    (-n_u(1)*(0.5d0*g0_uu(2,2)*dc_J)
     &                     +n_u(2)*(g0_uu(1,2)*dc_J))
     &              )
                cd_J_ll(2,3)=-kappa_cd*
     &              (
     &               n_l(2)*g0_uu(1,2)*dc_J
     &              +n_l(3)*g0_uu(1,3)*dc_J
     &              -(1+rho_cd)*g0_ll(2,2)*
     &                    (-n_u(1)*(g0_uu(2,3)*dc_J)
     &                     +n_u(2)*(g0_uu(1,3)*dc_J)
     &                     +n_u(3)*(g0_uu(1,2)*dc_J))
     &              )
                cd_J_ll(2,4)=-kappa_cd*
     &              (
     &               n_l(2)*g0_uu(1,2)*dc_J
     &              +n_l(4)*g0_uu(1,4)*dc_J
     &              -(1+rho_cd)*g0_ll(2,2)*
     &                    (-n_u(1)*(g0_uu(2,4)*dc_J)
     &                     +n_u(2)*(g0_uu(1,4)*dc_J)
     &                     +n_u(4)*(g0_uu(1,2)*dc_J))
     &              )
                cd_J_ll(3,3)=-kappa_cd*
     &              (
     &               2*n_l(3)*g0_uu(1,3)*dc_J
     &              -(1+rho_cd)*g0_ll(3,3)*
     &                    (-n_u(1)*(0.5d0*g0_uu(3,3)*dc_J)
     &                     +n_u(3)*(g0_uu(1,3)*dc_J))
     &              )
                cd_J_ll(3,4)=-kappa_cd*
     &              (
     &               n_l(3)*g0_uu(1,3)*dc_J
     &              +n_l(4)*g0_uu(1,4)*dc_J
     &              -(1+rho_cd)*g0_ll(3,3)*
     &                    (-n_u(1)*(g0_uu(3,4)*dc_J)
     &                     +n_u(3)*(g0_uu(1,4)*dc_J)
     &                     +n_u(4)*(g0_uu(1,3)*dc_J))
     &              )
!!!!!!2+1 version!!!!!!!
!                cd_J_ll(4,4)=-kappa_cd*
!     &              (
!     &              -(1+rho_cd)*g0_ll(4,4)*n_u(1)*
!     &                    (-0.5d0*g0_uu(4,4)*dc_J*y0**2)
!     &              )
!!!!!!!!!!!!!!!!!!!
                cd_J_ll(4,4)=-kappa_cd*
     &              (
     &               2*n_l(4)*g0_uu(1,4)*dc_J
     &              -(1+rho_cd)*g0_ll(4,4)*
     &                    (-n_u(1)*(0.5d0*g0_uu(4,4)*dc_J)
     &                     +n_u(4)*(g0_uu(1,4)*dc_J))
     &              )

                if (kappa_cd.ne.0) then
                  efe(1,1)=efe(1,1)+cd_ll(1,1)
                  efe(1,2)=efe(1,2)+cd_ll(1,2)
                  efe(1,3)=efe(1,3)+cd_ll(1,3)
                  efe(1,4)=efe(1,4)+cd_ll(1,4)
                  efe(2,2)=efe(2,2)+cd_ll(2,2)
                  efe(2,3)=efe(2,3)+cd_ll(2,3)
                  efe(2,4)=efe(2,4)+cd_ll(2,4)
                  efe(3,3)=efe(3,3)+cd_ll(3,3)
                  efe(3,4)=efe(3,4)+cd_ll(3,4)
                  efe(4,4)=efe(4,4)+cd_ll(4,4)
                  efe_J(1,1)=efe_J(1,1)+cd_J_ll(1,1)
                  efe_J(1,2)=efe_J(1,2)+cd_J_ll(1,2)
                  efe_J(1,3)=efe_J(1,3)+cd_J_ll(1,3)
                  efe_J(1,4)=efe_J(1,4)+cd_J_ll(1,4)
                  efe_J(2,2)=efe_J(2,2)+cd_J_ll(2,2)
                  efe_J(2,3)=efe_J(2,3)+cd_J_ll(2,3)
                  efe_J(2,4)=efe_J(2,4)+cd_J_ll(2,4)
                  efe_J(3,3)=efe_J(3,3)+cd_J_ll(3,3)
                  efe_J(3,4)=efe_J(3,4)+cd_J_ll(3,4)
                  efe_J(4,4)=efe_J(4,4)+cd_J_ll(4,4)
                end if

                ! update gbars 
                if (background.eq.0) then 
                  if (is_nan(efe(1,1)).or.is_nan(efe_J(1,1)).or.
     &              efe_J(1,1).eq.0) then
                    dump=.true.
                  else
                    gb_tt_np1(i,j,k)=gb_tt_np1(i,j,k)
     &                               -efe(1,1)/efe_J(1,1)
                  end if
 
                  if (is_nan(efe(1,2)).or.is_nan(efe_J(1,2)).or.
     &              efe_J(1,2).eq.0) then
                    dump=.true.
                  else
                    gb_tx_np1(i,j,k)=gb_tx_np1(i,j,k)
     &                               -efe(1,2)/efe_J(1,2)
                  end if

                  if (is_nan(efe(1,3)).or.is_nan(efe_J(1,3)).or.
     &              efe_J(1,3).eq.0) then
                    dump=.true.
                  else
                    gb_ty_np1(i,j,k)=gb_ty_np1(i,j,k)
     &                               -efe(1,3)/efe_J(1,3)
                  end if

                  if (is_nan(efe(1,4)).or.is_nan(efe_J(1,4)).or.
     &              efe_J(1,4).eq.0) then
                    dump=.true.
                  else
                    gb_tz_np1(i,j,k)=gb_tz_np1(i,j,k)
     &                               -efe(1,4)/efe_J(1,4)
                  end if


                  if (is_nan(efe(2,2)).or.is_nan(efe_J(2,2)).or.
     &             efe_J(2,2).eq.0) then
                    dump=.true.
                  else
                    gb_xx_np1(i,j,k)=gb_xx_np1(i,j,k)
     &                               -efe(2,2)/efe_J(2,2)
                  end if

                  if (is_nan(efe(2,3)).or.is_nan(efe_J(2,3)).or.
     &              efe_J(2,3).eq.0) then
                    dump=.true.
                  else
                    gb_xy_np1(i,j,k)=gb_xy_np1(i,j,k)
     &                               -efe(2,3)/efe_J(2,3)
                  end if

                  if (is_nan(efe(2,4)).or.is_nan(efe_J(2,4)).or.
     &              efe_J(2,4).eq.0) then
                    dump=.true.
                  else
                    gb_xz_np1(i,j,k)=gb_xz_np1(i,j,k)
     &                               -efe(2,4)/efe_J(2,4)
                  end if

                  if (is_nan(efe(3,3)).or.is_nan(efe_J(3,3)).or.
     &              efe_J(3,3).eq.0) then
                    dump=.true.
                  else
                    gb_yy_np1(i,j,k)=gb_yy_np1(i,j,k)
     &                               -efe(3,3)/efe_J(3,3)
                  end if

                  if (is_nan(efe(3,4)).or.is_nan(efe_J(3,4)).or.
     &              efe_J(3,4).eq.0) then
                    dump=.true.
                  else
                    gb_yz_np1(i,j,k)=gb_yz_np1(i,j,k)
     &                               -efe(3,4)/efe_J(3,4)
                  end if

                  if (is_nan(efe(4,4)).or.is_nan(efe_J(4,4)).or.
     &              efe_J(4,4).eq.0) then
                    dump=.true.
                  else
                    gb_zz_np1(i,j,k)=gb_zz_np1(i,j,k)
     &                               -efe(4,4)/efe_J(4,4)
                  end if
                end if

                ! update phi1
                if (is_nan(phi1_res).or.is_nan(phi1_J)) then
                  dump=.true.
                else
                  phi1_np1(i,j,k)=phi1_np1(i,j,k)-phi1_res/phi1_J
                end if

                gb_res(i,j,k) = 
     &            max(abs(efe(1,1)/efe_J(1,1)),
     &                abs(efe(1,2)/efe_J(1,2)),
     &                abs(efe(1,3)/efe_J(1,3)),
     &                abs(efe(1,4)/efe_J(1,4)),
     &                abs(efe(2,2)/efe_J(2,2)),
     &                abs(efe(2,3)/efe_J(2,3)),
     &                abs(efe(2,4)/efe_J(2,4)),
     &                abs(efe(3,3)/efe_J(3,3)),
     &                abs(efe(3,4)/efe_J(3,4)),
     &                abs(efe(4,4)/efe_J(4,4)))
                kg_res(i,j,k)=abs(phi1_res/phi1_J)

                ! save pointwise max of constraint violation
                cl_res(i,j,k)=
     &            max(abs(c_l(1)),abs(c_l(2)),abs(c_l(3)),
     &                abs(c_l(4)))

                if (dump.and.first_nan.or.
     &              (ltrace.and.abs(x(i)).lt.0.1.and.abs(y(j)).lt.0.1)
     &             ) then
                  first_nan=.false.
                  write(*,*)
                  write(*,*) 'g_evo_opt: Nan/zero at i,j,k,Nx,Ny,Nz,dx='
     &                                              ,i,j,k,Nx,Ny,Nz,dx
                  write(*,*) 'x,y,z=',x(i),y(j),z(k)
                  write(*,*) 'dt,dx,dy,dz=',dt,dx,dy,dz
                  write(*,*) 'x0,y0,z0,rho0=',x0,y0,z0,rho0


                  write(*,*) ' at tn:'
                  write(*,*) ' gb_tt np1,n,nm1:',gb_tt_np1(i,j,k),
     &                   gb_tt_n(i,j,k),gb_tt_nm1(i,j,k)
                  write(*,*) ' gb_tx np1,n,nm1:',gb_tx_np1(i,j,k),
     &                   gb_tx_n(i,j,k),gb_tx_nm1(i,j,k)
                  write(*,*) ' gb_ty np1,n,nm1:',gb_ty_np1(i,j,k),
     &                   gb_ty_n(i,j,k),gb_ty_nm1(i,j,k)
                  write(*,*) ' gb_tz np1,n,nm1:',gb_tz_np1(i,j,k),
     &                   gb_tz_n(i,j,k),gb_tz_nm1(i,j,k)
                  write(*,*) ' gb_xx np1,n,nm1:',gb_xx_np1(i,j,k),
     &                   gb_xx_n(i,j,k),gb_xx_nm1(i,j,k)
                  write(*,*) ' gb_xy np1,n,nm1:',gb_xy_np1(i,j,k),
     &                   gb_xy_n(i,j,k),gb_xy_nm1(i,j,k)
                  write(*,*) ' gb_xz np1,n,nm1:',gb_xz_np1(i,j,k),
     &                   gb_xz_n(i,j,k),gb_xz_nm1(i,j,k)
                  write(*,*) ' gb_yy np1,n,nm1:',gb_yy_np1(i,j,k),
     &                   gb_yy_n(i,j,k),gb_yy_nm1(i,j,k)
                  write(*,*) ' gb_yz np1,n,nm1:',gb_yz_np1(i,j,k),
     &                   gb_yz_n(i,j,k),gb_yz_nm1(i,j,k)
                  write(*,*) ' gb_zz np1,n,nm1:',gb_zz_np1(i,j,k),
     &                   gb_zz_n(i,j,k),gb_zz_nm1(i,j,k)

                  write(*,*) ' gads_tt :',gads_ll(1,1)
                  write(*,*) ' gads_tx :',gads_ll(1,2)
                  write(*,*) ' gads_ty :',gads_ll(1,3)
                  write(*,*) ' gads_tz :',gads_ll(1,4)
                  write(*,*) ' gads_xx :',gads_ll(2,2)
                  write(*,*) ' gads_xy :',gads_ll(2,3)
                  write(*,*) ' gads_xz :',gads_ll(2,4)
                  write(*,*) ' gads_yy :',gads_ll(3,3)
                  write(*,*) ' gads_yz :',gads_ll(3,4)
                  write(*,*) ' gads_zz:',gads_ll(4,4)
                  write(*,*) ' h0_tt :',h0_ll(1,1)
                  write(*,*) ' h0_tx :',h0_ll(1,2)
                  write(*,*) ' h0_ty :',h0_ll(1,3)
                  write(*,*) ' h0_tz :',h0_ll(1,4)
                  write(*,*) ' h0_xx :',h0_ll(2,2)
                  write(*,*) ' h0_xy :',h0_ll(2,3)
                  write(*,*) ' h0_xz :',h0_ll(2,4)
                  write(*,*) ' h0_yy :',h0_ll(3,3)
                  write(*,*) ' h0_yz :',h0_ll(3,4)
                  write(*,*) ' h0_zz:',h0_ll(4,4)
                  write(*,*) ' g0_tt :',g0_ll(1,1)
                  write(*,*) ' g0_tx :',g0_ll(1,2)
                  write(*,*) ' g0_ty :',g0_ll(1,3)
                  write(*,*) ' g0_tz :',g0_ll(1,4)
                  write(*,*) ' g0_xx :',g0_ll(2,2)
                  write(*,*) ' g0_xy :',g0_ll(2,3)
                  write(*,*) ' g0_xz :',g0_ll(2,4)
                  write(*,*) ' g0_yy :',g0_ll(3,3)
                  write(*,*) ' g0_yz :',g0_ll(3,4)
                  write(*,*) ' g0_zz:',g0_ll(4,4)
                  write(*,*) ' g0u_tt :',g0_uu(1,1)
                  write(*,*) ' g0u_tx :',g0_uu(1,2)
                  write(*,*) ' g0u_ty :',g0_uu(1,3)
                  write(*,*) ' g0u_tz :',g0_uu(1,4)
                  write(*,*) ' g0u_xx :',g0_uu(2,2)
                  write(*,*) ' g0u_xy :',g0_uu(2,3)
                  write(*,*) ' g0u_xz :',g0_uu(2,4)
                  write(*,*) ' g0u_yy :',g0_uu(3,3)
                  write(*,*) ' g0u_yz :',g0_uu(3,4)
                  write(*,*) ' g0u_zz:',g0_uu(4,4)
                  write(*,*) ' cd_tt:',cd_ll(1,1)
                  write(*,*) ' cd_tx:',cd_ll(1,2)
                  write(*,*) ' cd_ty:',cd_ll(1,3)
                  write(*,*) ' cd_tz:',cd_ll(1,4)
                  write(*,*) ' cd_xx:',cd_ll(2,2)
                  write(*,*) ' cd_xy:',cd_ll(2,3)
                  write(*,*) ' cd_xz:',cd_ll(2,4)
                  write(*,*) ' cd_yy:',cd_ll(3,3)
                  write(*,*) ' cd_yz:',cd_ll(3,4)
                  write(*,*) ' cd_zz:',cd_ll(4,4)
                  write(*,*) ' phi1:',phi1_n(i,j,k)
                  write(*,*) ' phi np1,n,nm1:',phi1_np1(i,j,k),
     &                     phi1_n(i,j,k),phi1_nm1(i,j,k)
                  write(*,*) ' res J:'
                  write(*,*) ' tt:',efe(1,1),efe_J(1,1)
                  write(*,*) ' tx:',efe(1,2),efe_J(1,2)
                  write(*,*) ' ty:',efe(1,3),efe_J(1,3)
                  write(*,*) ' tz:',efe(1,4),efe_J(1,4)
                  write(*,*) ' xx:',efe(2,2),efe_J(2,2)
                  write(*,*) ' xy:',efe(2,3),efe_J(2,3)
                  write(*,*) ' xz:',efe(2,4),efe_J(2,4)
                  write(*,*) ' yy:',efe(3,3),efe_J(3,3)
                  write(*,*) ' yz:',efe(3,4),efe_J(3,4)
                  write(*,*) ' zz:',efe(4,4),efe_J(4,4)
                  write(*,*) ' phi1:',phi1_res,phi1_J
                end if

              ! (REGION) next-to-ads-bdy points; set by linear interpolation
              else if (chr2(i,j,k).eq.ex) then
!        write (*,*) 'INTERP_FROM_ADS_BDY'
!        write (*,*) 'i,j,k,x(i),y(j)=',i,j,k,x(i),y(j)
                call interp_from_ads_bdy(gb_tt_np1,x,y,z,L,i,j,k,chr,ex,
     &                 Nx,Ny,Nz)
                call interp_from_ads_bdy(gb_tx_np1,x,y,z,L,i,j,k,chr,ex,
     &                 Nx,Ny,Nz)
                call interp_from_ads_bdy(gb_ty_np1,x,y,z,L,i,j,k,chr,ex,
     &                 Nx,Ny,Nz)
                call interp_from_ads_bdy(gb_tz_np1,x,y,z,L,i,j,k,chr,ex,
     &                 Nx,Ny,Nz)
                call interp_from_ads_bdy(gb_xx_np1,x,y,z,L,i,j,k,chr,ex,
     &                 Nx,Ny,Nz)
                call interp_from_ads_bdy(gb_xy_np1,x,y,z,L,i,j,k,chr,ex,
     &                 Nx,Ny,Nz)
                call interp_from_ads_bdy(gb_xz_np1,x,y,z,L,i,j,k,chr,ex,
     &                 Nx,Ny,Nz)
                call interp_from_ads_bdy(gb_yy_np1,x,y,z,L,i,j,k,chr,ex,
     &                 Nx,Ny,Nz)
                call interp_from_ads_bdy(gb_yz_np1,x,y,z,L,i,j,k,chr,ex,
     &                 Nx,Ny,Nz)
                call interp_from_ads_bdy(gb_zz_np1,x,y,z,L,i,j,k,chr,ex,
     &                 Nx,Ny,Nz)
                call interp_from_ads_bdy(Hb_t_np1,x,y,z,L,i,j,k,chr,ex,
     &                 Nx,Ny,Nz)
                call interp_from_ads_bdy(Hb_x_np1,x,y,z,L,i,j,k,chr,ex,
     &                 Nx,Ny,Nz)
                call interp_from_ads_bdy(Hb_y_np1,x,y,z,L,i,j,k,chr,ex,
     &                 Nx,Ny,Nz)
                call interp_from_ads_bdy(Hb_z_np1,x,y,z,L,i,j,k,chr,ex,
     &                 Nx,Ny,Nz)
                call interp_from_ads_bdy(phi1_np1,x,y,z,L,i,j,k,
     &                    chr,ex,Nx,Ny,Nz)
!                 gb_tx_np1(i,j,k)=0
!                 gb_ty_np1(i,j,k)=0
!                 gb_tz_np1(i,j,k)=0
!                 gb_xx_np1(i,j,k)=0
!                 gb_xy_np1(i,j,k)=0
!                 gb_xz_np1(i,j,k)=0
!                 gb_yy_np1(i,j,k)=0
!                 gb_yz_np1(i,j,k)=0
!                 gb_zz_np1(i,j,k)=0
!                 Hb_t_np1(i,j,k)=0
!                 Hb_x_np1(i,j,k)=0
!                 Hb_y_np1(i,j,k)=0
!                 Hb_z_np1(i,j,k)=0
!                 phi1_np1(i,j,k)=0
!                gb_tt_np1(i,j,k)=gb_xx_np1(i,j,k)+gb_yy_np1(i,j,k)
!     &                           +gb_zz_np1(i,j,k)  !CHECK

              ! (REGION) non-interior points; set to zero in prior to applying bcs 
              else 
                gb_tt_np1(i,j,k) = 0
                gb_tx_np1(i,j,k) = 0
                gb_ty_np1(i,j,k) = 0
                gb_tz_np1(i,j,k) = 0
                gb_xx_np1(i,j,k) = 0
                gb_xy_np1(i,j,k) = 0
                gb_xz_np1(i,j,k) = 0
                gb_yy_np1(i,j,k) = 0
                gb_yz_np1(i,j,k) = 0
                gb_zz_np1(i,j,k) = 0 
                phi1_np1(i,j,k) = 0 
                gb_res(i,j,k) = 0

              endif ! (near start of main loop)

            end do
          end do
         end do
        end do

        ! (REGION) y=0 axis; impose Neumann bcs by 2-pt regularization 
!        call axi_reg_g(gb_tt_np1,gb_tx_np1,gb_ty_np1,gb_xx_np1,
!     &                 gb_xy_np1,gb_yy_np1,psi_np1,tfunction,chr,ex,
!     &                 L,x,y,z,Nx,Ny,Nz,regtype)
!        call axi_reg_Hb(Hb_t_np1,Hb_x_np1,Hb_y_np1,
!     &                  chr,ex,L,x,y,z,Nx,Ny,Nz,regtype)
!        call axi_reg_phi(phi1_np1,chr,ex,L,x,y,z,Nx,Ny,Nz,regtype)

        return
        end

c----------------------------------------------------------------------
        logical function is_nan(x)
        implicit none
        real*8 x

        integer is_a_nan

        call check_nan(x,is_a_nan)

        if (is_a_nan.eq.0) then
           is_nan=.false.
        else
           is_nan=.true.
        end if

        return
        end
c----------------------------------------------------------------------

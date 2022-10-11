c----------------------------------------------------------------------
c routine for computing independent residuals 
c and Generalized Harmonic Constraints C_a:=g_{ab} (H^a-\Box x^a)
c----------------------------------------------------------------------
        subroutine ires_ghconstr(efe_all_ires,
     &                  efe_tt_ires,efe_tx_ires,efe_ty_ires,
     &                  efe_tz_ires,
     &                  efe_xx_ires,efe_xy_ires,
     &                  efe_xz_ires,
     &                  efe_yy_ires,
     &                  efe_yz_ires,
     &                  efe_zz_ires,
     &                  kg_ires,
     &                  ghconstr_all,
     &                  ghconstr_t,
     &                  ghconstr_x,
     &                  ghconstr_y,
     &                  ghconstr_z,
     &                  gb_tt_np1,gb_tt_n,gb_tt_nm1,
     &                  gb_tx_np1,gb_tx_n,gb_tx_nm1,
     &                  gb_ty_np1,gb_ty_n,gb_ty_nm1,
     &                  gb_tz_np1,gb_tz_n,gb_tz_nm1,
     &                  gb_xx_np1,gb_xx_n,gb_xx_nm1,
     &                  gb_xy_np1,gb_xy_n,gb_xy_nm1,
     &                  gb_xz_np1,gb_xz_n,gb_xz_nm1,
     &                  gb_yy_np1,gb_yy_n,gb_yy_nm1,
     &                  gb_yz_np1,gb_yz_n,gb_yz_nm1,
     &                  gb_zz_np1,gb_zz_n,gb_zz_nm1,
     &                  Hb_t_np1,Hb_t_n,Hb_t_nm1,
     &                  Hb_x_np1,Hb_x_n,Hb_x_nm1,
     &                  Hb_y_np1,Hb_y_n,Hb_y_nm1,
     &                  Hb_z_np1,Hb_z_n,Hb_z_nm1,
     &                  phi1_np1,phi1_n,phi1_nm1,
     &                  x,y,z,dt,chr,L,ex,Nx,Ny,Nz,phys_bdy,ghost_width,
     &                  ief_bh_r0,a_rot,kerrads_background)

        implicit none
        real*8 ief_bh_r0,a_rot
        integer kerrads_background
        logical calc_der,calc_adv_quant
        data calc_der/.true./
        data calc_adv_quant/.false./
        integer Nx,Ny,Nz
        integer i,j,k
        integer phys_bdy(6),ghost_width(6)
        real*8 efe_all_ires(Nx,Ny,Nz)
        real*8 efe_tt_ires(Nx,Ny,Nz),efe_tx_ires(Nx,Ny,Nz)
        real*8 efe_ty_ires(Nx,Ny,Nz)
        real*8 efe_tz_ires(Nx,Ny,Nz)
        real*8 efe_xx_ires(Nx,Ny,Nz),efe_xy_ires(Nx,Ny,Nz)
        real*8 efe_xz_ires(Nx,Ny,Nz)
        real*8 efe_yy_ires(Nx,Ny,Nz)
        real*8 efe_yz_ires(Nx,Ny,Nz)
        real*8 efe_zz_ires(Nx,Ny,Nz)
        real*8 kg_ires(Nx,Ny,Nz)
        real*8 ghconstr_all(Nx,Ny,Nz)
        real*8 ghconstr_t(Nx,Ny,Nz)
        real*8 ghconstr_x(Nx,Ny,Nz)
        real*8 ghconstr_y(Nx,Ny,Nz)
        real*8 ghconstr_z(Nx,Ny,Nz)
        real*8 chr(Nx,Ny,Nz),ex
        real*8 x(Nx),y(Ny),z(Nz),dt,L
        real*8 lambda4
        real*8 phi1_np1(Nx,Ny,Nz),phi1_n(Nx,Ny,Nz),phi1_nm1(Nx,Ny,Nz)
        real*8 gb_tt_np1(Nx,Ny,Nz),gb_tx_np1(Nx,Ny,Nz)
        real*8 gb_ty_np1(Nx,Ny,Nz)
        real*8 gb_tz_np1(Nx,Ny,Nz)
        real*8 gb_xx_np1(Nx,Ny,Nz),gb_xy_np1(Nx,Ny,Nz)
        real*8 gb_xz_np1(Nx,Ny,Nz)
        real*8 gb_yy_np1(Nx,Ny,Nz)
        real*8 gb_yz_np1(Nx,Ny,Nz)
        real*8 gb_zz_np1(Nx,Ny,Nz)
        real*8 gb_tt_n(Nx,Ny,Nz),gb_tx_n(Nx,Ny,Nz),gb_ty_n(Nx,Ny,Nz)
        real*8 gb_tz_n(Nx,Ny,Nz)
        real*8 gb_xx_n(Nx,Ny,Nz),gb_xy_n(Nx,Ny,Nz),gb_yy_n(Nx,Ny,Nz)
        real*8 gb_xz_n(Nx,Ny,Nz)
        real*8 gb_yz_n(Nx,Ny,Nz)
        real*8 gb_zz_n(Nx,Ny,Nz)
        real*8 gb_tt_nm1(Nx,Ny,Nz),gb_tx_nm1(Nx,Ny,Nz)
        real*8 gb_ty_nm1(Nx,Ny,Nz)
        real*8 gb_tz_nm1(Nx,Ny,Nz)
        real*8 gb_xx_nm1(Nx,Ny,Nz),gb_xy_nm1(Nx,Ny,Nz)
        real*8 gb_xz_nm1(Nx,Ny,Nz)
        real*8 gb_yy_nm1(Nx,Ny,Nz)
        real*8 gb_yz_nm1(Nx,Ny,Nz)
        real*8 gb_zz_nm1(Nx,Ny,Nz)

        real*8 Hb_t_np1(Nx,Ny,Nz),Hb_t_n(Nx,Ny,Nz),Hb_t_nm1(Nx,Ny,Nz)
        real*8 Hb_x_np1(Nx,Ny,Nz),Hb_x_n(Nx,Ny,Nz),Hb_x_nm1(Nx,Ny,Nz)
        real*8 Hb_y_np1(Nx,Ny,Nz),Hb_y_n(Nx,Ny,Nz),Hb_y_nm1(Nx,Ny,Nz)
        real*8 Hb_z_np1(Nx,Ny,Nz),Hb_z_n(Nx,Ny,Nz),Hb_z_nm1(Nx,Ny,Nz)

        integer is,ie,js,je,ks,ke

        integer i1,j1,k1,a,b,c,d,e,f,g,h
        integer ic,jc,kc
        real*8 efe_ires(4,4)

        real*8 PI
        parameter (PI=3.141592653589793d0)

        real*8 dx,dy,dz
        real*8 x0,y0,z0,rho0        

        real*8 boxx_u(4),boxx_l(4) 

        !--------------------------------------------------------------
        ! variables for tensor manipulations 
        !(indices are t,x,y,theta,phi)
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
        real*8 phi10,phi10_x(4),phi10_xx(4,4)
        real*8 phi_x(4),phi_xx(4,4)

        !--------------------------------------------------------------
        ! variables for outward null expansion in spherical symmetry
        !--------------------------------------------------------------
        real*8 n_l(4),s_l(4)
        real*8 n_u(4),s_u(4)
        real*8 n_l_x(4,4),n_u_x(4,4),s_l_x(4,4)
        real*8 f0_x(4)
        real*8 f0_xx(4,4)
        real*8 gam_uu(4,4),sig_uu(4,4)
        real*8 gam_uu_x(4,4,4)
        real*8 normsusq

        real*8 g0gamfx(4)
        real*8 nufx,nuxfx(4),gamxfxfx(4)

        real*8 theta(Nx,Ny,Nz)

!!!!!!!DEBUGGING!!!!!!
        integer max_i,max_j,max_k
        real*8  max_efe_all_ires
!!!!!!!!!!!!!!!!!!!!!1

        ! initialize fixed-size variables
        data i,j,k,is,ie,js,je,ks,ke/0,0,0,0,0,0,0,0,0/
        data ic,jc,kc/0,0,0/
        data i1,j1,k1,a,b,c,d,e/0,0,0,0,0,0,0,0/

        data dx,dy,dz/0.0,0.0,0.0/
        data x0,y0,rho0/0.0,0.0,0.0/    

        data g0_ll,g0_uu/16*0.0,16*0.0/
        data gads_ll,gads_uu/16*0.0,16*0.0/
        data h0_ll,h0_uu/16*0.0,16*0.0/
        data gamma_ull/64*0.0/
        data gamma_ull_x/256*0.0/

        data g0_ll_x,g0_uu_x/64*0.0,64*0.0/
        data gads_ll_x,gads_uu_x/64*0.0,64*0.0/
        data h0_ll_x,h0_uu_x/64*0.0,64*0.0/

        data g0_ll_xx/256*0.0/
        data gads_ll_xx/256*0.0/
        data h0_ll_xx/256*0.0/

        data ricci/0.0/
        data ricci_ll,ricci_lu/16*0.0,16*0.0/
        data einstein_ll,set_ll/16*0.0,16*0.0/
        data riemann_ulll/256*0.0/

        data A_l,Hads_l/4*0.0,4*0.0/
        data A_l_x/16*0.0/

        data phi10_x/4*0.0/
        data phi10_xx/16*0.0/
        data phi_x/4*0.0/
        data phi_xx/16*0.0/

        data boxx_u,boxx_l/4*0.0,4*0.0/

!----------------------------------------------------------------------


        dx=(x(2)-x(1))
        dy=(y(2)-y(1))
        dz=(z(2)-z(1))

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

        ! adjust index bounds to compensate for ghost_width
        if (ghost_width(1).gt.0) is=is+ghost_width(1)-1
        if (ghost_width(2).gt.0) ie=ie-(ghost_width(2)-1)
        if (ghost_width(3).gt.0) js=js+ghost_width(3)-1
        if (ghost_width(4).gt.0) je=je-(ghost_width(4)-1)
        if (ghost_width(5).gt.0) ks=ks+ghost_width(5)-1
        if (ghost_width(6).gt.0) ke=ke-(ghost_width(6)-1)

        ! (MAIN LOOP) loop through spacetime points x(i),y(j)
        do i=is,ie
          do j=js,je
           do k=ks,ke

            x0=x(i)
            y0=y(j)
            z0=z(k)
            rho0=sqrt(x0**2+y0**2+z0**2)

            if (chr(i,j,k).ne.ex) then

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

!       write (*,*) 'L,i,j,k,x0,y0,z0,rho0=',L,i,j,k,x0,y0,z0,rho0
!       write (*,*) 'g0_ll_x(1,1,1)=',g0_ll_x(1,1,1)
!       write (*,*) 'g0_ll_x(1,1,2)=',g0_ll_x(1,1,2)
!       write (*,*) 'g0_ll_x(1,1,3)=',g0_ll_x(1,1,3)
!       write (*,*) 'g0_ll_x(1,1,4)=',g0_ll_x(1,1,4)
!       write (*,*) 'g0_ll_x(1,1,1,1)=',g0_ll_xx(1,1,1,1)
!       write (*,*) 'g0_ll_x(1,1,1,2)=',g0_ll_xx(1,1,1,2)
!       write (*,*) 'g0_ll_x(1,1,1,3)=',g0_ll_xx(1,1,1,3)
!       write (*,*) 'g0_ll_x(1,1,1,4)=',g0_ll_xx(1,1,1,4)
!       write (*,*) 'g0_ll_x(1,1,2,2)=',g0_ll_xx(1,1,2,2)
!       write (*,*) 'g0_ll_x(1,1,2,3)=',g0_ll_xx(1,1,2,3)
!       write (*,*) 'g0_ll_x(1,1,2,4)=',g0_ll_xx(1,1,2,4)
!       write (*,*) 'g0_ll_x(1,1,3,3)=',g0_ll_xx(1,1,3,3)
!       write (*,*) 'g0_ll_x(1,1,3,4)=',g0_ll_xx(1,1,3,4)
!       write (*,*) 'g0_ll_x(1,1,4,4)=',g0_ll_xx(1,1,4,4)
!
!       write (*,*) 'gb_tt_nm1(i,j,k),gb_tt_n(i,j,k),gb_tt_np1(i,j,k)='
!     &             ,gb_tt_nm1(i,j,k),gb_tt_n(i,j,k),gb_tt_np1(i,j,k)
!       write (*,*) 'gb_tx_nm1(i,j,k),gb_tx_n(i,j,k),gb_tx_np1(i,j,k)='
!     &             ,gb_tx_nm1(i,j,k),gb_tx_n(i,j,k),gb_tx_np1(i,j,k)
!       write (*,*) 'gb_ty_nm1(i,j,k),gb_ty_n(i,j,k),gb_ty_np1(i,j,k)='
!     &             ,gb_ty_nm1(i,j,k),gb_ty_n(i,j,k),gb_ty_np1(i,j,k)
!       write (*,*) 'gb_tz_nm1(i,j,k),gb_tz_n(i,j,k),gb_tz_np1(i,j,k)='
!     &             ,gb_tz_nm1(i,j,k),gb_tz_n(i,j,k),gb_tz_np1(i,j,k)
!       write (*,*) 'gb_xx_nm1(i,j,k),gb_xx_n(i,j,k),gb_xx_np1(i,j,k)='
!     &             ,gb_xx_nm1(i,j,k),gb_xx_n(i,j,k),gb_xx_np1(i,j,k)
!       write (*,*) 'gb_xy_nm1(i,j,k),gb_xy_n(i,j,k),gb_xy_np1(i,j,k)='
!     &             ,gb_xy_nm1(i,j,k),gb_xy_n(i,j,k),gb_xy_np1(i,j,k)
!       write (*,*) 'gb_xz_nm1(i,j,k),gb_xz_n(i,j,k),gb_xz_np1(i,j,k)='
!     &             ,gb_xz_nm1(i,j,k),gb_xz_n(i,j,k),gb_xz_np1(i,j,k)
!       write (*,*) 'gb_yy_nm1(i,j,k),gb_yy_n(i,j,k),gb_yy_np1(i,j,k)='
!     &             ,gb_yy_nm1(i,j,k),gb_yy_n(i,j,k),gb_yy_np1(i,j,k)
!       write (*,*) 'gb_yz_nm1(i,j,k),gb_yz_n(i,j,k),gb_yz_np1(i,j,k)='
!     &             ,gb_yz_nm1(i,j,k),gb_yz_n(i,j,k),gb_yz_np1(i,j,k)
!       write (*,*) 'gb_zz_nm1(i,j,k),gb_zz_n(i,j,k),gb_zz_np1(i,j,k)='
!     &             ,gb_zz_nm1(i,j,k),gb_zz_n(i,j,k),gb_zz_np1(i,j,k)

              ! calculates efe_ires functions at point i,j
              !(efe_ires_ab=G_ab+lambda4*g_ab-8*PI*T_ab)
              do a=1,4
                do b=a,4
                  efe_ires(a,b)=einstein_ll(a,b)+lambda4*g0_ll(a,b)
     &                                          -8*PI*set_ll(a,b)
                end do
              end do
              efe_tt_ires(i,j,k)=efe_ires(1,1) 
              efe_tx_ires(i,j,k)=efe_ires(1,2) 
              efe_ty_ires(i,j,k)=efe_ires(1,3) 
              efe_tz_ires(i,j,k)=efe_ires(1,4)
              efe_xx_ires(i,j,k)=efe_ires(2,2)
              efe_xy_ires(i,j,k)=efe_ires(2,3)
              efe_xz_ires(i,j,k)=efe_ires(2,4)
              efe_yy_ires(i,j,k)=efe_ires(3,3)
              efe_yz_ires(i,j,k)=efe_ires(3,4)
              efe_zz_ires(i,j,k)=efe_ires(4,4)

              ! calculate efe_all_ires function at point i,j,k
              efe_all_ires(i,j,k)=
     &         max(abs(efe_tt_ires(i,j,k)),
     &            abs(efe_tx_ires(i,j,k)),
     &            abs(efe_ty_ires(i,j,k)),
     &            abs(efe_tz_ires(i,j,k)),
     &            abs(efe_xx_ires(i,j,k)),
     &            abs(efe_xy_ires(i,j,k)),
     &            abs(efe_xz_ires(i,j,k)),
     &            abs(efe_yy_ires(i,j,k)),
     &            abs(efe_yz_ires(i,j,k)),
     &            abs(efe_zz_ires(i,j,k)))

              ! calculate boxx^c at point i,j,k
              ! (boxx^c = -g^ab gamma^c_ab)
              do c=1,4
                boxx_u(c)=-( gamma_ull(c,1,1)*g0_uu(1,1)+
     &                       gamma_ull(c,2,2)*g0_uu(2,2)+
     &                       gamma_ull(c,3,3)*g0_uu(3,3)+
     &                       gamma_ull(c,4,4)*g0_uu(4,4)+
     &                    2*(gamma_ull(c,1,2)*g0_uu(1,2)+
     &                       gamma_ull(c,1,3)*g0_uu(1,3)+
     &                       gamma_ull(c,1,4)*g0_uu(1,4)+
     &                       gamma_ull(c,2,3)*g0_uu(2,3)+
     &                       gamma_ull(c,2,4)*g0_uu(2,4)+
     &                       gamma_ull(c,3,4)*g0_uu(3,4)) )
              end do

              ! calculate boxx_a at point i,j,k
              ! (boxx_a = g_ab boxx^b)
              do a=1,4
                boxx_l(a)=boxx_u(1)*g0_ll(a,1)+
     &                    boxx_u(2)*g0_ll(a,2)+
     &                    boxx_u(3)*g0_ll(a,3)+
     &                    boxx_u(4)*g0_ll(a,4)
              end do

              !calculate generalized harmonic constraints at point i,j,k
              ghconstr_t(i,j,k)=Hads_l(1)+A_l(1)-boxx_l(1)
              ghconstr_x(i,j,k)=Hads_l(2)+A_l(2)-boxx_l(2)
              ghconstr_y(i,j,k)=Hads_l(3)+A_l(3)-boxx_l(3)
              ghconstr_z(i,j,k)=Hads_l(4)+A_l(4)-boxx_l(4)

              ghconstr_all(i,j,k)=
     &        max(abs(ghconstr_t(i,j,k)),
     &            abs(ghconstr_x(i,j,k)),
     &            abs(ghconstr_y(i,j,k)),
     &            abs(ghconstr_z(i,j,k)))


!              ! define unit time-like vector n, normal to t=const
!              ! surfaces
!              n_l(1)=-1/sqrt(-g0_uu(1,1))
!              do a=1,4
!                n_u(a)=n_l(1)*g0_uu(a,1)+
!     &                 n_l(2)*g0_uu(a,2)+
!     &                 n_l(3)*g0_uu(a,3)+
!     &                 n_l(4)*g0_uu(a,4)
!              end do
!              do b=1,4
!                n_l_x(1,b)=-1/2.0d0/sqrt(-g0_uu(1,1))**3*g0_uu_x(1,1,b)
!              end do
!              do a=1,4
!                do b=1,4
!                  n_u_x(a,b)=n_l_x(1,b)*g0_uu(a,1)+
!     &                       n_l_x(2,b)*g0_uu(a,2)+
!     &                       n_l_x(3,b)*g0_uu(a,3)+
!     &                       n_l_x(4,b)*g0_uu(a,4)+
!     &                       n_l(1)*g0_uu_x(a,1,b)+
!     &                       n_l(2)*g0_uu_x(a,2,b)+
!     &                       n_l(3)*g0_uu_x(a,3,b)+
!     &                       n_l(4)*g0_uu_x(a,4,b)
!                end do
!              end do
!
!              ! define gradients of the flow field f=r-AH_R(chi,phi) 
!              ! NOTE: CHECK THESE WITH MATHEMATICA GIVEN
!              ! f(x,y)=sqrt(x^2+y^2)
!              f0_x(1)=0
!              f0_x(2)=x0/rho0
!              f0_x(3)=y0/rho0
!              f0_x(4)=z0/rho0
!              f0_xx(1,1)=0
!              f0_xx(1,2)=0
!              f0_xx(1,3)=0
!              f0_xx(1,4)=0
!              f0_xx(2,2)=(y0**2+z0**2)/rho0**3
!              f0_xx(2,3)=-x0*y0/rho0**3
!              f0_xx(2,4)=-x0*z0/rho0**3
!              f0_xx(3,3)=(x0**2+z0**2)/rho0**3
!              f0_xx(3,4)=-y0*z0/rho0**3
!              f0_xx(4,4)=(x0**2+y0**2)/rho0**3
!
!              do a=1,3
!                do b=a+1,4
!                  f0_xx(b,a)=f0_xx(a,b)
!                end do
!              end do
!
!              ! define metric on codimension-1 surfaces
!              do a=1,4
!                do b=1,4
!                  gam_uu(a,b)=g0_uu(a,b)+n_u(a)*n_u(b)
!                end do
!              end do
!              do a=1,4
!                do b=1,4
!                  do c=1,4
!                    gam_uu_x(a,b,c)=g0_uu_x(a,b,c)
!     &                             +n_u_x(a,c)*n_u(b)
!     &                             +n_u(a)*n_u_x(b,c)
!                  end do
!                end do
!              end do
!
!              ! define unit space-like vector s, orthogonal to n and
!              ! projected gradient of the flow field f
!              do a=1,4
!                s_u(a)=0.0d0
!                do b=1,4
!                  s_u(a)=s_u(a)+gam_uu(a,b)*f0_x(b)
!                end do
!              end do
!              normsusq=0.0d0
!              do a=1,4
!                do b=1,4
!                  normsusq=normsusq+gam_uu(a,b)*f0_x(a)*f0_x(b)
!                end do
!              end do
!              do a=1,4
!                s_u(a)=s_u(a)/sqrt(normsusq)
!              end do
!              do a=1,4
!                s_l(a)=s_u(1)*g0_ll(a,1)+
!     &                 s_u(2)*g0_ll(a,2)+
!     &                 s_u(3)*g0_ll(a,3)+
!     &                 s_u(4)*g0_ll(a,4)
!              end do
!
!              nufx=0
!              do a=1,4
!                nufx=nufx
!     &              +n_u(a)*f0_x(a)
!                nuxfx(a)=0
!                gamxfxfx(a)=0
!                do c=1,4
!                  nuxfx(a)=nuxfx(a)
!     &                    +n_u_x(c,a)*f0_x(c)
!                  do d=1,4
!                    gamxfxfx(a)=gamxfxfx(a)
!     &                        +gam_uu_x(c,d,a)*f0_x(c)*f0_x(d)
!     &                        +gam_uu(c,d)*f0_xx(c,a)*f0_x(d)
!     &                        +gam_uu(c,d)*f0_x(c)*f0_xx(d,a)
!                  end do
!                end do
!              end do
!              do a=1,4
!                do b=1,4
!                  s_l_x(a,b)=
!     &                     (f0_xx(a,b)+n_l_x(a,b)*nufx+n_l(a)*nuxfx(b))
!     &                     /sqrt(normsusq)
!     &                    -(f0_x(a)+n_l(a)*nufx)*gamxfxfx(b)
!     &                     /2.0d0/sqrt(normsusq)**3
!                end do
!              end do
!
!              ! define metric on codimension-2 surfaces
!              do a=1,4
!                do b=1,4
!                  sig_uu(a,b)=g0_uu(a,b)+n_u(a)*n_u(b)-s_u(a)*s_u(b)
!                end do
!              end do
!
!              ! for theta: outward null expansion
!              theta(i,j,k)=0.0d0
!              do c=1,4
!                do d=1,4
!                  theta(i,j,k)=theta(i,j,k)
!     &                   +sig_uu(c,d)*(n_l_x(c,d)+s_l_x(c,d))
!                  do e=1,4
!                    theta(i,j,k)=theta(i,j,k)
!     &                   -sig_uu(c,d)*gamma_ull(e,c,d)*(n_l(e)+s_l(e))
!                  end do
!                end do
!              end do

!              efe_tt_ires(i,j,k)=!Hads_l(1)+A_l(1)-boxx_l(1)
!     &           sqrt((-(-1+rho0**2)**6*gb_xy_n(i,j,k)**2+
!     &                (-4+(-1+rho0**2)**3*gb_xx_n(i,j,k))*
!     &                (-4+(-1+rho0**2)**3*gb_yy_n(i,j,k)))*
!     &                (-4+(-1+rho0**2)**3*gb_zz_n(i,j,k))**2)!/(-1+rho0**2)**8*y0**4
!              efe_tx_ires(i,j,k)=g0_ll(1,2)!Hads_l(2)+A_l(2)-boxx_l(2)
!              efe_ty_ires(i,j,k)=g0_ll(1,3)!Hads_l(3)+A_l(3)-boxx_l(3)
!              efe_xx_ires(i,j,k)=theta(i,j,k)
!              efe_xy_ires(i,j,k)=boxx_l(2)
!              efe_yy_ires(i,j,k)=-1/g0_uu(1,1)*(1-x(i)**2-y(j)**2)**2
!              efe_zz_ires(i,j,k)=0
!              do a=1,4
!                do b=1,4
!                  do c=1,4
!                    do d=1,4
!                      do i1=1,4
!                        do j1=1,4
!                          do k1=1,4
!                            do e=1,4
!                              efe_zz_ires(i,j,k)=efe_zz_ires(i,j,k)+
!     &                                          g0_ll(a,i1)*
!     &                                          g0_uu(b,j1)*
!     &                                          g0_uu(c,k1)*
!     &                                          g0_uu(d,e)*
!     &                                          riemann_ulll(a,b,c,d)*
!     &                                          riemann_ulll(i1,j1,k1,e)
!                            end do
!                          end do
!                        end do
!                      end do
!                    end do
!                  end do
!                end do
!              end do


                !--------------------------------------------------------------------------
                ! phi_res = phi,ab g^ab + phi,b g^ab,a + phi,c g^cb gamma^a_ab
                !         (= g^ab phi,ab - g^ab gamma^c_ab phi,c) 
                !--------------------------------------------------------------------------
!NO NEED TO RECOMPUTE DERIVATIVES OF FULL SCALAR FIELD PHI. THESE ARE ALREADY COMPUTED IN tensor_init
!                phi10=phi1_n(i,j,k)
!
!                phi_x(1)=(1-rho0**2)**2*phi10_x(1)
!                phi_x(2)=(1-rho0**2)**2*phi10_x(2)
!     &             -4*x0*(1-rho0**2)*phi10
!                phi_x(3)=(1-rho0**2)**2*phi10_x(3)
!     &             -4*y0*(1-rho0**2)*phi10
!                phi_x(2)=(1-rho0**2)**2*phi10_x(4)
!     &             -4*z0*(1-rho0**2)*phi10
!
!                phi_xx(1,1)=(1-rho0**2)**2*phi10_xx(1,1)
!                phi_xx(1,2)=(1-rho0**2)**2*phi10_xx(1,2)
!     &             -4*x0*(1-rho0**2)*phi10_x(1)
!                phi_xx(1,3)=(1-rho0**2)**2*phi10_xx(1,3)
!     &             -4*y0*(1-rho0**2)*phi10_x(1)
!                phi_xx(1,4)=(1-rho0**2)**2*phi10_xx(1,4)
!     &             -4*z0*(1-rho0**2)*phi10_x(1)
!                phi_xx(2,2)=(1-rho0**2)**2*phi10_xx(2,2)
!     &             -2*4*x0*(1-rho0**2)*phi10_x(2)
!     &             -4*(1-rho0**2-2*x0**2)*phi10
!                phi_xx(2,3)=(1-rho0**2)**2*phi10_xx(2,3)
!     &             -4*y0*(1-rho0**2)*phi10_x(2)           
!     &             -4*x0*(1-rho0**2)*phi10_x(3)
!     &             +8*x0*y0*phi10
!                phi_xx(2,4)=(1-rho0**2)**2*phi10_xx(2,4)
!     &             -4*z0*(1-rho0**2)*phi10_x(2)           
!     &             -4*x0*(1-rho0**2)*phi10_x(4)
!     &             +8*x0*z0*phi10
!                phi_xx(3,3)=(1-rho0**2)**2*phi10_xx(3,3)
!     &             -2*4*y0*(1-rho0**2)*phi10_x(3)
!     &             -4*(1-rho0**2-2*y0**2)*phi10
!                phi_xx(3,4)=(1-rho0**2)**2*phi10_xx(3,4)
!     &             -4*z0*(1-rho0**2)*phi10_x(3)           
!     &             -4*y0*(1-rho0**2)*phi10_x(4)
!     &             +8*y0*z0*phi10
!                phi_xx(4,4)=(1-rho0**2)**2*phi10_xx(4,4)
!     &             -2*4*z0*(1-rho0**2)*phi10_x(4)
!     &             -4*(1-rho0**2-2*z0**2)*phi10
!
!                do a=1,3
!                 do b=a+1,4
!                   phi_xx(b,a)=phi_xx(a,b)
!                 end do
!               end do

               kg_ires(i,j,k)=0.0d0
               do a=1,4
                do b=1,4
                 kg_ires(i,j,k)=
     &               kg_ires(i,j,k)
     &                +g0_uu(a,b)*phi10_xx(a,b)
                 do c=1,4
                   kg_ires(i,j,k)=
     &               kg_ires(i,j,k)
     &                -g0_uu(a,b)*gamma_ull(c,a,b)*phi10_x(c)
                 end do
                end do
               end do

            else
               efe_tt_ires(i,j,k)=0.0d0
               efe_tx_ires(i,j,k)=0.0d0
               efe_ty_ires(i,j,k)=0.0d0
               efe_tz_ires(i,j,k)=0.0d0
               efe_xx_ires(i,j,k)=0.0d0
               efe_xy_ires(i,j,k)=0.0d0
               efe_xz_ires(i,j,k)=0.0d0
               efe_yy_ires(i,j,k)=0.0d0
               efe_yz_ires(i,j,k)=0.0d0
               efe_zz_ires(i,j,k)=0.0d0
               kg_ires(i,j,k)=0.0d0
               ghconstr_t(i,j,k)=0.0d0
               ghconstr_x(i,j,k)=0.0d0
               ghconstr_y(i,j,k)=0.0d0
               ghconstr_z(i,j,k)=0.0d0
               ghconstr_all(i,j,k)=0.0d0
            end if



           end do
          end do
        end do

!!!!!DEBUGGING!!!!!!!
!        max_efe_all_ires=0.0d0
!        max_i=0
!        max_j=0
!        max_k=0
!
!        do i=is,ie
!          do j=js,je
!           do k=ks,ke
!
!            if (chr(i,j,k).ne.ex) then
!             if (efe_all_ires(i,j,k).gt.max_efe_all_ires) then
!              max_efe_all_ires=efe_all_ires(i,j,k)
!              max_i=i
!              max_j=j
!              max_k=k
!             end if
!            end if
!           end do
!          end do
!         end do
!
!        write(*,*) 'max_i,max_j,max_k,x(max_i),y(max_j),z(max_k),rho,
!     &              max_efe_all_ires='
!     &             ,max_i,max_j,max_k,x(max_i),y(max_j),z(max_k),
!     &              x(max_i)**2+y(max_j)**2+z(max_k)**2,
!     &              max_efe_all_ires 
!!!!!!!!!!!!!!!!!!!!!!!

        return
        end




c-----------------------------------------------------------------------
c calculate Kretschmann scalar and Riemann cube scalar
c-----------------------------------------------------------------------
        subroutine kretsch_riemanncube(kretsch_n,
     &                  kretschcentregrid,
     &                  riemanncube_n,
     &                  gb_tt_np1,gb_tt_n,gb_tt_nm1,
     &                  gb_tx_np1,gb_tx_n,gb_tx_nm1,
     &                  gb_ty_np1,gb_ty_n,gb_ty_nm1,
     &                  gb_tz_np1,gb_tz_n,gb_tz_nm1,
     &                  gb_xx_np1,gb_xx_n,gb_xx_nm1,
     &                  gb_xy_np1,gb_xy_n,gb_xy_nm1,
     &                  gb_xz_np1,gb_xz_n,gb_xz_nm1,
     &                  gb_yy_np1,gb_yy_n,gb_yy_nm1,
     &                  gb_yz_np1,gb_yz_n,gb_yz_nm1,
     &                  gb_zz_np1,gb_zz_n,gb_zz_nm1,
     &                  Hb_t_np1,Hb_t_n,Hb_t_nm1,
     &                  Hb_x_np1,Hb_x_n,Hb_x_nm1,
     &                  Hb_y_np1,Hb_y_n,Hb_y_nm1,
     &                  Hb_z_np1,Hb_z_n,Hb_z_nm1,
     &                  phi1_np1,phi1_n,phi1_nm1,
     &                  x,y,z,dt,chr,L,ex,Nx,Ny,Nz,phys_bdy,ghost_width,
     &                  ief_bh_r0,a_rot,kerrads_background,
     &                  output_kretsch,
     &                  output_riemanncube)

        implicit none
        real*8 ief_bh_r0,a_rot
        integer kerrads_background
        logical calc_der,calc_adv_quant
        data calc_der/.true./
        data calc_adv_quant/.false./
        integer output_kretsch
        integer output_riemanncube
        integer Nx,Ny,Nz
        integer i,j,k
        integer phys_bdy(6),ghost_width(6)
        real*8 chr(Nx,Ny,Nz),ex
        real*8 x(Nx),y(Ny),z(Nz),dt,L
        real*8 lambda4
        real*8 phi1_np1(Nx,Ny,Nz),phi1_n(Nx,Ny,Nz),phi1_nm1(Nx,Ny,Nz)
        real*8 gb_tt_np1(Nx,Ny,Nz),gb_tx_np1(Nx,Ny,Nz)
        real*8 gb_ty_np1(Nx,Ny,Nz)
        real*8 gb_tz_np1(Nx,Ny,Nz)
        real*8 gb_xx_np1(Nx,Ny,Nz),gb_xy_np1(Nx,Ny,Nz)
        real*8 gb_xz_np1(Nx,Ny,Nz)
        real*8 gb_yy_np1(Nx,Ny,Nz)
        real*8 gb_yz_np1(Nx,Ny,Nz)
        real*8 gb_zz_np1(Nx,Ny,Nz)
        real*8 gb_tt_n(Nx,Ny,Nz),gb_tx_n(Nx,Ny,Nz),gb_ty_n(Nx,Ny,Nz)
        real*8 gb_tz_n(Nx,Ny,Nz)
        real*8 gb_xx_n(Nx,Ny,Nz),gb_xy_n(Nx,Ny,Nz),gb_yy_n(Nx,Ny,Nz)
        real*8 gb_xz_n(Nx,Ny,Nz)
        real*8 gb_yz_n(Nx,Ny,Nz)
        real*8 gb_zz_n(Nx,Ny,Nz)
        real*8 gb_tt_nm1(Nx,Ny,Nz),gb_tx_nm1(Nx,Ny,Nz)
        real*8 gb_ty_nm1(Nx,Ny,Nz)
        real*8 gb_tz_nm1(Nx,Ny,Nz)
        real*8 gb_xx_nm1(Nx,Ny,Nz),gb_xy_nm1(Nx,Ny,Nz)
        real*8 gb_xz_nm1(Nx,Ny,Nz)
        real*8 gb_yy_nm1(Nx,Ny,Nz)
        real*8 gb_yz_nm1(Nx,Ny,Nz)
        real*8 gb_zz_nm1(Nx,Ny,Nz)

        real*8 Hb_t_np1(Nx,Ny,Nz),Hb_t_n(Nx,Ny,Nz),Hb_t_nm1(Nx,Ny,Nz)
        real*8 Hb_x_np1(Nx,Ny,Nz),Hb_x_n(Nx,Ny,Nz),Hb_x_nm1(Nx,Ny,Nz)
        real*8 Hb_y_np1(Nx,Ny,Nz),Hb_y_n(Nx,Ny,Nz),Hb_y_nm1(Nx,Ny,Nz)
        real*8 Hb_z_np1(Nx,Ny,Nz),Hb_z_n(Nx,Ny,Nz),Hb_z_nm1(Nx,Ny,Nz)

        integer is,ie,js,je,ks,ke

        integer i1,j1,k1,a,b,c,d,e,f,g,h,p,q,r
        integer ic,jc,kc
        real*8 efe_ires(4,4)

        real*8 PI
        parameter (PI=3.141592653589793d0)

        real*8 dx,dy,dz
        real*8 x0,y0,z0,rho0        

        real*8 boxx_u(4),boxx_l(4) 

        !--------------------------------------------------------------
        ! variables for tensor manipulations 
        !(indices are t,x,y,theta,phi)
        !--------------------------------------------------------------
        real*8 g0_ll(4,4),g0_uu(4,4)
        real*8 g0_ll_x(4,4,4),g0_uu_x(4,4,4),g0_ll_xx(4,4,4,4)
        real*8 gads_ll(4,4),gads_uu(4,4)
        real*8 gads_ll_x(4,4,4),gads_uu_x(4,4,4),gads_ll_xx(4,4,4,4)
        real*8 gammaads_ull(4,4,4)
        real*8 gammaads_ull_x(4,4,4,4)
        real*8 riemannads_ulll(4,4,4,4)
        real*8 riemannads_llll(4,4,4,4)
        real*8 riemannads_uuuu(4,4,4,4)
        real*8 riemannads_uull(4,4,4,4)
        real*8 kretschads,riemanncubeads

        real*8 h0_ll(4,4),h0_uu(4,4)
        real*8 h0_ll_x(4,4,4),h0_uu_x(4,4,4),h0_ll_xx(4,4,4,4)
        real*8 gamma_ull(4,4,4),gamma_ull_x(4,4,4,4)
        real*8 riemann_ulll(4,4,4,4)
        real*8 riemann_llll(4,4,4,4),riemann_uuuu(4,4,4,4)
        real*8 riemann_uull(4,4,4,4)
        real*8 ricci_ll(4,4),ricci_lu(4,4),ricci
        real*8 einstein_ll(4,4),set_ll(4,4)
        real*8 Hads_l(4),A_l(4),A_l_x(4,4)
        real*8 phi10_x(4),phi10_xx(4,4)

        !--------------------------------------------------------------
        ! variables for outward null expansion in spherical symmetry
        !--------------------------------------------------------------
        real*8 n_l(4),s_l(4)
        real*8 n_u(4),s_u(4)
        real*8 n_l_x(4,4),n_u_x(4,4),s_l_x(4,4)
        real*8 f0_x(4)
        real*8 f0_xx(4,4)
        real*8 gam_uu(4,4),sig_uu(4,4)
        real*8 gam_uu_x(4,4,4)
        real*8 normsusq

        real*8 g0gamfx(4)
        real*8 nufx,nuxfx(4),gamxfxfx(4)

        real*8 theta(Nx,Ny,Nz)

        real*8 kretsch_n(Nx,Ny,Nz)
        real*8 kretschcentregrid
        real*8 riemanncube_n(Nx,Ny,Nz)

!!!!!!!DEBUGGING!!!!!!
        integer max_i,max_j,max_k
        real*8  max_efe_all_ires
!!!!!!!!!!!!!!!!!!!!!1

        ! initialize fixed-size variables
        data i,j,k,is,ie,js,je,ks,ke/0,0,0,0,0,0,0,0,0/
        data ic,jc,kc/0,0,0/
        data i1,j1,k1,a,b,c,d,e,p,q,r/0,0,0,0,0,0,0,0,0,0,0/

        data dx,dy,dz/0.0,0.0,0.0/
        data x0,y0,rho0/0.0,0.0,0.0/    

        data g0_ll,g0_uu/16*0.0,16*0.0/
        data gads_ll,gads_uu/16*0.0,16*0.0/
        data h0_ll,h0_uu/16*0.0,16*0.0/
        data gamma_ull/64*0.0/
        data gamma_ull_x/256*0.0/

        data g0_ll_x,g0_uu_x/64*0.0,64*0.0/
        data gads_ll_x,gads_uu_x/64*0.0,64*0.0/
        data h0_ll_x,h0_uu_x/64*0.0,64*0.0/

        data g0_ll_xx/256*0.0/
        data gads_ll_xx/256*0.0/
        data h0_ll_xx/256*0.0/

        data ricci/0.0/
        data ricci_ll,ricci_lu/16*0.0,16*0.0/
        data einstein_ll,set_ll/16*0.0,16*0.0/
        data riemann_ulll/256*0.0/
        data riemann_llll/256*0.0/
        data riemann_uuuu/256*0.0/
        data riemann_uull/256*0.0/

        data A_l,Hads_l/4*0.0,4*0.0/
        data A_l_x/16*0.0/

        data phi10_x/4*0.0/
        data phi10_xx/16*0.0/

        data boxx_u,boxx_l/4*0.0,4*0.0/

!----------------------------------------------------------------------


        dx=(x(2)-x(1))
        dy=(y(2)-y(1))
        dz=(z(2)-z(1))

        ! set index bounds for main loop
        is=2
        ie=Nx-1
        js=2
        je=Ny-1
        ks=2
        ke=Nz-1

        ! adjust index bounds to compensate for ghost_width
        if (ghost_width(1).gt.0) is=is+ghost_width(1)-1
        if (ghost_width(2).gt.0) ie=ie-(ghost_width(2)-1)
        if (ghost_width(3).gt.0) js=js+ghost_width(3)-1
        if (ghost_width(4).gt.0) je=je-(ghost_width(4)-1)
        if (ghost_width(5).gt.0) ks=ks+ghost_width(5)-1
        if (ghost_width(6).gt.0) ke=ke-(ghost_width(6)-1)


        ! (MAIN LOOP) loop through spacetime points x(i),y(j)
        do i=is,ie
          do j=js,je
           do k=ks,ke

              x0=x(i)
              y0=y(j)
              z0=z(k)
              rho0=sqrt(x0**2+y0**2+z0**2)

            if (chr(i,j,k).ne.ex) then

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
                  do d=1,4
                   riemann_llll(a,b,c,d)=0.0d0
                   riemann_uuuu(a,b,c,d)=0.0d0
                   riemann_uull(a,b,c,d)=0.0d0
                   do e=1,4
                     riemann_llll(a,b,c,d)=
     &                     riemann_llll(a,b,c,d)
     &                     +riemann_ulll(e,b,c,d)*g0_ll(e,a)
                     if (output_riemanncube.eq.1) then
                      riemann_uull(a,b,c,d)=riemann_uull(a,b,c,d)
     &                             +riemann_ulll(a,e,c,d)*g0_uu(e,b)
                     end if
                    do f=1,4
                     do g=1,4
                     riemann_uuuu(a,b,c,d)=
     &                          riemann_uuuu(a,b,c,d)
     &                          +riemann_ulll(a,e,f,g)*g0_uu(e,b)
     &                            *g0_uu(f,c)*g0_uu(g,d)
                     end do
                    end do
                   end do

                  end do
                 end do
                end do
               end do


               kretsch_n(i,j,k)=0.0d0
               riemanncube_n(i,j,k)=0.0d0
               if ((output_kretsch.eq.1)
     &          .or.(output_riemanncube.eq.1)) then
               do a=1,4
                do b=1,4
                 do c=1,4
                  do d=1,4
                    if (output_kretsch.eq.1) then
                       kretsch_n(i,j,k)=
     &                         kretsch_n(i,j,k)
     &                         +riemann_uuuu(a,b,c,d)
     &                         *riemann_llll(a,b,c,d)
                    end if

                    if (output_riemanncube.eq.1) then
                     do e=1,4
                      do f=1,4
                       riemanncube_n(i,j,k)=
     &                         riemanncube_n(i,j,k)
     &                         +riemann_uuuu(a,b,c,d)
     &                         *riemann_llll(c,d,e,f)
     &                         *riemann_uull(e,f,a,b)

                      end do
                     end do
                    end if
       
                  end do
                 end do
                end do
               end do
               end if


            else
               kretsch_n(i,j,k)=0.0d0
               riemanncube_n(i,j,k)=0.0d0
            end if

!find the indices denoting the point at the centre of the grid. Needed to compute kretschcentregrid
               if (     (abs(x0).lt.10.0d0**(-10))
     &             .and.(abs(y0).lt.10.0d0**(-10))
     &             .and.(abs(z0).lt.10.0d0**(-10)) ) then

                  ic=i
                  jc=j
                  kc=k
               end if


           end do
          end do
        end do

        kretschcentregrid=0.0d0

        if ((ic.gt.0).and.(jc.gt.0).and.(kc.gt.0)) then !this condition is activated only if the processor calling ires contains the centre of the grid (where ic,jc and kc are set to a positive number by the cycle above)
             kretschcentregrid=kretsch_n(ic,jc,kc)
        else
             kretschcentregrid=0.0d0
        end if

!           write(*,*) "ex,chr(ic,jc,kc)",ex,chr(ic,jc,kc)
!           write(*,*) "ic,jc,kc=",ic,jc,kc
!           write(*,*) "kretschcentregrid=",kretschcentregrid


        return
        end

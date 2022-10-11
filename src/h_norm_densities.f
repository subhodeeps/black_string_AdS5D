c-----------------------------------------------------------------------
c Calculate square root of the integrand of the H_AdS^(0,sp), H_AdS^(1,sp), H_AdS^(2,sp) norms 
c (defined in arxiv:1103.0710v2 and used in arxiv:1110.6794v2) of the scalar field phi.
c We replace the coordinates of arxiv:1110.6794v2 with uncompactified quasi-spherical Kerr-Schild coords.
c We write the (Euclidean) volume measure dV=R^2sin(theta)dR dtheta dphi as 
c R^2/rho^2dR_drho rho^2 sin(theta)drho dtheta dphi=R^2/rho^2 dR_drho dxdydz.
c In the definition of sqrth0/1/2spnormdensity_f we only include the factors R^2/rho^2 dR_drho. Therefore, in order
c to obtain the H_AdS^(0/1/2,sp) norm for quasi-spherical coordinates, i.e., ||f||^2_{H_AdS^(0/1/2,sp)},
c sqrth0/1/2spnormdensity_f needs to be squared and integrated in dxdydz.
c This can be done (up to a numerical factor coming from the integral over the entire space), for example, 
c by using the L^2-norm calculation in DV and then squaring the result.
c We do not need to square the L^2-norm if we need ||f||_{H_AdS^(0/1/2,sp)} (instead of ||f||^2_{H_AdS^(0/1/2,sp)}).
c The variable f1 represents the evoluton variable \barphi. 
c The variable f represents, alternatively, the functions whose norm we want to compute. 
c The exact function that f represents is determined by the value of hnorm_argtype:
c - hnorm_argtype=0 corresponds to f=(1-rho^2)^2*f1,
c - hnorm_argtype=1 corresponds to f=(1-rho^2)^2*f1_t,
c - hnorm_argtype=2 corresponds to f=(1-rho^2)^2*f1_tt,
c - hnorm_argtype=3 corresponds to f=(1-rho^2)^2*m_1[f1], where m_1=d/dphi is 
c the first generator of the asymptotic rotation symmetry group,
c - hnorm_argtype=4 corresponds to f=(1-rho^2)^2*m_2[f1], where m_2=-sin(phi)d/dtheta-cot(theta)cos(phi)d/dphi is 
c the second generator of the asymptotic rotation symmetry group,
c - hnorm_argtype=5 corresponds to f=(1-rho^2)^2*m_3[f1], where m_3=cos(phi)d/dtheta-cot(theta)sin(phi)d/dphi is 
c the third generator of the asymptotic rotation symmetry group,
c-----------------------------------------------------------------------


c-----------------------------------------------------------------------
c Compute square root of H_AdS^(0,sp) norm density of f
c-----------------------------------------------------------------------
        subroutine sqrth0spnormdensity_func(
     &                  sqrth0spnormdensity_f,
     &                  sp,hnorm_argtype,
     &                  f1_np1,f1_n,f1_nm1,
     &                  AH_R,AH_xc,
     &                  AH_Ntheta,AH_Nphi,
     &                  x,y,z,dt,ct,chr,L,ex,Nx,Ny,Nz,
     &                  phys_bdy,ghost_width,
     &                  ief_bh_r0,a_rot,kerrads_background)

        implicit none

        logical is_nan
        real*8  ief_bh_r0,a_rot,M0,M0_min
        real*8  kerrads_background
        integer Nx,Ny,Nz
        integer i,j,k
        integer phys_bdy(6),ghost_width(6)
        real*8 chr(Nx,Ny,Nz),ex
        real*8 x(Nx),y(Ny),z(Nz),dt,ct,L
        real*8 lambda4
        integer AH_Ntheta,AH_Nphi
        real*8 AH_R(AH_Ntheta,AH_Nphi)
        real*8 AH_xc(3)
        real*8 xp,yp,zp,rhop,thetap,phip
        integer i0,j0
        real*8 dahtheta,dahphi
        real*8 ft,fp
        real*8 AH_Rp,rhostop
        real*8 f1_np1(Nx,Ny,Nz),f1_n(Nx,Ny,Nz),f1_nm1(Nx,Ny,Nz)
        real*8 sqrth0spnormdensity_f(Nx,Ny,Nz)
        real*8 h0spnormdensity_f0
        integer sp,hnorm_argtype
        real*8 f10,f0
        real*8 df_drho,df_dtheta,df_dphi
        real*8 f1_t,f1_x,f1_y,f1_z
        real*8 f_t,f_x,f_y,f_z
        real*8 f1_tt,f1_tx,f1_ty,f1_tz
        real*8 f1_xx,f1_xy,f1_xz
        real*8 f1_yy,f1_yz,f1_zz
        real*8 f1_rho,f1_theta,f1_phi
        real*8 derf_dRad

        real*8 rblhor,Radhor,rhohor

        logical ltrace
        parameter (ltrace=.false.)


        integer is,ie,js,je,ks,ke

        integer i1,j1,k1,a,b,c,d,e,f,g,h,p,q,r
        integer ic,jc,kc
        real*8 efe_ires(4,4)

        real*8 PI
        parameter (PI=3.141592653589793d0)

        real*8 dx,dy,dz
        real*8 x0,y0,z0
        real*8 rho0,theta0,phi0   
        real*8 Rad0
        real*8 drho_dRad,dRad_drho



        ! initialize fixed-size variables
        data i,j,k,is,ie,js,je,ks,ke/0,0,0,0,0,0,0,0,0/
        data ic,jc,kc/0,0,0/
        data i1,j1,k1,a,b,c,d,e,p,q,r/0,0,0,0,0,0,0,0,0,0,0/

        data dx,dy,dz/0.0,0.0,0.0/
        data x0,y0,z0,rho0/0.0,0.0,0.0,0.0/    



!----------------------------------------------------------------------


      ! Black hole mass
        M0=ief_bh_r0/2
      ! Minimum black hole mass. For M0 below this value, there is a naked singularity
        M0_min=((2*(1 + a_rot**2/L**2) + Sqrt((1 + a_rot**2/L**2)**2 
     &   + (12*a_rot**2)/L**2))*Sqrt(-1 + Sqrt((1 + a_rot**2/L**2)**2 
     &   + (12*a_rot**2)/L**2) - a_rot**2/L**2))/(3.*Sqrt(6.))

        if (a_rot.ge.L) then
         write (*,*) "sqrth0spnormdensity_phi: ERROR in choice of Kerr-AdS initial parameters: 
     &   the rotation parameter a must be smaller than the AdS radius L"
          write (*,*) "a_rot,L=",a_rot,L
          stop
        end if

        if ((abs(M0).gt.10.0d0**(-10))
     &     .and.(M0.le.M0_min)) then
          write (*,*) "sqrth0spnormdensity_phi: ERROR in choice of Kerr-AdS initial parameters: 
     &      the black hole mass M0=2*r0 must be larger
     &      than the M0_min value"
          write (*,*) "M0,M0_min=",M0,M0_min
          stop
        end if


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


        ! (MAIN LOOP) loop through spacetime points x(i),y(j),z(k)
        do i=is,ie
          do j=js,je
           do k=ks,ke

              x0=x(i)
              y0=y(j)
              z0=z(k)
              rho0=sqrt(x0**2+y0**2+z0**2)
              theta0=acos(x0/rho0)
              if (z0.lt.0) then
               phi0=atan2(z0,y0)+2*PI
              else
               phi0=atan2(z0,y0)
              end if

              !uncompactified quasi-spherical Kerr-Schild radial coordinate
              Rad0=2*rho0/(1-rho0**2)
              drho_dRad=((1-rho0**2)**2)/(2*(1+rho0**2))
              dRad_drho=(-1 + rho0)**(-2) 
     &           + (1 + rho0)**(-2)

              !determine radius at which we stop computing the norms. This is the apparent horizon (AH) radius.

              !Cartesian coordinates of grid point p=(i,j,k) w.r.t. the centre of the AH
              xp=x0-AH_xc(1)
              yp=y0-AH_xc(2)
              zp=z0-AH_xc(3)

              !The center of the AH xp=yp=zp=0 will be excluded or included depending on the excision mask chr(i,j,k), 
              !so no need to force inclusion or exclusion via the value of rhostop. 
              !In fact, xp=yp=zp=0 is also the centre of the grid if there is no AH (also if the AH is centred at the origin), 
              !so we should include it if we want to compute these norms in the case of no horizon.
              !For this reason, it's best to leave the inclusion or exclusion of that point to the excision mask chr(i,j,k).
              if ((abs(xp).lt.10.0d0**(-12)).and.
     &            (abs(yp).lt.10.0d0**(-12)).and.  
     &            (abs(zp).lt.10.0d0**(-12)) ) then
                rhostop=-1
              else
                !corresponding spherical coordinates of point p
                !thetap,phip is in general NOT on the AH_Ntheta,AH_Nphi grid.
                rhop=sqrt(xp**2+yp**2+zp**2)
                thetap=acos(xp/rhop)
                if (zp.lt.0) then
                    phip=atan2(zp,yp)+2*PI
                else
                    phip=atan2(zp,yp)
                end if

                ! i0 is the index of the point with chi value closest to chip on the AH_Nchi grid (indices going from 0 to AH_Ntheta-1),
                ! and larger than chip (except if thetap=Pi, then it is smaller)
                ! similarly for j0
                dahtheta=PI/(AH_Ntheta-1)
                dahphi=2*PI/(AH_Nphi-1)
                i0=floor(thetap/dahtheta+1)
                if (AH_Ntheta-1<i0) i0=AH_Ntheta-1
                j0=floor(phip/dahphi+1)
                if (AH_Nphi-1<j0) j0=AH_Nphi-1
                ft=((thetap-(i0-1)*dahtheta))/dahtheta
                fp=((phip-(j0-1)*dahphi))/dahphi

                !bi-linear interpolated value of AH radius at (thetap,phip)
                AH_Rp=  AH_R(i0,j0)*(1-ft)*(1-fp)+
     &                  AH_R(i0+1,j0)*(ft)*(1-fp)+
     &                  AH_R(i0,j0+1)*(1-ft)*(fp)+
     &                  AH_R(i0+1,j0+1)*(ft)*(fp)

                if (abs(AH_Rp-1).gt.10.0d0**(-10)) then 
                    rhostop=AH_Rp
                else
                    rhostop=-1.0d0
                end if
              end if


            if (
     &          (rho0.gt.rhostop).and.
     &          (chr(i,j,k).ne.ex)
     &          ) then

                ! set f1 value
                f10=f1_n(i,j,k)

                !calculate R^s*f^2*(R/rho)^2*dR_drho
                !To be integrated in rho^2 sin(theta)drho dtheta dphi, or dx dy dz

                if (hnorm_argtype.eq.0) then 
                !compute norm for f=(1-rho^2)^2*f1

                 f0=(1-rho0**2)**2*f10

                else if (hnorm_argtype.eq.1) then
                !compute norm for f=(1-rho^2)^2*f1_t

                 call   df1_int(f1_np1,f1_n,f1_nm1,
     &                  f1_t,f1_x,f1_y,f1_z,
     &                  x,y,z,dt,i,j,k,chr,ex,Nx,Ny,Nz,'f1')

                 f0=(1-rho0**2)**2*f1_t

                else if (hnorm_argtype.eq.2) then
                !compute norm for f=(1-rho^2)^2*f1_tt

                 call df2_int(f1_np1,f1_n,f1_nm1,
     &                      f1_t,f1_x,f1_y,f1_z,
     &                      f1_tt,f1_tx,f1_ty,f1_tz,
     &                      f1_xx,f1_xy,f1_xz,
     &                      f1_yy,f1_yz,f1_zz,
     &                      x,y,z,dt,i,j,k,chr,ex,Nx,Ny,Nz,'f1')

                 f0=(1-rho0**2)**2*f1_tt

                else if (hnorm_argtype.eq.3) then
                !compute norm for f=m_1[(1-rho^2)^2*f1]=(1-rho^2)^2*m1[f1], 
                !where m_1=d/dphi is the first generator of the asymptotic rotation symmetry group

                 call derf_dxsph(f1_np1,f1_n,f1_nm1,
     &              f1_t,f1_rho,f1_theta,f1_phi,
     &              x,y,z,dt,i,j,k,chr,ex,Nx,Ny,Nz,'f1')
                 f0=(1-rho0**2)**2*f1_phi

                else if (hnorm_argtype.eq.4) then
                !compute norm for f=m_2[(1-rho^2)^2*f1]=(1-rho^2)^2*m2[f1], 
                !where m_2=-sin(phi)d/dtheta-cot(theta)cos(phi)d/dphi is 
                !the second generator of the asymptotic rotation symmetry group

                 call derf_dxsph(f1_np1,f1_n,f1_nm1,
     &              f1_t,f1_rho,f1_theta,f1_phi,
     &              x,y,z,dt,i,j,k,chr,ex,Nx,Ny,Nz,'f1')
                 f0=(1-rho0**2)**2*
     &                  (-sin(phi0)*f1_theta
     &           -(cos(theta0)/sin(theta0))*cos(phi0)*f1_phi)


                else if (hnorm_argtype.eq.5) then
                !compute norm for f=m_3[(1-rho^2)^2*f1]=(1-rho^2)^2*m3[f1], 
                !where m_3=cos(phi)d/dtheta-cot(theta)sin(phi)d/dphi is 
                !the third generator of the asymptotic rotation symmetry group

                 call derf_dxsph(f1_np1,f1_n,f1_nm1,
     &              f1_t,f1_rho,f1_theta,f1_phi,
     &              x,y,z,dt,i,j,k,chr,ex,Nx,Ny,Nz,'f1')

                 f0=(1-rho0**2)**2*
     &                  (cos(phi0)*f1_theta
     &           -(cos(theta0)/sin(theta0))*sin(phi0)*f1_phi)

                else
                 write(*,*) " WARNING - sqrth0spnormdensity_func:"//
     &  " invalid value of hnorm_argtype: hnorm_argtype=",hnorm_argtype
                 write(*,*) " hnorm_argtype can only be 0,1,2,3,4,5"
                 return
                end if

                h0spnormdensity_f0=f0**2

                h0spnormdensity_f0=
     &          (Rad0**sp)*h0spnormdensity_f0
     &          *Rad0**2/rho0**2
     &          *dRad_drho



!Take square root, so we can use DV calculation of L2-norm, and then square it, to compute the H_AdS(1,s) norm.
                sqrth0spnormdensity_f(i,j,k)=
     &             sqrt(h0spnormdensity_f0)

                 !the phi coordinate is not defined at y=z=0, i.e., theta=0,PI 
                !(not surprising since spherical coordinates are not valid everywhere on the sphere), 
                !hence we set the norm density to 0 at these points

                if ((abs(y0).lt.10.0d0**(-10)).and.
     &          (abs(z0).lt.10.0d0**(-10))) then
                    sqrth0spnormdensity_f(i,j,k)=0.0d0
                end if

                if (ltrace) then
!                 if ((abs(x0).gt.10.0d0**(-1.0d0)).and.
!     &               (abs(y0).gt.10.0d0**(-1.0d0)).and.
!     &               (abs(z0).gt.10.0d0**(-1.0d0)).and.
!     &                (rho0.lt.0.9)  ) then
                   write(*,*) 'sqrth0spnormdensity_func '
                   write(*,*) 'sp,hnorm_argtype= ',sp,hnorm_argtype
                   write(*,*) 'at i,j,k= ', i,j,k
                   write(*,*) 'i.e., x,y,z,rho=', x0,y0,z0,rho0
                   write(*,*) ' sqrth0spnormdensity_f = ',
     &                  sqrth0spnormdensity_f(i,j,k)
                   return
!                 end if
                end if


            else !i.e., the point is either excised or inside the apparent horizon
               sqrth0spnormdensity_f(i,j,k)=0.0d0
            end if

           end do
          end do
        end do


        return
        end



c-----------------------------------------------------------------------
c Compute square root of H_AdS^(1,sp) norm density of f
c-----------------------------------------------------------------------
        subroutine sqrth1spnormdensity_func(
     &                  sqrth1spnormdensity_f,
     &                  sp,hnorm_argtype,
     &                  f1_np1,f1_n,f1_nm1,
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
     &                  AH_R,AH_xc,
     &                  AH_Ntheta,AH_Nphi,
     &                  x,y,z,dt,ct,chr,L,ex,Nx,Ny,Nz,
     &                  phys_bdy,ghost_width,
     &                  ief_bh_r0,a_rot,kerrads_background)

        implicit none

        logical is_nan

        real*8  ief_bh_r0,a_rot,M0,M0_min
        real*8  kerrads_background

        logical calc_der,calc_adv_quant
        data calc_der/.false./
        data calc_adv_quant/.false./

        integer Nx,Ny,Nz
        integer i,j,k
        integer phys_bdy(6),ghost_width(6)
        real*8 chr(Nx,Ny,Nz),ex
        real*8 x(Nx),y(Ny),z(Nz),dt,ct,L
        real*8 lambda4
        integer AH_Ntheta,AH_Nphi
        real*8 AH_R(AH_Ntheta,AH_Nphi)
        real*8 AH_xc(3)
        real*8 xp,yp,zp,rhop,thetap,phip
        integer i0,j0
        real*8 dahtheta,dahphi
        real*8 ft,fp
        real*8 AH_Rp,rhostop
        real*8 f1_np1(Nx,Ny,Nz),f1_n(Nx,Ny,Nz),f1_nm1(Nx,Ny,Nz)
        real*8 sqrth1spnormdensity_f(Nx,Ny,Nz)
        real*8 h1spnormdensity_f0
        integer sp,hnorm_argtype
        real*8 f10,f0
        real*8 df1_dxqssph(4),d2f1_dxqssphdxqssph(4,4)
        real*8 f1_t,f1_x,f1_y,f1_z
        real*8 f1_rho,f1_theta,f1_phi
        real*8 df_dRad,df_dtheta,df_dphi

        logical ltrace
        parameter (ltrace=.false.)


        integer is,ie,js,je,ks,ke

        integer i1,j1,k1,a,b,c,d,e,f,g,h,p,q,r
        integer ic,jc,kc
        real*8 efe_ires(4,4)

        real*8 PI
        parameter (PI=3.141592653589793d0)

        real*8 dx,dy,dz
        real*8 x0,y0,z0
        real*8 rho0,theta0,phi0   
        real*8 Rad0
        real*8 drho_dRad,dRad_drho

        real*8 dxcar_dxqssph(4,4)

        real*8 gb_tt_np1(Nx,Ny,Nz),gb_tt_n(Nx,Ny,Nz),gb_tt_nm1(Nx,Ny,Nz)
        real*8 gb_tx_np1(Nx,Ny,Nz),gb_tx_n(Nx,Ny,Nz),gb_tx_nm1(Nx,Ny,Nz)
        real*8 gb_ty_np1(Nx,Ny,Nz),gb_ty_n(Nx,Ny,Nz),gb_ty_nm1(Nx,Ny,Nz)
        real*8 gb_tz_np1(Nx,Ny,Nz),gb_tz_n(Nx,Ny,Nz),gb_tz_nm1(Nx,Ny,Nz)
        real*8 gb_xx_np1(Nx,Ny,Nz),gb_xx_n(Nx,Ny,Nz),gb_xx_nm1(Nx,Ny,Nz)
        real*8 gb_xy_np1(Nx,Ny,Nz),gb_xy_n(Nx,Ny,Nz),gb_xy_nm1(Nx,Ny,Nz)
        real*8 gb_xz_np1(Nx,Ny,Nz),gb_xz_n(Nx,Ny,Nz),gb_xz_nm1(Nx,Ny,Nz)
        real*8 gb_yy_np1(Nx,Ny,Nz),gb_yy_n(Nx,Ny,Nz),gb_yy_nm1(Nx,Ny,Nz)
        real*8 gb_yz_np1(Nx,Ny,Nz),gb_yz_n(Nx,Ny,Nz),gb_yz_nm1(Nx,Ny,Nz)
        real*8 gb_zz_np1(Nx,Ny,Nz),gb_zz_n(Nx,Ny,Nz),gb_zz_nm1(Nx,Ny,Nz)
        real*8 Hb_t_np1(Nx,Ny,Nz),Hb_t_n(Nx,Ny,Nz),Hb_t_nm1(Nx,Ny,Nz)
        real*8 Hb_x_np1(Nx,Ny,Nz),Hb_x_n(Nx,Ny,Nz),Hb_x_nm1(Nx,Ny,Nz)
        real*8 Hb_y_np1(Nx,Ny,Nz),Hb_y_n(Nx,Ny,Nz),Hb_y_nm1(Nx,Ny,Nz)
        real*8 Hb_z_np1(Nx,Ny,Nz),Hb_z_n(Nx,Ny,Nz),Hb_z_nm1(Nx,Ny,Nz)
        real*8 phi1_np1(Nx,Ny,Nz),phi1_n(Nx,Ny,Nz),phi1_nm1(Nx,Ny,Nz)

        !--------------------------------------------------------------
        ! variables for tensor manipulations 
        !(indices are t,x,y,theta,phi)
        !--------------------------------------------------------------
        real*8 gads_ll(4,4),gads_uu(4,4)
        real*8 gads_ll_x(4,4,4),gads_uu_x(4,4,4)
        real*8 gads_ll_xx(4,4,4,4)
        real*8 Hads_l(4)
        real*8 phi1ads, phi1ads_x(4)

        real*8 g0_ll(4,4),g0_uu(4,4),detg0
        real*8 g0_ll_x(4,4,4),g0_uu_x(4,4,4),g0_ll_xx(4,4,4,4)
        real*8 h0_ll(4,4),h0_uu(4,4)
        real*8 h0_ll_x(4,4,4),h0_uu_x(4,4,4),h0_ll_xx(4,4,4,4)
        real*8 gamma_ull(4,4,4),gamma_ull_x(4,4,4,4)
        real*8 riemann_ulll(4,4,4,4)
        real*8 ricci_ll(4,4),ricci_lu(4,4),ricci
        real*8 einstein_ll(4,4),set_ll(4,4)
        real*8 A_l(4),A_l_x(4,4)
        real*8 phi10_x(4),phi10_xx(4,4)

        real*8 detg0_qssph0
        real*8 g0_ll_qssph(4,4),g0_uu_qssph(4,4)



        ! initialize fixed-size variables
        data i,j,k,is,ie,js,je,ks,ke/0,0,0,0,0,0,0,0,0/
        data ic,jc,kc/0,0,0/
        data i1,j1,k1,a,b,c,d,e,p,q,r/0,0,0,0,0,0,0,0,0,0,0/

        data dx,dy,dz/0.0,0.0,0.0/
        data x0,y0,z0,rho0/0.0,0.0,0.0,0.0/    

        data gads_ll,gads_uu/16*0.0,16*0.0/
        data gads_ll_x,gads_uu_x/64*0.0,64*0.0/
        data gads_ll_xx/256*0.0/
        data Hads_l/4*0.0/
        data phi1ads_x/4*0.0/

        data g0_ll,g0_uu/16*0.0,16*0.0/
        data g0_ll_x,g0_uu_x/64*0.0,64*0.0/
        data g0_ll_xx/256*0.0/
        data h0_ll,h0_uu/16*0.0,16*0.0/
        data h0_ll_x,h0_uu_x/64*0.0,64*0.0/
        data h0_ll_xx/256*0.0/
        data gamma_ull/64*0.0/
        data gamma_ull_x/256*0.0/
        data riemann_ulll/256*0.0/
        data ricci_ll,ricci_lu/16*0.0,16*0.0/
        data einstein_ll,set_ll/16*0.0,16*0.0/
        data A_l,A_l_x/4*0.0,16*0.0/
        data phi10_x,phi10_xx/4*0.0,16*0.0/
        data dxcar_dxqssph/16*0.0/

        data g0_ll_qssph/16*0.0/
        data g0_uu_qssph/16*0.0/



!----------------------------------------------------------------------


      ! Black hole mass
        M0=ief_bh_r0/2
      ! Minimum black hole mass. For M0 below this value, there is a naked singularity
        M0_min=((2*(1 + a_rot**2/L**2) + Sqrt((1 + a_rot**2/L**2)**2 
     &   + (12*a_rot**2)/L**2))*Sqrt(-1 + Sqrt((1 + a_rot**2/L**2)**2 
     &   + (12*a_rot**2)/L**2) - a_rot**2/L**2))/(3.*Sqrt(6.))

        if (a_rot.ge.L) then
         write (*,*) "sqrth1spnormdensity_phi: ERROR in choice of Kerr-AdS initial parameters: 
     &   the rotation parameter a must be smaller than the AdS radius L"
          write (*,*) "a_rot,L=",a_rot,L
          stop
        end if

        if ((abs(M0).gt.10.0d0**(-10))
     &     .and.(M0.le.M0_min)) then
          write (*,*) "sqrth1spnormdensity_phi: ERROR in choice of Kerr-AdS initial parameters: 
     &      the black hole mass M0=2*r0 must be larger
     &      than the M0_min value"
          write (*,*) "M0,M0_min=",M0,M0_min
          stop
        end if


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
              theta0=acos(x0/rho0)
              if (z0.lt.0) then
               phi0=atan2(z0,y0)+2*PI
              else
               phi0=atan2(z0,y0)
              end if
              !uncompactified quasi-spherical Kerr-Schild radial coordinate
              Rad0=2*rho0/(1-rho0**2)
              drho_dRad=((1-rho0**2)**2)/(2*(1+rho0**2))
              dRad_drho=(-1 + rho0)**(-2) 
     &           + (1 + rho0)**(-2)

              !determine radius at which we stop computing the norms. This is the apparent horizon (AH) radius.

              !Cartesian coordinates of grid point p=(i,j,k) w.r.t. the centre of the AH
              xp=x0-AH_xc(1)
              yp=y0-AH_xc(2)
              zp=z0-AH_xc(3)

              !The center of the AH xp=yp=zp=0 will be excluded or included depending on the excision mask chr(i,j,k), 
              !so no need to force inclusion or exclusion via the value of rhostop. 
              !In fact, xp=yp=zp=0 is also the centre of the grid if there is no AH (also if the AH is centred at the origin), 
              !so we should include it if we want to compute these norms in the case of no horizon.
              !For this reason, it's best to leave the inclusion or exclusion of that point to the excision mask chr(i,j,k).
              if ((abs(xp).lt.10.0d0**(-12)).and.
     &            (abs(yp).lt.10.0d0**(-12)).and.  
     &            (abs(zp).lt.10.0d0**(-12)) ) then
                rhostop=-1
              else
                !corresponding spherical coordinates of point p
                !thetap,phip is in general NOT on the AH_Ntheta,AH_Nphi grid.
                rhop=sqrt(xp**2+yp**2+zp**2)
                thetap=acos(xp/rhop)
                if (zp.lt.0) then
                    phip=atan2(zp,yp)+2*PI
                else
                    phip=atan2(zp,yp)
                end if

                ! i0 is the index of the point with chi value closest to chip on the AH_Nchi grid (indices going from 0 to AH_Ntheta-1),
                ! and larger than chip (except if thetap=Pi, then it is smaller)
                ! similarly for j0
                dahtheta=PI/(AH_Ntheta-1)
                dahphi=2*PI/(AH_Nphi-1)
                i0=floor(thetap/dahtheta+1)
                if (AH_Ntheta-1<i0) i0=AH_Ntheta-1
                j0=floor(phip/dahphi+1)
                if (AH_Nphi-1<j0) j0=AH_Nphi-1
                ft=((thetap-(i0-1)*dahtheta))/dahtheta
                fp=((phip-(j0-1)*dahphi))/dahphi

                !bi-linear interpolated value of AH radius at (thetap,phip)
                AH_Rp=  AH_R(i0,j0)*(1-ft)*(1-fp)+
     &                  AH_R(i0+1,j0)*(ft)*(1-fp)+
     &                  AH_R(i0,j0+1)*(1-ft)*(fp)+
     &                  AH_R(i0+1,j0+1)*(ft)*(fp)

                if (AH_Rp.lt.1) then 
                    rhostop=AH_Rp
                else
                    rhostop=-1
                end if
              end if



            if ((rho0.gt.rhostop).and.
     &         (chr(i,j,k).ne.ex)) then

                ! set f1 value
                f10=f1_n(i,j,k)

    !compute metric and its derivatives 
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
     &                gads_ll,gads_uu,
     &                gads_ll_x,gads_uu_x,gads_ll_xx,
     &                h0_ll,h0_uu,h0_ll_x,h0_uu_x,h0_ll_xx,
     &                A_l,A_l_x,Hads_l,
     &                gamma_ull,gamma_ull_x,
     &                riemann_ulll,ricci_ll,ricci_lu,ricci,
     &                einstein_ll,set_ll,
     &                phi10_x,phi10_xx,
     &                x,y,z,dt,chr,L,ex,Nx,Ny,Nz,i,j,k,
     &                ief_bh_r0,a_rot,kerrads_background,
     &                calc_der,calc_adv_quant)



!define transformation matrix between Cartesian coordinates to compactified (quasi-)spherical coordinates, 
!        !e.g. dxcar_dxqssph(3,2)=dy/drho=sin(theta0)*cos(phi0)

             dxcar_dxqssph(1,1)=1
             dxcar_dxqssph(1,2)=0
             dxcar_dxqssph(1,3)=0
             dxcar_dxqssph(1,4)=0
             dxcar_dxqssph(2,1)=0
             dxcar_dxqssph(2,2)=cos(theta0)
             dxcar_dxqssph(2,3)=-rho0*sin(theta0)
             dxcar_dxqssph(2,4)=0
             dxcar_dxqssph(3,1)=0
             dxcar_dxqssph(3,2)=sin(theta0)*cos(phi0)
             dxcar_dxqssph(3,3)=rho0*cos(theta0)*cos(phi0)
             dxcar_dxqssph(3,4)=-rho0*sin(theta0)*sin(phi0)
             dxcar_dxqssph(4,1)=0
             dxcar_dxqssph(4,2)=sin(theta0)*sin(phi0)
             dxcar_dxqssph(4,3)=rho0*cos(theta0)*sin(phi0)
             dxcar_dxqssph(4,4)=rho0*sin(theta0)*cos(phi0)

             do a=1,4
               do b=1,4
                g0_ll_qssph(a,b)=0.0d0
                do c=1,4
                 do d=1,4
                  g0_ll_qssph(a,b)=g0_ll_qssph(a,b)
     &                        +dxcar_dxqssph(c,a)*dxcar_dxqssph(d,b)
     &                          *g0_ll(c,d)
                 end do
                end do
               end do
             end do

                call calc_g0uu(g0_ll_qssph(1,1),
     &              g0_ll_qssph(1,2),
     &              g0_ll_qssph(1,3),
     &              g0_ll_qssph(1,4),
     &              g0_ll_qssph(2,2),
     &              g0_ll_qssph(2,3),
     &              g0_ll_qssph(2,4),
     &              g0_ll_qssph(3,3),
     &              g0_ll_qssph(3,4),
     &              g0_ll_qssph(4,4),
     &              g0_uu_qssph(1,1),
     &              g0_uu_qssph(1,2),
     &              g0_uu_qssph(1,3),
     &              g0_uu_qssph(1,4),
     &              g0_uu_qssph(2,2),
     &              g0_uu_qssph(2,3),
     &              g0_uu_qssph(2,4),
     &              g0_uu_qssph(3,3),
     &              g0_uu_qssph(3,4),
     &              g0_uu_qssph(4,4),
     &              detg0_qssph0)


             do a=1,3
               do b=a+1,4
                 g0_uu_qssph(b,a)=g0_uu_qssph(a,b) 
               end do
             end do



!calculate R^s*[R^2*(df_dR)^2+|\barnabla f|^2+f^2]*(R/rho)^2*dR_drho
!To be integrated in rho^2 sin(theta)drho dtheta dphi, or dx dy dz

            if (hnorm_argtype.eq.2) then

             write(*,*) "WARNING - sqrth1spnormdensity_func:"//
     & " invalid value of hnorm_argtype: hnorm_argtype=", hnorm_argtype
             write(*,*) " hnorm_argtype can only be 0,1,3,4,5"
             write(*,*) "H_AdS^(1,sp) norm density for "//
     &                      "f=(1-rho^2)^2*f1_tt not implemented"//
     &                      " (third derivatives needed)"
             write(*,*) "Not computing the norm"
             return
            end if

            if (hnorm_argtype.eq.0) then 
            !compute norm for f=(1-rho^2)^2*f1

             f0=(1-rho0**2)**2*f10

             !compute derivatives of f
             call derf_dxsph(f1_np1,f1_n,f1_nm1,
     &         f1_t,f1_rho,f1_theta,f1_phi,
     &         x,y,z,dt,i,j,k,chr,ex,Nx,Ny,Nz,'f1')
             df1_dxqssph(1)=f1_t
             df1_dxqssph(2)=f1_rho
             df1_dxqssph(3)=f1_theta
             df1_dxqssph(4)=f1_phi

             df_dRad=
     &        -(1-rho0**2)**3*(4*rho0*f10
     &          -(1-rho0**2)*f1_rho)/(2*(1+rho0**2))
             df_dtheta=(1-rho0**2)**2*f1_theta
             df_dphi  =(1-rho0**2)**2*f1_phi

            else if (hnorm_argtype.eq.1) then
            !compute norm for f=(1-rho^2)^2*f1_t

             !compute derivatives of f
             call der2f_dxsphdxsph(f1_np1,f1_n,f1_nm1,
     &                df1_dxqssph,d2f1_dxqssphdxqssph,
     &                x,y,z,dt,i,j,k,chr,ex,Nx,Ny,Nz,'f1')

             f0=(1-rho0**2)**2*df1_dxqssph(1)
             
             df_dRad=
     &        -(1-rho0**2)**3*(4*rho0*df1_dxqssph(1)
     &          -(1-rho0**2)*d2f1_dxqssphdxqssph(1,2))/(2*(1+rho0**2))
             df_dtheta=(1-rho0**2)**2*d2f1_dxqssphdxqssph(1,3)
             df_dphi  =(1-rho0**2)**2*d2f1_dxqssphdxqssph(1,4)

            else if (hnorm_argtype.eq.3) then
            !compute norm for f=m_1[(1-rho^2)^2*f1]=(1-rho^2)^2*m1[f1], 
            !where m_1=d/dphi is the first generator of the asymptotic rotation symmetry group

             !compute derivatives of f
             call der2f_dxsphdxsph(f1_np1,f1_n,f1_nm1,
     &                df1_dxqssph,d2f1_dxqssphdxqssph,
     &                x,y,z,dt,i,j,k,chr,ex,Nx,Ny,Nz,'f1')

             f0=(1-rho0**2)**2*df1_dxqssph(4)
             
             df_dRad=
     &        -(1-rho0**2)**3*(4*rho0*df1_dxqssph(4)
     &          -(1-rho0**2)*d2f1_dxqssphdxqssph(2,4))
     &          /(2*(1+rho0**2))
             df_dtheta=(1-rho0**2)**2*d2f1_dxqssphdxqssph(3,4)
             df_dphi  =(1-rho0**2)**2*d2f1_dxqssphdxqssph(4,4)

            else if (hnorm_argtype.eq.4) then
            !compute norm for f=m_2[(1-rho^2)^2*f1]=(1-rho^2)^2*m_2[f1], 
            !where m_2=-sin(phi)d/dtheta-cot(theta)cos(phi)d/dphi is 
            !the second generator of the asymptotic rotation symmetry group

             !compute derivatives of f
             call der2f_dxsphdxsph(f1_np1,f1_n,f1_nm1,
     &                df1_dxqssph,d2f1_dxqssphdxqssph,
     &                x,y,z,dt,i,j,k,chr,ex,Nx,Ny,Nz,'f1')

             f0=(1-rho0**2)**2*
     &                  (-sin(phi0)*df1_dxqssph(3)
     &           -(cos(theta0)/sin(theta0))*cos(phi0)*df1_dxqssph(4))
             
             df_dRad=
     &        -(1-rho0**2)**3*(4*rho0*
     &          (-sin(phi0)*df1_dxqssph(3)
     &           -(cos(theta0)/sin(theta0))*cos(phi0)*df1_dxqssph(4))
     &          -(1-rho0**2)*
     &         (-sin(phi0)*d2f1_dxqssphdxqssph(2,3)
     &          -(cos(theta0)/sin(theta0))*cos(phi0)*
     &          d2f1_dxqssphdxqssph(2,4)))
     &       /(2*(1+rho0**2))

             df_dtheta=(1-rho0**2)**2*
     &          ((cos(phi0)/sin(theta0)**2)*df1_dxqssph(4)
     &         -(cos(theta0)/sin(theta0))*cos(phi0)
     &              *d2f1_dxqssphdxqssph(3,4)
     &         -sin(phi0)*d2f1_dxqssphdxqssph(3,3) )
             df_dphi=(1-rho0**2)**2*
     &          ((cos(theta0)/sin(theta0))*sin(phi0)
     &                      *df1_dxqssph(4)
     &         -(cos(theta0)/sin(theta0))*cos(phi0)
     &              *d2f1_dxqssphdxqssph(4,4)
     &         -cos(phi0)*df1_dxqssph(3)
     &         -sin(phi0)*d2f1_dxqssphdxqssph(3,4) )


            else if (hnorm_argtype.eq.5) then
            !compute norm for f=m_3[(1-rho^2)^2*f1]=(1-rho^2)^2*m_3[f1], 
            !where m_3=cos(phi)d/dtheta-cot(theta)sin(phi)d/dphi is 
            !the third generator of the asymptotic rotation symmetry group

             !compute derivatives of f
             call der2f_dxsphdxsph(f1_np1,f1_n,f1_nm1,
     &                df1_dxqssph,d2f1_dxqssphdxqssph,
     &                x,y,z,dt,i,j,k,chr,ex,Nx,Ny,Nz,'f1')

             f0=(1-rho0**2)**2*
     &             (cos(phi0)*df1_dxqssph(3)
     &         -(cos(theta0)/sin(theta0))*sin(phi0)*df1_dxqssph(4))
             
             df_dRad=
     &        -(1-rho0**2)**3*(4*rho0*(
     &             cos(phi0)*df1_dxqssph(3)
     &           -(cos(theta0)/sin(theta0))*sin(phi0)*df1_dxqssph(4))
     &          -(1-rho0**2)*(
     &             cos(phi0)*d2f1_dxqssphdxqssph(2,3)
     &           -(cos(theta0)/sin(theta0))*sin(phi0)*
     &             d2f1_dxqssphdxqssph(2,4)))
     &          /(2*(1+rho0**2))
             df_dtheta=(1-rho0**2)**2*
     &          ((sin(phi0)/sin(theta0)**2)*df1_dxqssph(4)
     &         -(cos(theta0)/sin(theta0))*sin(phi0)
     &              *d2f1_dxqssphdxqssph(3,4)
     &         +cos(phi0)*d2f1_dxqssphdxqssph(3,3) )
             df_dphi  =(1-rho0**2)**2*
     &           (-(cos(theta0)/sin(theta0))*cos(phi0)
     &                      *df1_dxqssph(4)
     &         -(cos(theta0)/sin(theta0))*sin(phi0)
     &              *d2f1_dxqssphdxqssph(4,4)
     &         -sin(phi0)*df1_dxqssph(3)
     &         +cos(phi0)*d2f1_dxqssphdxqssph(3,4) )


            else
             write(*,*) "WARNING - sqrth1spnormdensity_func:"//
     &  " invalid value of hnorm_argtype: hnorm_argtype=", hnorm_argtype
             write(*,*) " hnorm_argtype can only be 0,1,3,4,5"
             return
            end if

              h1spnormdensity_f0=
     &             Rad0**2*df_dRad**2
     &            +(g0_uu_qssph(3,3)*
     &               df_dtheta**2+ 
     &             2*g0_uu_qssph(3,4)*
     &               df_dtheta*df_dphi+ 
     &             g0_uu_qssph(4,4)*
     &               df_dphi**2)
     &            +f0**2

            h1spnormdensity_f0=(Rad0**sp)*
     &            h1spnormdensity_f0*
     &            Rad0**2/rho0**2*
     &            dRad_drho

!TEST logarithm dependence
!        h1spnormdensity_f0=1/abs(log(ct+dt))


!Take square root, so we can use DV calculation of L2-norm, and then square it, to compute the H_AdS(1,s) norm.
                sqrth1spnormdensity_f(i,j,k)=
     &             sqrt(h1spnormdensity_f0)

                 !the phi coordinate is not defined at y=z=0, i.e., theta=0,PI 
                !(not surprising since spherical coordinates are not valid everywhere on the sphere), 
                !hence we set sqrth1spnormdensity_f(i,j,k) to 0 at these points

                if ((abs(y0).lt.10.0d0**(-10)).and.
     &          (abs(z0).lt.10.0d0**(-10))) then
                    sqrth1spnormdensity_f(i,j,k)=0.0d0
                end if

                if (ltrace) then
                 if ((abs(x0).gt.10.0d0**(-1.0d0)).and.
     &               (abs(y0).gt.10.0d0**(-1.0d0)).and.
     &               (abs(z0).gt.10.0d0**(-1.0d0)).and.
     &                (rho0.lt.0.9)  ) then
                   call df1_int(f1_np1,f1_n,f1_nm1,
     &                     f1_t,f1_x,f1_y,f1_z,
     &                     x,y,z,dt,i,j,k,chr,ex,Nx,Ny,Nz,'f1')
                   write(*,*) 'sqrth1spnormdensity_func '
                   write(*,*) 'sp,hnorm_argtype= ',sp,hnorm_argtype
                   write(*,*) 'at i,j,k= ', i,j,k
                   write(*,*) 'i.e., x,y,z=', x(i),y(j),z(k)
                   write(*,*) 'f10=', f10
                   write(*,*) 'f1_t=', f1_t
                   write(*,*) 'f1_x=', f1_x
                   write(*,*) 'f1_y=', f1_y
                   write(*,*) 'f1_z=', f1_z
                   write(*,*) 'f1_t=', f1_t
                   write(*,*) 'f1_rho=', df1_dxqssph(2)
                   write(*,*) 'f1_theta=', df1_dxqssph(3)
                   write(*,*) 'f1_phi=', df1_dxqssph(4)
            write(*,*) 'f1_rhorho=', d2f1_dxqssphdxqssph(2,2)
            write(*,*) 'f1_rhotheta=', d2f1_dxqssphdxqssph(2,3)
            write(*,*) 'f1_rhophi=', d2f1_dxqssphdxqssph(2,4)
            write(*,*) 'f1_thetatheta=', d2f1_dxqssphdxqssph(3,3)
            write(*,*) 'f1_thetaphi=', d2f1_dxqssphdxqssph(3,4)
            write(*,*) 'f1_phiphi=', d2f1_dxqssphdxqssph(4,4)
                   write(*,*) 'f0=', f0
                   write(*,*) 'df_dRad=',df_dRad
                   write(*,*) 'df_dtheta=',df_dtheta
                   write(*,*) 'df_dphi=',df_dphi
                   write(*,*) 'nablafq=',
     &              (g0_uu_qssph(3,3)*
     &               df_dtheta**2+ 
     &             2*g0_uu_qssph(3,4)*
     &               df_dtheta*df_dphi+ 
     &             g0_uu_qssph(4,4)*
     &               df_dphi**2)
                   write(*,*) ' sqrth1spnormdensity_f = ',
     &                  sqrth1spnormdensity_f(i,j,k)
                   return
                 end if
                end if


           else !i.e., the point is either excised or inside the apparent horizon
               sqrth1spnormdensity_f(i,j,k)=0.0d0
           end if


          end do
         end do
        end do


        return
        end



c-----------------------------------------------------------------------
c Compute H_AdS^(2,sp) of scalar field f=(1-rho^2)^2*f1
c-----------------------------------------------------------------------
        subroutine sqrth2spnormdensity_func(
     &                  sqrth2spnormdensity_f,
     &                  sp,hnorm_argtype,
     &                  f1_np1,f1_n,f1_nm1,
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
     &                  AH_R,AH_xc,
     &                  AH_Ntheta,AH_Nphi,
     &                  x,y,z,dt,ct,chr,L,ex,Nx,Ny,Nz,
     &                  phys_bdy,ghost_width,
     &                  ief_bh_r0,a_rot,kerrads_background)

        implicit none

        logical is_nan
        real*8 ief_bh_r0,a_rot,M0,M0_min
        real*8 kerrads_background

        logical calc_der,calc_adv_quant
        data calc_der/.true./
        data calc_adv_quant/.false./

        integer Nx,Ny,Nz
        integer i,j,k
        integer phys_bdy(6),ghost_width(6)
        real*8 chr(Nx,Ny,Nz),ex
        real*8 x(Nx),y(Ny),z(Nz),dt,ct,L
        real*8 lambda4
        integer AH_Ntheta,AH_Nphi
        real*8 AH_R(AH_Ntheta,AH_Nphi)
        real*8 AH_xc(3)
        real*8 xp,yp,zp,rhop,thetap,phip
        integer i0,j0
        real*8 dahtheta,dahphi
        real*8 ft,fp
        real*8 AH_Rp,rhostop
        real*8 f1_np1(Nx,Ny,Nz),f1_n(Nx,Ny,Nz),f1_nm1(Nx,Ny,Nz)
        real*8 sqrth2spnormdensity_f(Nx,Ny,Nz)
        real*8 h2spnormdensity_f0
        integer sp,hnorm_argtype
        real*8 f10,f0
        real*8 df1_dxqssph(4),d2f1_dxqssphdxqssph(4,4)
        real*8 f1_t,f1_rho,f1_theta,f1_phi
        real*8 f1_tt,f1_trho,f1_ttheta,f1_tphi
        real*8 f1_rhorho,f1_rhotheta,f1_rhophi
        real*8 f1_thetatheta,f1_thetaphi
        real*8 f1_phiphi
        real*8 df_dRad,df_dtheta,df_dphi
        real*8 d2f_dRaddRad,d2f_dRaddtheta,d2f_dRaddphi
        real*8 gamma2dim_ull_qssph(2,2,2)
        real*8 nablanablaf_qssph_trhocons_xx(2,2)

        logical ltrace
        parameter (ltrace=.false.)

        integer is,ie,js,je,ks,ke

        integer a,b,ap,bp,c,d,e,f
        integer as,bs,aps,bps,cs,ds
        integer ic,jc,kc
        real*8 efe_ires(4,4)

        real*8 PI
        parameter (PI=3.141592653589793d0)

        real*8 dx,dy,dz
        real*8 x0,y0,z0
        real*8 rho0,theta0,phi0   
        real*8 Rad0
        real*8 drho_dRad,dRad_drho

        real*8 dxcar_dxqssph(4,4)
        real*8 d2xcar_dxqssphdxqssph(4,4,4)

        real*8 gb_tt_np1(Nx,Ny,Nz),gb_tt_n(Nx,Ny,Nz),gb_tt_nm1(Nx,Ny,Nz)
        real*8 gb_tx_np1(Nx,Ny,Nz),gb_tx_n(Nx,Ny,Nz),gb_tx_nm1(Nx,Ny,Nz)
        real*8 gb_ty_np1(Nx,Ny,Nz),gb_ty_n(Nx,Ny,Nz),gb_ty_nm1(Nx,Ny,Nz)
        real*8 gb_tz_np1(Nx,Ny,Nz),gb_tz_n(Nx,Ny,Nz),gb_tz_nm1(Nx,Ny,Nz)
        real*8 gb_xx_np1(Nx,Ny,Nz),gb_xx_n(Nx,Ny,Nz),gb_xx_nm1(Nx,Ny,Nz)
        real*8 gb_xy_np1(Nx,Ny,Nz),gb_xy_n(Nx,Ny,Nz),gb_xy_nm1(Nx,Ny,Nz)
        real*8 gb_xz_np1(Nx,Ny,Nz),gb_xz_n(Nx,Ny,Nz),gb_xz_nm1(Nx,Ny,Nz)
        real*8 gb_yy_np1(Nx,Ny,Nz),gb_yy_n(Nx,Ny,Nz),gb_yy_nm1(Nx,Ny,Nz)
        real*8 gb_yz_np1(Nx,Ny,Nz),gb_yz_n(Nx,Ny,Nz),gb_yz_nm1(Nx,Ny,Nz)
        real*8 gb_zz_np1(Nx,Ny,Nz),gb_zz_n(Nx,Ny,Nz),gb_zz_nm1(Nx,Ny,Nz)
        real*8 Hb_t_np1(Nx,Ny,Nz),Hb_t_n(Nx,Ny,Nz),Hb_t_nm1(Nx,Ny,Nz)
        real*8 Hb_x_np1(Nx,Ny,Nz),Hb_x_n(Nx,Ny,Nz),Hb_x_nm1(Nx,Ny,Nz)
        real*8 Hb_y_np1(Nx,Ny,Nz),Hb_y_n(Nx,Ny,Nz),Hb_y_nm1(Nx,Ny,Nz)
        real*8 Hb_z_np1(Nx,Ny,Nz),Hb_z_n(Nx,Ny,Nz),Hb_z_nm1(Nx,Ny,Nz)
        real*8 phi1_np1(Nx,Ny,Nz),phi1_n(Nx,Ny,Nz),phi1_nm1(Nx,Ny,Nz)

        !--------------------------------------------------------------
        ! variables for tensor manipulations 
        !(indices are t,x,y,theta,phi)
        !--------------------------------------------------------------
        real*8 gads_ll(4,4),gads_uu(4,4)
        real*8 gads_ll_x(4,4,4),gads_uu_x(4,4,4)
        real*8 gads_ll_xx(4,4,4,4)
        real*8 Hads_l(4)
        real*8 phi1ads, phi1ads_x(4)

        real*8 g0_ll(4,4),g0_uu(4,4),detg0
        real*8 g0_ll_x(4,4,4),g0_uu_x(4,4,4),g0_ll_xx(4,4,4,4)
        real*8 h0_ll(4,4),h0_uu(4,4)
        real*8 h0_ll_x(4,4,4),h0_uu_x(4,4,4),h0_ll_xx(4,4,4,4)
        real*8 gamma_ull(4,4,4),gamma_ull_x(4,4,4,4)
        real*8 riemann_ulll(4,4,4,4)
        real*8 ricci_ll(4,4),ricci_lu(4,4),ricci
        real*8 einstein_ll(4,4),set_ll(4,4)
        real*8 A_l(4),A_l_x(4,4)
        real*8 phi10_x(4),phi10_xx(4,4)

        real*8 detg0_qssph0
        real*8 g0_ll_qssph(4,4),g0_uu_qssph(4,4)
        real*8 g0_ll_qssph_x(4,4,4)




        ! initialize fixed-size variables
        data i,j,k,is,ie,js,je,ks,ke/0,0,0,0,0,0,0,0,0/
        data ic,jc,kc/0,0,0/
        data a,b,ap,bp,c,d,e,f/0,0,0,0,0,0,0,0/
        data as,bs,aps,bps,cs,ds/0,0,0,0,0,0/

        data dx,dy,dz/0.0,0.0,0.0/
        data x0,y0,z0,rho0/0.0,0.0,0.0,0.0/    


        data gads_ll,gads_uu/16*0.0,16*0.0/
        data gads_ll_x,gads_uu_x/64*0.0,64*0.0/
        data gads_ll_xx/256*0.0/
        data Hads_l/4*0.0/
        data phi1ads_x/4*0.0/

        data g0_ll,g0_uu/16*0.0,16*0.0/
        data g0_ll_x,g0_uu_x/64*0.0,64*0.0/
        data g0_ll_xx/256*0.0/
        data h0_ll,h0_uu/16*0.0,16*0.0/
        data h0_ll_x,h0_uu_x/64*0.0,64*0.0/
        data h0_ll_xx/256*0.0/
        data gamma_ull/64*0.0/
        data gamma_ull_x/256*0.0/
        data riemann_ulll/256*0.0/
        data ricci_ll,ricci_lu/16*0.0,16*0.0/
        data einstein_ll,set_ll/16*0.0,16*0.0/
        data A_l,A_l_x/4*0.0,16*0.0/
        data phi10_x,phi10_xx/4*0.0,16*0.0/
        data dxcar_dxqssph/16*0.0/

        data g0_ll_qssph/16*0.0/
        data g0_uu_qssph/16*0.0/
        data g0_ll_qssph_x/64*0.0/


        data df1_dxqssph,d2f1_dxqssphdxqssph/4*0.0,16*0.0/
        data dxcar_dxqssph,d2xcar_dxqssphdxqssph/16*0.0,64*0.0/
        data nablanablaf_qssph_trhocons_xx/4*0.0/
        data gamma2dim_ull_qssph/8*0.0/



!----------------------------------------------------------------------


      ! Black hole mass
        M0=ief_bh_r0/2
      ! Minimum black hole mass. For M0 below this value, there is a naked singularity
        M0_min=((2*(1 + a_rot**2/L**2) + Sqrt((1 + a_rot**2/L**2)**2 
     &   + (12*a_rot**2)/L**2))*Sqrt(-1 + Sqrt((1 + a_rot**2/L**2)**2 
     &   + (12*a_rot**2)/L**2) - a_rot**2/L**2))/(3.*Sqrt(6.))

        if (a_rot.ge.L) then
         write (*,*) "sqrth1spnormdensity_phi: ERROR in choice of Kerr-AdS initial parameters: 
     &   the rotation parameter a must be smaller than the AdS radius L"
          write (*,*) "a_rot,L=",a_rot,L
          stop
        end if

        if ((abs(M0).gt.10.0d0**(-10))
     &     .and.(M0.le.M0_min)) then
          write (*,*) "sqrth1spnormdensity_phi: ERROR in choice of Kerr-AdS initial parameters: 
     &      the black hole mass M0=2*r0 must be larger
     &      than the M0_min value"
          write (*,*) "M0,M0_min=",M0,M0_min
          stop
        end if


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
              theta0=acos(x0/rho0)
              if (z0.lt.0) then
               phi0=atan2(z0,y0)+2*PI
              else
               phi0=atan2(z0,y0)
              end if
              !uncompactified quasi-sherical Kerr-Schild radial coordinate
              Rad0=2*rho0/(1-rho0**2)
              drho_dRad=((1-rho0**2)**2)/(2*(1+rho0**2))
              dRad_drho=(-1 + rho0)**(-2) 
     &           + (1 + rho0)**(-2)

              !determine radius at which we stop computing the norms. This is the apparent horizon (AH) radius.

              !Cartesian coordinates of grid point p=(i,j,k) w.r.t. the centre of the AH
              xp=x0-AH_xc(1)
              yp=y0-AH_xc(2)
              zp=z0-AH_xc(3)

              !The center of the AH xp=yp=zp=0 will be excluded or included depending on the excision mask chr(i,j,k), 
              !so no need to force inclusion or exclusion via the value of rhostop. 
              !In fact, xp=yp=zp=0 is also the centre of the grid if there is no AH (also if the AH is centred at the origin), 
              !so we should include it if we want to compute these norms in the case of no horizon.
              !For this reason, it's best to leave the inclusion or exclusion of that point to the excision mask chr(i,j,k).
              if ((abs(xp).lt.10.0d0**(-12)).and.
     &            (abs(yp).lt.10.0d0**(-12)).and.  
     &            (abs(zp).lt.10.0d0**(-12)) ) then
                rhostop=-1
              else
                !corresponding spherical coordinates of point p
                !thetap,phip is in general NOT on the AH_Ntheta,AH_Nphi grid.
                rhop=sqrt(xp**2+yp**2+zp**2)
                thetap=acos(xp/rhop)
                if (zp.lt.0) then
                    phip=atan2(zp,yp)+2*PI
                else
                    phip=atan2(zp,yp)
                end if

                ! i0 is the index of the point with chi value closest to chip on the AH_Nchi grid (indices going from 0 to AH_Ntheta-1),
                ! and larger than chip (except if thetap=Pi, then it is smaller)
                ! similarly for j0
                dahtheta=PI/(AH_Ntheta-1)
                dahphi=2*PI/(AH_Nphi-1)
                i0=floor(thetap/dahtheta+1)
                if (AH_Ntheta-1<i0) i0=AH_Ntheta-1
                j0=floor(phip/dahphi+1)
                if (AH_Nphi-1<j0) j0=AH_Nphi-1
                ft=((thetap-(i0-1)*dahtheta))/dahtheta
                fp=((phip-(j0-1)*dahphi))/dahphi

                !bi-linear interpolated value of AH radius at (thetap,phip)
                AH_Rp=  AH_R(i0,j0)*(1-ft)*(1-fp)+
     &                  AH_R(i0+1,j0)*(ft)*(1-fp)+
     &                  AH_R(i0,j0+1)*(1-ft)*(fp)+
     &                  AH_R(i0+1,j0+1)*(ft)*(fp)

                if (AH_Rp.lt.1) then 
                    rhostop=AH_Rp
                else
                    rhostop=-1
                end if
              end if



            if ((rho0.gt.rhostop).and.
     &         (chr(i,j,k).ne.ex)) then

                ! set f1 value
                f10=f1_n(i,j,k)

    !compute metric and its derivatives 
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
     &                gads_ll,gads_uu,
     &                gads_ll_x,gads_uu_x,gads_ll_xx,
     &                h0_ll,h0_uu,h0_ll_x,h0_uu_x,h0_ll_xx,
     &                A_l,A_l_x,Hads_l,
     &                gamma_ull,gamma_ull_x,
     &                riemann_ulll,ricci_ll,ricci_lu,ricci,
     &                einstein_ll,set_ll,
     &                phi10_x,phi10_xx,
     &                x,y,z,dt,chr,L,ex,Nx,Ny,Nz,i,j,k,
     &                ief_bh_r0,a_rot,kerrads_background,
     &                calc_der,calc_adv_quant)

!define transformation matrix between Cartesian coordinates to compactified (quasi-)spherical coordinates, 
!        !e.g. dxcar_dxqssph(3,2)=dy/drho=sin(theta0)*cos(phi0)

             dxcar_dxqssph(1,1)=1
             dxcar_dxqssph(1,2)=0
             dxcar_dxqssph(1,3)=0
             dxcar_dxqssph(1,4)=0
             dxcar_dxqssph(2,1)=0
             dxcar_dxqssph(2,2)=cos(theta0)
             dxcar_dxqssph(2,3)=-rho0*sin(theta0)
             dxcar_dxqssph(2,4)=0
             dxcar_dxqssph(3,1)=0
             dxcar_dxqssph(3,2)=sin(theta0)*cos(phi0)
             dxcar_dxqssph(3,3)=rho0*cos(theta0)*cos(phi0)
             dxcar_dxqssph(3,4)=-rho0*sin(theta0)*sin(phi0)
             dxcar_dxqssph(4,1)=0
             dxcar_dxqssph(4,2)=sin(theta0)*sin(phi0)
             dxcar_dxqssph(4,3)=rho0*cos(theta0)*sin(phi0)
             dxcar_dxqssph(4,4)=rho0*sin(theta0)*cos(phi0)


             d2xcar_dxqssphdxqssph(1,1,1)=0
             d2xcar_dxqssphdxqssph(1,1,2)=0
             d2xcar_dxqssphdxqssph(1,1,3)=0
             d2xcar_dxqssphdxqssph(1,1,4)=0
             d2xcar_dxqssphdxqssph(1,2,1)=0
             d2xcar_dxqssphdxqssph(1,2,2)=0
             d2xcar_dxqssphdxqssph(1,2,3)=0
             d2xcar_dxqssphdxqssph(1,2,4)=0
             d2xcar_dxqssphdxqssph(1,3,1)=0
             d2xcar_dxqssphdxqssph(1,3,2)=0
             d2xcar_dxqssphdxqssph(1,3,3)=0
             d2xcar_dxqssphdxqssph(1,3,4)=0
             d2xcar_dxqssphdxqssph(1,4,1)=0
             d2xcar_dxqssphdxqssph(1,4,2)=0
             d2xcar_dxqssphdxqssph(1,4,3)=0
             d2xcar_dxqssphdxqssph(1,4,4)=0

             d2xcar_dxqssphdxqssph(2,1,1)=0
             d2xcar_dxqssphdxqssph(2,1,2)=0
             d2xcar_dxqssphdxqssph(2,1,3)=0
             d2xcar_dxqssphdxqssph(2,1,4)=0
             d2xcar_dxqssphdxqssph(2,2,1)=0
             d2xcar_dxqssphdxqssph(2,2,2)=0
             d2xcar_dxqssphdxqssph(2,2,3)=-sin(theta0)
             d2xcar_dxqssphdxqssph(2,2,4)=0
             d2xcar_dxqssphdxqssph(2,3,1)=0
             d2xcar_dxqssphdxqssph(2,3,2)=-sin(theta0)
             d2xcar_dxqssphdxqssph(2,3,3)=-rho0*cos(theta0)
             d2xcar_dxqssphdxqssph(2,3,4)=0
             d2xcar_dxqssphdxqssph(2,4,1)=0
             d2xcar_dxqssphdxqssph(2,4,2)=0
             d2xcar_dxqssphdxqssph(2,4,3)=0
             d2xcar_dxqssphdxqssph(2,4,4)=0

             d2xcar_dxqssphdxqssph(3,1,1)=0
             d2xcar_dxqssphdxqssph(3,1,2)=0
             d2xcar_dxqssphdxqssph(3,1,3)=0
             d2xcar_dxqssphdxqssph(3,1,4)=0
             d2xcar_dxqssphdxqssph(3,2,1)=0
             d2xcar_dxqssphdxqssph(3,2,2)=0
             d2xcar_dxqssphdxqssph(3,2,3)=cos(theta0)*cos(phi0)
             d2xcar_dxqssphdxqssph(3,2,4)=-sin(theta0)*sin(phi0)
             d2xcar_dxqssphdxqssph(3,3,1)=0
             d2xcar_dxqssphdxqssph(3,3,2)=cos(theta0)*cos(phi0)
             d2xcar_dxqssphdxqssph(3,3,3)=-rho0*sin(theta0)*cos(phi0)
             d2xcar_dxqssphdxqssph(3,3,4)=-rho0*cos(theta0)*sin(phi0)
             d2xcar_dxqssphdxqssph(3,4,1)=0
             d2xcar_dxqssphdxqssph(3,4,2)=-sin(theta0)*sin(phi0)
             d2xcar_dxqssphdxqssph(3,4,3)=-rho0*cos(theta0)*sin(phi0)
             d2xcar_dxqssphdxqssph(3,4,4)=-rho0*sin(theta0)*cos(phi0)

             d2xcar_dxqssphdxqssph(4,1,1)=0
             d2xcar_dxqssphdxqssph(4,1,2)=0
             d2xcar_dxqssphdxqssph(4,1,3)=0
             d2xcar_dxqssphdxqssph(4,1,4)=0
             d2xcar_dxqssphdxqssph(4,2,1)=0
             d2xcar_dxqssphdxqssph(4,2,2)=0
             d2xcar_dxqssphdxqssph(4,2,3)=cos(theta0)*sin(phi0)
             d2xcar_dxqssphdxqssph(4,2,4)=sin(theta0)*cos(phi0)
             d2xcar_dxqssphdxqssph(4,3,1)=0
             d2xcar_dxqssphdxqssph(4,3,2)=cos(theta0)*sin(phi0)
             d2xcar_dxqssphdxqssph(4,3,3)=-rho0*sin(theta0)*sin(phi0)
             d2xcar_dxqssphdxqssph(4,3,4)=rho0*cos(theta0)*cos(phi0)
             d2xcar_dxqssphdxqssph(4,4,1)=0
             d2xcar_dxqssphdxqssph(4,4,2)=sin(theta0)*cos(phi0)
             d2xcar_dxqssphdxqssph(4,4,3)=rho0*cos(theta0)*cos(phi0)
             d2xcar_dxqssphdxqssph(4,4,4)=-rho0*sin(theta0)*sin(phi0)

             do a=1,4
               do b=1,4
                g0_ll_qssph(a,b)=0.0d0
                do c=1,4
                 g0_ll_qssph_x(a,b,c)=0.0d0
                 do d=1,4
                  g0_ll_qssph(a,b)=g0_ll_qssph(a,b)
     &                        +dxcar_dxqssph(c,a)*dxcar_dxqssph(d,b)
     &                          *g0_ll(c,d)
                  do e=1,4
                   g0_ll_qssph_x(a,b,c)=g0_ll_qssph_x(a,b,c)
     &              +d2xcar_dxqssphdxqssph(d,c,a)*dxcar_dxqssph(e,b)
     &              *g0_ll(d,e)
     &              +dxcar_dxqssph(d,a)*d2xcar_dxqssphdxqssph(e,c,b)
     &              *g0_ll(d,e)
                   do f=1,4
                    g0_ll_qssph_x(a,b,c)=
     &               g0_ll_qssph_x(a,b,c)
     &               +dxcar_dxqssph(d,a)*
     &                dxcar_dxqssph(e,b)*
     &                dxcar_dxqssph(f,c)*
     &                g0_ll_x(d,e,f)
                   end do
                  end do
                 end do
                end do
               end do
             end do

                call calc_g0uu(g0_ll_qssph(1,1),
     &              g0_ll_qssph(1,2),
     &              g0_ll_qssph(1,3),
     &              g0_ll_qssph(1,4),
     &              g0_ll_qssph(2,2),
     &              g0_ll_qssph(2,3),
     &              g0_ll_qssph(2,4),
     &              g0_ll_qssph(3,3),
     &              g0_ll_qssph(3,4),
     &              g0_ll_qssph(4,4),
     &              g0_uu_qssph(1,1),
     &              g0_uu_qssph(1,2),
     &              g0_uu_qssph(1,3),
     &              g0_uu_qssph(1,4),
     &              g0_uu_qssph(2,2),
     &              g0_uu_qssph(2,3),
     &              g0_uu_qssph(2,4),
     &              g0_uu_qssph(3,3),
     &              g0_uu_qssph(3,4),
     &              g0_uu_qssph(4,4),
     &              detg0_qssph0)


        do a=1,3
          do b=a+1,4
            g0_uu_qssph(b,a)=g0_uu_qssph(a,b) 
          end do
        end do

        if (hnorm_argtype.ne.0) then
         write(*,*) "WARNING - sqrth2spnormdensity_func:"//
     & " invalid value of hnorm_argtype: hnorm_argtype=", hnorm_argtype
         write(*,*) " hnorm_argtype can only be 0"
         write(*,*) "H_AdS^(2,sp) norm density for "//
     &     "derivatives of f=(1-rho^2)^2*f1 not implemented "//
     &     "(third derivatives needed)"
         write(*,*) "Not computing the norm"
         return
        end if


!compute norm for f=(1-rho^2)^2*f1

        f0=(1-rho0**2)**2*f10

        call der2f_dxsphdxsph(f1_np1,f1_n,f1_nm1,
     &                df1_dxqssph,d2f1_dxqssphdxqssph,
     &                x,y,z,dt,i,j,k,chr,ex,Nx,Ny,Nz,'f1')

            f1_rho         = df1_dxqssph(2)
            f1_theta       = df1_dxqssph(3)
            f1_phi         = df1_dxqssph(4)

            f1_rhorho      = d2f1_dxqssphdxqssph(2,2)
            f1_rhotheta    = d2f1_dxqssphdxqssph(2,3)
            f1_rhophi      = d2f1_dxqssphdxqssph(2,4)
            f1_thetatheta  = d2f1_dxqssphdxqssph(3,3)
            f1_thetaphi    = d2f1_dxqssphdxqssph(3,4)
            f1_phiphi      = d2f1_dxqssphdxqssph(4,4)



            !compute derivatives of f recalling f=(1-rho^2)^2*f1
            df_dRad=
     &       -(1-rho0**2)**3*(4*rho0*f10
     &        -(1-rho0**2)*f1_rho)/(2*(1+rho0**2))
            df_dtheta=(1-rho0**2)**2*f1_theta
             df_dphi  =(1-rho0**2)**2*f1_phi
            d2f_dRaddRad=   
     &       (((1-rho0**2)**4)/(4*(1+rho0**2)**3))
     &       *(4*(-1+8*rho0**2+5*rho0**4)*f10
     &         -(1-rho0**2)*
     &          (2*rho0*(7+5*rho0**2)*f1_rho
     &           -(1-rho0**4)*f1_rhorho)
     &        )
            d2f_dRaddtheta=
     &       -(1-rho0**2)**3*(4*rho0*f1_theta
     &        -(1-rho0**2)*f1_rhotheta)/(2*(1+rho0**2))
            d2f_dRaddphi=
     &       -(1-rho0**2)**3*(4*rho0*f1_phi
     &        -(1-rho0**2)*f1_rhophi)/(2*(1+rho0**2))

!compute second covariant derivative of f for surfaces of constant t,rho

            do a=3,4
             do b=3,4
              do c=3,4
               as=a-2
               bs=b-2
               cs=c-2
               gamma2dim_ull_qssph(as,bs,cs)=0
               do d=3,4
                gamma2dim_ull_qssph(as,bs,cs)=
     &            gamma2dim_ull_qssph(as,bs,cs)+
     &              0.5d0*g0_uu_qssph(a,d)
     &                   *(g0_ll_qssph_x(b,d,c)
     &                    -g0_ll_qssph_x(b,c,d)
     &                    +g0_ll_qssph_x(d,c,b))
               end do
              end do
             end do
            end do

            do a=3,4
             do b=3,4
              as=a-2
              bs=b-2
              nablanablaf_qssph_trhocons_xx(as,bs)=
     &          (1-rho0**2)**2*d2f1_dxqssphdxqssph(a,b)
              do c=3,4
                cs=c-2
                nablanablaf_qssph_trhocons_xx(as,bs)=
     &           nablanablaf_qssph_trhocons_xx(as,bs)
     &              -gamma2dim_ull_qssph(cs,as,bs)
     &              *(1-rho0**2)**2*df1_dxqssph(c)
              end do
             end do
            end do

!calculate R^s*[R^2*(df_dR)^2+|\barnabla f|^2+f^2
!               +R^4*(d2f_dRdR)^2+R^2*|\barnabla df_dR|^2+|\barnabla\barnabla f|^2)]*(R/rho)^2*dR_drho
!To be integrated in rho^2 sin(theta)drho dtheta dphi, or dx dy dz

              h2spnormdensity_f0=
     &             Rad0**2*df_dRad**2
     &            +(g0_uu_qssph(3,3)*
     &               df_dtheta**2+ 
     &             2*g0_uu_qssph(3,4)*
     &               df_dtheta*df_dphi+ 
     &             g0_uu_qssph(4,4)*
     &               df_dphi**2)
     &            +f0**2

            h2spnormdensity_f0=
     &       h2spnormdensity_f0+
     &          Rad0**4*d2f_dRaddRad**2
     &          +Rad0**2*(
     &           g0_uu_qssph(3,3)*d2f_dRaddtheta**2+
     &         2*g0_uu_qssph(3,4)*d2f_dRaddtheta*d2f_dRaddphi+
     &           g0_uu_qssph(4,4)*d2f_dRaddphi**2)
         do a=3,4
          do b=3,4
           do ap=3,4
            do bp=3,4
                as=a-2
                bs=b-2
                aps=ap-2
                bps=bp-2
                 h2spnormdensity_f0=
     &            h2spnormdensity_f0
     &            +g0_uu_qssph(a,ap)*
     &             g0_uu_qssph(b,bp)*
     &            (nablanablaf_qssph_trhocons_xx(as,bs)*
     &             nablanablaf_qssph_trhocons_xx(aps,bps) )
            end do
           end do
          end do
         end do

            h2spnormdensity_f0=(Rad0**sp)*
     &          h2spnormdensity_f0*
     &          Rad0**2/rho0**2*
     &          dRad_drho


!TEST logarithm dependence
!        h2spnormdensity_f0=1/abs(log(ct+dt))


!Take square root, so we can use DV calculation of L2-norm, and then square it, to compute the H_AdS(2,s) norm.
                sqrth2spnormdensity_f(i,j,k)=
     &             sqrt(h2spnormdensity_f0)

                 !the phi coordinate is not defined at y=z=0, i.e., theta=0,PI 
                !(not surprising since spherical coordinates are not valid everywhere on the sphere), 
                !hence we set sqrth2spnormdensity_f(i,j,k) to 0 at these points

                if ((abs(y0).lt.10.0d0**(-10)).and.
     &          (abs(z0).lt.10.0d0**(-10))) then
                    sqrth2spnormdensity_f(i,j,k)=0.0d0
                end if

                if (ltrace) then
                 if ((abs(x0).gt.10.0d0**(-1.0d0)).and.
     &               (abs(y0).gt.10.0d0**(-1.0d0)).and.
     &               (abs(z0).gt.10.0d0**(-1.0d0)).and.
     &                (rho0.lt.0.9)  ) then
                   write(*,*) 'sqrth2spnormdensity_func '
                   write(*,*) 'sp,hnorm_argtype= ',sp,hnorm_argtype
                   write(*,*) 'at i,j,k= ', i,j,k
                   write(*,*) 'i.e., x,y,z=', x(i),y(j),z(k)
                   write(*,*) 'f10=', f10
                   write(*,*) 'f1_t=', df1_dxqssph(1)
                   write(*,*) 'f1_rho=', df1_dxqssph(2)
                   write(*,*) 'f1_theta=', df1_dxqssph(3)
                   write(*,*) 'f1_phi=', df1_dxqssph(4)
                   write(*,*) 'f0=', f0
                   write(*,*) 'df_dRad=',df_dRad
                   write(*,*) 'd2f_dRaddRad=',d2f_dRaddRad
                   write(*,*) 'd2f_dRaddtheta=',d2f_dRaddtheta
                   write(*,*) 'd2f_dRaddphi=',d2f_dRaddphi
                   write(*,*) 'g0_ll_qssph_x(3,3,3)=',
     &                          g0_ll_qssph_x(3,3,3)
                   write(*,*) 'g0_ll_qssph_x(3,3,4)=',
     &                          g0_ll_qssph_x(3,3,4)
                   write(*,*) 'g0_ll_qssph_x(3,4,3)=',
     &                          g0_ll_qssph_x(3,4,3)
                   write(*,*) 'g0_ll_qssph_x(3,4,4)=',
     &                          g0_ll_qssph_x(3,4,4)
                   write(*,*) 'g0_ll_qssph_x(4,3,3)=',
     &                          g0_ll_qssph_x(4,3,3)
                   write(*,*) 'g0_ll_qssph_x(4,3,4)=',
     &                          g0_ll_qssph_x(4,3,4)
                   write(*,*) 'g0_ll_qssph_x(4,4,3)=',
     &                          g0_ll_qssph_x(4,4,3)
                   write(*,*) 'g0_ll_qssph_x(4,4,4)=',
     &                          g0_ll_qssph_x(4,4,4)
                   write(*,*) 'gamma2dim_ull_qssph(1,1,1)=',
     &                          gamma2dim_ull_qssph(1,1,1)
                   write(*,*) 'gamma2dim_ull_qssph(1,1,2)=',
     &                          gamma2dim_ull_qssph(1,1,2)
                   write(*,*) 'gamma2dim_ull_qssph(1,2,1)=',
     &                          gamma2dim_ull_qssph(1,2,1)
                   write(*,*) 'gamma2dim_ull_qssph(1,2,2)=',
     &                          gamma2dim_ull_qssph(1,2,2)
                   write(*,*) 'gamma2dim_ull_qssph(2,1,1)=',
     &                          gamma2dim_ull_qssph(2,1,1)
                   write(*,*) 'gamma2dim_ull_qssph(2,1,2)=',
     &                          gamma2dim_ull_qssph(2,1,2)
                   write(*,*) 'gamma2dim_ull_qssph(2,2,1)=',
     &                          gamma2dim_ull_qssph(2,2,1)
                   write(*,*) 'gamma2dim_ull_qssph(2,2,2)=',
     &                          gamma2dim_ull_qssph(2,2,2)
                   write(*,*) 'nablanablaf_qssph_trhocons_xx(1,1)=',
     &                         nablanablaf_qssph_trhocons_xx(1,1)
                   write(*,*) 'nablanablaf_qssph_trhocons_xx(1,2)=',
     &                         nablanablaf_qssph_trhocons_xx(1,2)
                   write(*,*) 'nablanablaf_qssph_trhocons_xx(2,1)=',
     &                         nablanablaf_qssph_trhocons_xx(2,1)
                   write(*,*) 'nablanablaf_qssph_trhocons_xx(2,2)=',
     &                         nablanablaf_qssph_trhocons_xx(2,2)
                   write(*,*) ' sqrth2spnormdensity_f = ',
     &                  sqrth2spnormdensity_f(i,j,k)
                   return
                 end if
                end if


            else !i.e., the point is either excised or inside the apparent horizon
               sqrth2spnormdensity_f(i,j,k)=0.0d0
            end if

           end do
          end do
        end do


        return
        end


c-----------------------------------------------------------------------
c Compute square root of density of the first energy functional, E_1[f], for f=(1-rho**2)**2*f1.
c E_1[f] is defined in eq. 27 of arxiv:1110.6794v2
c-----------------------------------------------------------------------
        subroutine sqrten1density_func(
     &                  sqrten1density_f,
     &                  f1_np1,f1_n,f1_nm1,
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
     &                  AH_R,AH_xc,
     &                  AH_Ntheta,AH_Nphi,
     &                  x,y,z,dt,ct,chr,L,ex,Nx,Ny,Nz,
     &                  phys_bdy,ghost_width,
     &                  ief_bh_r0,a_rot,kerrads_background)

        implicit none

        logical is_nan
        real*8  ief_bh_r0,a_rot,M0,M0_min
        real*8 kerrads_background
        integer Nx,Ny,Nz
        integer i,j,k
        integer phys_bdy(6),ghost_width(6)
        real*8 chr(Nx,Ny,Nz),ex
        real*8 x(Nx),y(Ny),z(Nz),dt,ct,L
        real*8 lambda4
        integer AH_Ntheta,AH_Nphi
        real*8 AH_R(AH_Ntheta,AH_Nphi)
        real*8 AH_xc(3)
        real*8 xp,yp,zp,rhop,thetap,phip
        integer i0,j0
        real*8 dahtheta,dahphi
        real*8 ft,fp
        real*8 AH_Rp,rhostop
        real*8 f1_np1(Nx,Ny,Nz),f1_n(Nx,Ny,Nz),f1_nm1(Nx,Ny,Nz)
        real*8 sqrten1density_f(Nx,Ny,Nz)
        real*8 en1dens_f0
        integer sp,hnorm_argtype
        real*8 sqrth10normdensity_f(Nx,Ny,Nz)
        real*8 sqrth0m2normdensity_dfdt(Nx,Ny,Nz)

        logical ltrace
        parameter (ltrace=.false.)


        integer is,ie,js,je,ks,ke

        integer i1,j1,k1,a,b,c,d,e,f,g,h,p,q,r
        integer ic,jc,kc
        real*8 efe_ires(4,4)

        real*8 PI
        parameter (PI=3.141592653589793d0)

        real*8 dx,dy,dz
        real*8 x0,y0,z0
        real*8 rho0,theta0,phi0   
        real*8 Rad0
        real*8 drho_dRad,dRad_drho

        real*8 gb_tt_np1(Nx,Ny,Nz),gb_tt_n(Nx,Ny,Nz),gb_tt_nm1(Nx,Ny,Nz)
        real*8 gb_tx_np1(Nx,Ny,Nz),gb_tx_n(Nx,Ny,Nz),gb_tx_nm1(Nx,Ny,Nz)
        real*8 gb_ty_np1(Nx,Ny,Nz),gb_ty_n(Nx,Ny,Nz),gb_ty_nm1(Nx,Ny,Nz)
        real*8 gb_tz_np1(Nx,Ny,Nz),gb_tz_n(Nx,Ny,Nz),gb_tz_nm1(Nx,Ny,Nz)
        real*8 gb_xx_np1(Nx,Ny,Nz),gb_xx_n(Nx,Ny,Nz),gb_xx_nm1(Nx,Ny,Nz)
        real*8 gb_xy_np1(Nx,Ny,Nz),gb_xy_n(Nx,Ny,Nz),gb_xy_nm1(Nx,Ny,Nz)
        real*8 gb_xz_np1(Nx,Ny,Nz),gb_xz_n(Nx,Ny,Nz),gb_xz_nm1(Nx,Ny,Nz)
        real*8 gb_yy_np1(Nx,Ny,Nz),gb_yy_n(Nx,Ny,Nz),gb_yy_nm1(Nx,Ny,Nz)
        real*8 gb_yz_np1(Nx,Ny,Nz),gb_yz_n(Nx,Ny,Nz),gb_yz_nm1(Nx,Ny,Nz)
        real*8 gb_zz_np1(Nx,Ny,Nz),gb_zz_n(Nx,Ny,Nz),gb_zz_nm1(Nx,Ny,Nz)
        real*8 Hb_t_np1(Nx,Ny,Nz),Hb_t_n(Nx,Ny,Nz),Hb_t_nm1(Nx,Ny,Nz)
        real*8 Hb_x_np1(Nx,Ny,Nz),Hb_x_n(Nx,Ny,Nz),Hb_x_nm1(Nx,Ny,Nz)
        real*8 Hb_y_np1(Nx,Ny,Nz),Hb_y_n(Nx,Ny,Nz),Hb_y_nm1(Nx,Ny,Nz)
        real*8 Hb_z_np1(Nx,Ny,Nz),Hb_z_n(Nx,Ny,Nz),Hb_z_nm1(Nx,Ny,Nz)
        real*8 phi1_np1(Nx,Ny,Nz),phi1_n(Nx,Ny,Nz),phi1_nm1(Nx,Ny,Nz)



        ! initialize fixed-size variables
        data i,j,k,is,ie,js,je,ks,ke/0,0,0,0,0,0,0,0,0/
        data ic,jc,kc/0,0,0/
        data i1,j1,k1,a,b,c,d,e,p,q,r/0,0,0,0,0,0,0,0,0,0,0/

        data dx,dy,dz/0.0,0.0,0.0/
        data x0,y0,z0,rho0/0.0,0.0,0.0,0.0/    

!----------------------------------------------------------------------


      ! Black hole mass
        M0=ief_bh_r0/2
      ! Minimum black hole mass. For M0 below this value, there is a naked singularity
        M0_min=((2*(1 + a_rot**2/L**2) + Sqrt((1 + a_rot**2/L**2)**2 
     &   + (12*a_rot**2)/L**2))*Sqrt(-1 + Sqrt((1 + a_rot**2/L**2)**2 
     &   + (12*a_rot**2)/L**2) - a_rot**2/L**2))/(3.*Sqrt(6.))

        if (a_rot.ge.L) then
         write (*,*) "sqrth0spnormdensity_phi: ERROR in choice of Kerr-AdS initial parameters: 
     &   the rotation parameter a must be smaller than the AdS radius L"
          write (*,*) "a_rot,L=",a_rot,L
          stop
        end if

        if ((abs(M0).gt.10.0d0**(-10))
     &     .and.(M0.le.M0_min)) then
          write (*,*) "sqrth0spnormdensity_phi: ERROR in choice of Kerr-AdS initial parameters: 
     &      the black hole mass M0=2*r0 must be larger
     &      than the M0_min value"
          write (*,*) "M0,M0_min=",M0,M0_min
          stop
        end if


        sp=0
        hnorm_argtype=0
        call sqrth1spnormdensity_func(
     &          sqrth10normdensity_f,
     &          sp,hnorm_argtype,
     &          f1_np1,f1_n,f1_nm1,
     &          gb_tt_np1,gb_tt_n,gb_tt_nm1,
     &          gb_tx_np1,gb_tx_n,gb_tx_nm1,
     &          gb_ty_np1,gb_ty_n,gb_ty_nm1,
     &          gb_tz_np1,gb_tz_n,gb_tz_nm1,
     &          gb_xx_np1,gb_xx_n,gb_xx_nm1,
     &          gb_xy_np1,gb_xy_n,gb_xy_nm1,
     &          gb_xz_np1,gb_xz_n,gb_xz_nm1,
     &          gb_yy_np1,gb_yy_n,gb_yy_nm1,
     &          gb_yz_np1,gb_yz_n,gb_yz_nm1,
     &          gb_zz_np1,gb_zz_n,gb_zz_nm1,
     &          Hb_t_np1,Hb_t_n,Hb_t_nm1,
     &          Hb_x_np1,Hb_x_n,Hb_x_nm1,
     &          Hb_y_np1,Hb_y_n,Hb_y_nm1,
     &          Hb_z_np1,Hb_z_n,Hb_z_nm1,
     &          phi1_np1,phi1_n,phi1_nm1,
     &          AH_R,AH_xc,
     &          AH_Ntheta,AH_Nphi,
     &          x,y,z,dt,ct,chr,L,ex,Nx,Ny,Nz,
     &          phys_bdy,ghost_width,
     &          ief_bh_r0,a_rot,kerrads_background)

        sp=-2
        hnorm_argtype=1
        call sqrth0spnormdensity_func(
     &          sqrth0m2normdensity_dfdt,
     &          sp,hnorm_argtype,
     &          f1_np1,f1_n,f1_nm1,
     &          AH_R,AH_xc,
     &          AH_Ntheta,AH_Nphi,
     &          x,y,z,dt,ct,chr,L,ex,Nx,Ny,Nz,
     &          phys_bdy,ghost_width,
     &          ief_bh_r0,a_rot,kerrads_background)



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
              theta0=acos(x0/rho0)
              if (z0.lt.0) then
               phi0=atan2(z0,y0)+2*PI
              else
               phi0=atan2(z0,y0)
              end if
              !uncompactified quasi-sherical Kerr-Schild radial coordinate
              Rad0=2*rho0/(1-rho0**2)
              drho_dRad=((1-rho0**2)**2)/(2*(1+rho0**2))
              dRad_drho=(-1 + rho0)**(-2) 
     &           + (1 + rho0)**(-2)

              !determine radius at which we stop computing the norms. This is the apparent horizon (AH) radius.

              !Cartesian coordinates of grid point p=(i,j,k) w.r.t. the centre of the AH
              xp=x0-AH_xc(1)
              yp=y0-AH_xc(2)
              zp=z0-AH_xc(3)

              !The center of the AH xp=yp=zp=0 will be excluded or included depending on the excision mask chr(i,j,k), 
              !so no need to force inclusion or exclusion via the value of rhostop. 
              !In fact, xp=yp=zp=0 is also the centre of the grid if there is no AH (also if the AH is centred at the origin), 
              !so we should include it if we want to compute these norms in the case of no horizon.
              !For this reason, it's best to leave the inclusion or exclusion of that point to the excision mask chr(i,j,k).
              if ((abs(xp).lt.10.0d0**(-12)).and.
     &            (abs(yp).lt.10.0d0**(-12)).and.  
     &            (abs(zp).lt.10.0d0**(-12)) ) then
                rhostop=-1
              else
                !corresponding spherical coordinates of point p
                !thetap,phip is in general NOT on the AH_Ntheta,AH_Nphi grid.
                rhop=sqrt(xp**2+yp**2+zp**2)
                thetap=acos(xp/rhop)
                if (zp.lt.0) then
                    phip=atan2(zp,yp)+2*PI
                else
                    phip=atan2(zp,yp)
                end if

                ! i0 is the index of the point with chi value closest to chip on the AH_Nchi grid (indices going from 0 to AH_Ntheta-1),
                ! and larger than chip (except if thetap=Pi, then it is smaller)
                ! similarly for j0
                dahtheta=PI/(AH_Ntheta-1)
                dahphi=2*PI/(AH_Nphi-1)
                i0=floor(thetap/dahtheta+1)
                if (AH_Ntheta-1<i0) i0=AH_Ntheta-1
                j0=floor(phip/dahphi+1)
                if (AH_Nphi-1<j0) j0=AH_Nphi-1
                ft=((thetap-(i0-1)*dahtheta))/dahtheta
                fp=((phip-(j0-1)*dahphi))/dahphi

                !bi-linear interpolated value of AH radius at (thetap,phip)
                AH_Rp=  AH_R(i0,j0)*(1-ft)*(1-fp)+
     &                  AH_R(i0+1,j0)*(ft)*(1-fp)+
     &                  AH_R(i0,j0+1)*(1-ft)*(fp)+
     &                  AH_R(i0+1,j0+1)*(ft)*(fp)

                if (AH_Rp.lt.1) then 
                    rhostop=AH_Rp
                else
                    rhostop=-1
                end if
              end if



            if ((rho0.gt.rhostop).and.
     &         (chr(i,j,k).ne.ex)) then

                en1dens_f0=sqrth10normdensity_f(i,j,k)**2
     &              +sqrth0m2normdensity_dfdt(i,j,k)**2


!Take square root, so we can use DV calculation of L2-norm, and then square it, to compute the H_AdS(1,s) norm.
                sqrten1density_f(i,j,k)=sqrt(en1dens_f0)

                 !the phi coordinate is not defined at y=z=0, i.e., theta=0,PI 
                !(not surprising since spherical coordinates are not valid everywhere on the sphere), 
                !hence we set the norm density to 0 at these points

                if ((abs(y0).lt.10.0d0**(-10)).and.
     &          (abs(z0).lt.10.0d0**(-10))) then
                    sqrten1density_f(i,j,k)=0.0d0
                end if

                if (ltrace) then
                 if ((abs(x0).gt.10.0d0**(-1.0d0)).and.
     &               (abs(y0).gt.10.0d0**(-1.0d0)).and.
     &               (abs(z0).gt.10.0d0**(-1.0d0)).and.
     &                (rho0.lt.0.9)  ) then
                   write(*,*) 'sqrten1density_func '
                   write(*,*) 'at i,j,k= ', i,j,k
                   write(*,*) 'i.e., x,y,z=', x(i),y(j),z(k)
                   write(*,*) ' sqrth10normdensity_f = ',
     &                  sqrth10normdensity_f(i,j,k)
                   write(*,*) ' sqrth0m2normdensity_dfdt = ',
     &                  sqrth0m2normdensity_dfdt(i,j,k)
                   write(*,*) ' sqrten1density_f = ',
     &                  sqrten1density_f(i,j,k)
                   return
                 end if
                end if


            else !i.e., the point is either excised or inside the apparent horizon
               sqrten1density_f(i,j,k)=0.0d0
            end if

           end do
          end do
        end do


        return
        end


c-----------------------------------------------------------------------
c Compute square root of density of the second energy functional, E_2[f], for f=(1-rho**2)**2*f1.
c E_2[f] is defined in eq. 27 of arxiv:1110.6794v2
c-----------------------------------------------------------------------
        subroutine sqrten2density_func(
     &                  sqrten2density_f,
     &                  f1_np1,f1_n,f1_nm1,
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
     &                  AH_R,AH_xc,
     &                  AH_Ntheta,AH_Nphi,
     &                  x,y,z,dt,ct,chr,L,ex,Nx,Ny,Nz,
     &                  phys_bdy,ghost_width,
     &                  ief_bh_r0,a_rot,kerrads_background)

        implicit none

        logical is_nan
        real*8  ief_bh_r0,a_rot,M0,M0_min
        real*8 kerrads_background
        integer Nx,Ny,Nz
        integer i,j,k
        integer phys_bdy(6),ghost_width(6)
        real*8 chr(Nx,Ny,Nz),ex
        real*8 x(Nx),y(Ny),z(Nz),dt,ct,L
        real*8 lambda4
        integer AH_Ntheta,AH_Nphi
        real*8 AH_R(AH_Ntheta,AH_Nphi)
        real*8 AH_xc(3)
        real*8 xp,yp,zp,rhop,thetap,phip
        integer i0,j0
        real*8 dahtheta,dahphi
        real*8 ft,fp
        real*8 AH_Rp,rhostop
        real*8 f1_np1(Nx,Ny,Nz),f1_n(Nx,Ny,Nz),f1_nm1(Nx,Ny,Nz)
        real*8 sqrten2density_f(Nx,Ny,Nz)
        real*8 en2dens_f0
        integer sp,hnorm_argtype
        real*8 sqrth20normdensity_f(Nx,Ny,Nz)
        real*8 sqrth10normdensity_dfdt(Nx,Ny,Nz)
        real*8 sqrth10normdensity_m1f(Nx,Ny,Nz)
        real*8 sqrth10normdensity_m2f(Nx,Ny,Nz)
        real*8 sqrth10normdensity_m3f(Nx,Ny,Nz)
        real*8 sqrth0m2normdensity_d2fdtdt(Nx,Ny,Nz)

        logical ltrace
        parameter (ltrace=.false.)


        integer is,ie,js,je,ks,ke

        integer i1,j1,k1,a,b,c,d,e,f,g,h,p,q,r
        integer ic,jc,kc
        real*8 efe_ires(4,4)

        real*8 PI
        parameter (PI=3.141592653589793d0)

        real*8 dx,dy,dz
        real*8 x0,y0,z0
        real*8 rho0,theta0,phi0   
        real*8 Rad0
        real*8 drho_dRad,dRad_drho

        real*8 gb_tt_np1(Nx,Ny,Nz),gb_tt_n(Nx,Ny,Nz),gb_tt_nm1(Nx,Ny,Nz)
        real*8 gb_tx_np1(Nx,Ny,Nz),gb_tx_n(Nx,Ny,Nz),gb_tx_nm1(Nx,Ny,Nz)
        real*8 gb_ty_np1(Nx,Ny,Nz),gb_ty_n(Nx,Ny,Nz),gb_ty_nm1(Nx,Ny,Nz)
        real*8 gb_tz_np1(Nx,Ny,Nz),gb_tz_n(Nx,Ny,Nz),gb_tz_nm1(Nx,Ny,Nz)
        real*8 gb_xx_np1(Nx,Ny,Nz),gb_xx_n(Nx,Ny,Nz),gb_xx_nm1(Nx,Ny,Nz)
        real*8 gb_xy_np1(Nx,Ny,Nz),gb_xy_n(Nx,Ny,Nz),gb_xy_nm1(Nx,Ny,Nz)
        real*8 gb_xz_np1(Nx,Ny,Nz),gb_xz_n(Nx,Ny,Nz),gb_xz_nm1(Nx,Ny,Nz)
        real*8 gb_yy_np1(Nx,Ny,Nz),gb_yy_n(Nx,Ny,Nz),gb_yy_nm1(Nx,Ny,Nz)
        real*8 gb_yz_np1(Nx,Ny,Nz),gb_yz_n(Nx,Ny,Nz),gb_yz_nm1(Nx,Ny,Nz)
        real*8 gb_zz_np1(Nx,Ny,Nz),gb_zz_n(Nx,Ny,Nz),gb_zz_nm1(Nx,Ny,Nz)
        real*8 Hb_t_np1(Nx,Ny,Nz),Hb_t_n(Nx,Ny,Nz),Hb_t_nm1(Nx,Ny,Nz)
        real*8 Hb_x_np1(Nx,Ny,Nz),Hb_x_n(Nx,Ny,Nz),Hb_x_nm1(Nx,Ny,Nz)
        real*8 Hb_y_np1(Nx,Ny,Nz),Hb_y_n(Nx,Ny,Nz),Hb_y_nm1(Nx,Ny,Nz)
        real*8 Hb_z_np1(Nx,Ny,Nz),Hb_z_n(Nx,Ny,Nz),Hb_z_nm1(Nx,Ny,Nz)
        real*8 phi1_np1(Nx,Ny,Nz),phi1_n(Nx,Ny,Nz),phi1_nm1(Nx,Ny,Nz)



        ! initialize fixed-size variables
        data i,j,k,is,ie,js,je,ks,ke/0,0,0,0,0,0,0,0,0/
        data ic,jc,kc/0,0,0/
        data i1,j1,k1,a,b,c,d,e,p,q,r/0,0,0,0,0,0,0,0,0,0,0/

        data dx,dy,dz/0.0,0.0,0.0/
        data x0,y0,z0,rho0/0.0,0.0,0.0,0.0/    

!----------------------------------------------------------------------


      ! Black hole mass
        M0=ief_bh_r0/2
      ! Minimum black hole mass. For M0 below this value, there is a naked singularity
        M0_min=((2*(1 + a_rot**2/L**2) + Sqrt((1 + a_rot**2/L**2)**2 
     &   + (12*a_rot**2)/L**2))*Sqrt(-1 + Sqrt((1 + a_rot**2/L**2)**2 
     &   + (12*a_rot**2)/L**2) - a_rot**2/L**2))/(3.*Sqrt(6.))

        if (a_rot.ge.L) then
         write (*,*) "sqrth0spnormdensity_phi: ERROR in choice of Kerr-AdS initial parameters: 
     &   the rotation parameter a must be smaller than the AdS radius L"
          write (*,*) "a_rot,L=",a_rot,L
          stop
        end if

        if ((abs(M0).gt.10.0d0**(-10))
     &     .and.(M0.le.M0_min)) then
          write (*,*) "sqrth0spnormdensity_phi: ERROR in choice of Kerr-AdS initial parameters: 
     &      the black hole mass M0=2*r0 must be larger
     &      than the M0_min value"
          write (*,*) "M0,M0_min=",M0,M0_min
          stop
        end if

        sp=0
        hnorm_argtype=0
        call sqrth2spnormdensity_func(
     &          sqrth20normdensity_f,
     &          sp,hnorm_argtype,
     &          f1_np1,f1_n,f1_nm1,
     &          gb_tt_np1,gb_tt_n,gb_tt_nm1,
     &          gb_tx_np1,gb_tx_n,gb_tx_nm1,
     &          gb_ty_np1,gb_ty_n,gb_ty_nm1,
     &          gb_tz_np1,gb_tz_n,gb_tz_nm1,
     &          gb_xx_np1,gb_xx_n,gb_xx_nm1,
     &          gb_xy_np1,gb_xy_n,gb_xy_nm1,
     &          gb_xz_np1,gb_xz_n,gb_xz_nm1,
     &          gb_yy_np1,gb_yy_n,gb_yy_nm1,
     &          gb_yz_np1,gb_yz_n,gb_yz_nm1,
     &          gb_zz_np1,gb_zz_n,gb_zz_nm1,
     &          Hb_t_np1,Hb_t_n,Hb_t_nm1,
     &          Hb_x_np1,Hb_x_n,Hb_x_nm1,
     &          Hb_y_np1,Hb_y_n,Hb_y_nm1,
     &          Hb_z_np1,Hb_z_n,Hb_z_nm1,
     &          phi1_np1,phi1_n,phi1_nm1,
     &          AH_R,AH_xc,
     &          AH_Ntheta,AH_Nphi,
     &          x,y,z,dt,ct,chr,L,ex,Nx,Ny,Nz,
     &          phys_bdy,ghost_width,
     &          ief_bh_r0,a_rot,kerrads_background)

        sp=0
        hnorm_argtype=1
        call sqrth1spnormdensity_func(
     &          sqrth10normdensity_dfdt,
     &          sp,hnorm_argtype,
     &          f1_np1,f1_n,f1_nm1,
     &          gb_tt_np1,gb_tt_n,gb_tt_nm1,
     &          gb_tx_np1,gb_tx_n,gb_tx_nm1,
     &          gb_ty_np1,gb_ty_n,gb_ty_nm1,
     &          gb_tz_np1,gb_tz_n,gb_tz_nm1,
     &          gb_xx_np1,gb_xx_n,gb_xx_nm1,
     &          gb_xy_np1,gb_xy_n,gb_xy_nm1,
     &          gb_xz_np1,gb_xz_n,gb_xz_nm1,
     &          gb_yy_np1,gb_yy_n,gb_yy_nm1,
     &          gb_yz_np1,gb_yz_n,gb_yz_nm1,
     &          gb_zz_np1,gb_zz_n,gb_zz_nm1,
     &          Hb_t_np1,Hb_t_n,Hb_t_nm1,
     &          Hb_x_np1,Hb_x_n,Hb_x_nm1,
     &          Hb_y_np1,Hb_y_n,Hb_y_nm1,
     &          Hb_z_np1,Hb_z_n,Hb_z_nm1,
     &          phi1_np1,phi1_n,phi1_nm1,
     &          AH_R,AH_xc,
     &          AH_Ntheta,AH_Nphi,
     &          x,y,z,dt,ct,chr,L,ex,Nx,Ny,Nz,
     &          phys_bdy,ghost_width,
     &          ief_bh_r0,a_rot,kerrads_background)

        sp=0
        hnorm_argtype=3
        call sqrth1spnormdensity_func(
     &          sqrth10normdensity_m1f,
     &          sp,hnorm_argtype,
     &          f1_np1,f1_n,f1_nm1,
     &          gb_tt_np1,gb_tt_n,gb_tt_nm1,
     &          gb_tx_np1,gb_tx_n,gb_tx_nm1,
     &          gb_ty_np1,gb_ty_n,gb_ty_nm1,
     &          gb_tz_np1,gb_tz_n,gb_tz_nm1,
     &          gb_xx_np1,gb_xx_n,gb_xx_nm1,
     &          gb_xy_np1,gb_xy_n,gb_xy_nm1,
     &          gb_xz_np1,gb_xz_n,gb_xz_nm1,
     &          gb_yy_np1,gb_yy_n,gb_yy_nm1,
     &          gb_yz_np1,gb_yz_n,gb_yz_nm1,
     &          gb_zz_np1,gb_zz_n,gb_zz_nm1,
     &          Hb_t_np1,Hb_t_n,Hb_t_nm1,
     &          Hb_x_np1,Hb_x_n,Hb_x_nm1,
     &          Hb_y_np1,Hb_y_n,Hb_y_nm1,
     &          Hb_z_np1,Hb_z_n,Hb_z_nm1,
     &          phi1_np1,phi1_n,phi1_nm1,
     &          AH_R,AH_xc,
     &          AH_Ntheta,AH_Nphi,
     &          x,y,z,dt,ct,chr,L,ex,Nx,Ny,Nz,
     &          phys_bdy,ghost_width,
     &          ief_bh_r0,a_rot,kerrads_background)

        sp=0
        hnorm_argtype=4
        call sqrth1spnormdensity_func(
     &          sqrth10normdensity_m2f,
     &          sp,hnorm_argtype,
     &          f1_np1,f1_n,f1_nm1,
     &          gb_tt_np1,gb_tt_n,gb_tt_nm1,
     &          gb_tx_np1,gb_tx_n,gb_tx_nm1,
     &          gb_ty_np1,gb_ty_n,gb_ty_nm1,
     &          gb_tz_np1,gb_tz_n,gb_tz_nm1,
     &          gb_xx_np1,gb_xx_n,gb_xx_nm1,
     &          gb_xy_np1,gb_xy_n,gb_xy_nm1,
     &          gb_xz_np1,gb_xz_n,gb_xz_nm1,
     &          gb_yy_np1,gb_yy_n,gb_yy_nm1,
     &          gb_yz_np1,gb_yz_n,gb_yz_nm1,
     &          gb_zz_np1,gb_zz_n,gb_zz_nm1,
     &          Hb_t_np1,Hb_t_n,Hb_t_nm1,
     &          Hb_x_np1,Hb_x_n,Hb_x_nm1,
     &          Hb_y_np1,Hb_y_n,Hb_y_nm1,
     &          Hb_z_np1,Hb_z_n,Hb_z_nm1,
     &          phi1_np1,phi1_n,phi1_nm1,
     &          AH_R,AH_xc,
     &          AH_Ntheta,AH_Nphi,
     &          x,y,z,dt,ct,chr,L,ex,Nx,Ny,Nz,
     &          phys_bdy,ghost_width,
     &          ief_bh_r0,a_rot,kerrads_background)

        sp=0
        hnorm_argtype=5
        call sqrth1spnormdensity_func(
     &          sqrth10normdensity_m3f,
     &          sp,hnorm_argtype,
     &          f1_np1,f1_n,f1_nm1,
     &          gb_tt_np1,gb_tt_n,gb_tt_nm1,
     &          gb_tx_np1,gb_tx_n,gb_tx_nm1,
     &          gb_ty_np1,gb_ty_n,gb_ty_nm1,
     &          gb_tz_np1,gb_tz_n,gb_tz_nm1,
     &          gb_xx_np1,gb_xx_n,gb_xx_nm1,
     &          gb_xy_np1,gb_xy_n,gb_xy_nm1,
     &          gb_xz_np1,gb_xz_n,gb_xz_nm1,
     &          gb_yy_np1,gb_yy_n,gb_yy_nm1,
     &          gb_yz_np1,gb_yz_n,gb_yz_nm1,
     &          gb_zz_np1,gb_zz_n,gb_zz_nm1,
     &          Hb_t_np1,Hb_t_n,Hb_t_nm1,
     &          Hb_x_np1,Hb_x_n,Hb_x_nm1,
     &          Hb_y_np1,Hb_y_n,Hb_y_nm1,
     &          Hb_z_np1,Hb_z_n,Hb_z_nm1,
     &          phi1_np1,phi1_n,phi1_nm1,
     &          AH_R,AH_xc,
     &          AH_Ntheta,AH_Nphi,
     &          x,y,z,dt,ct,chr,L,ex,Nx,Ny,Nz,
     &          phys_bdy,ghost_width,
     &          ief_bh_r0,a_rot,kerrads_background)

        sp=-2
        hnorm_argtype=2
        call sqrth0spnormdensity_func(
     &          sqrth0m2normdensity_d2fdtdt,
     &          sp,hnorm_argtype,
     &          f1_np1,f1_n,f1_nm1,
     &          AH_R,AH_xc,
     &          AH_Ntheta,AH_Nphi,
     &          x,y,z,dt,ct,chr,L,ex,Nx,Ny,Nz,
     &          phys_bdy,ghost_width,
     &          ief_bh_r0,a_rot,kerrads_background)


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
              theta0=acos(x0/rho0)
              if (z0.lt.0) then
               phi0=atan2(z0,y0)+2*PI
              else
               phi0=atan2(z0,y0)
              end if
              !uncompactified quasi-sherical Kerr-Schild radial coordinate
              Rad0=2*rho0/(1-rho0**2)
              drho_dRad=((1-rho0**2)**2)/(2*(1+rho0**2))
              dRad_drho=(-1 + rho0)**(-2) 
     &           + (1 + rho0)**(-2)

              !determine radius at which we stop computing the norms. This is the apparent horizon (AH) radius.

              !Cartesian coordinates of grid point p=(i,j,k) w.r.t. the centre of the AH
              xp=x0-AH_xc(1)
              yp=y0-AH_xc(2)
              zp=z0-AH_xc(3)

              !The center of the AH xp=yp=zp=0 will be excluded or included depending on the excision mask chr(i,j,k), 
              !so no need to force inclusion or exclusion via the value of rhostop. 
              !In fact, xp=yp=zp=0 is also the centre of the grid if there is no AH (also if the AH is centred at the origin), 
              !so we should include it if we want to compute these norms in the case of no horizon.
              !For this reason, it's best to leave the inclusion or exclusion of that point to the excision mask chr(i,j,k).
              if ((abs(xp).lt.10.0d0**(-12)).and.
     &            (abs(yp).lt.10.0d0**(-12)).and.  
     &            (abs(zp).lt.10.0d0**(-12)) ) then
                rhostop=-1
              else
                !corresponding spherical coordinates of point p
                !thetap,phip is in general NOT on the AH_Ntheta,AH_Nphi grid.
                rhop=sqrt(xp**2+yp**2+zp**2)
                thetap=acos(xp/rhop)
                if (zp.lt.0) then
                    phip=atan2(zp,yp)+2*PI
                else
                    phip=atan2(zp,yp)
                end if

                ! i0 is the index of the point with chi value closest to chip on the AH_Nchi grid (indices going from 0 to AH_Ntheta-1),
                ! and larger than chip (except if thetap=Pi, then it is smaller)
                ! similarly for j0
                dahtheta=PI/(AH_Ntheta-1)
                dahphi=2*PI/(AH_Nphi-1)
                i0=floor(thetap/dahtheta+1)
                if (AH_Ntheta-1<i0) i0=AH_Ntheta-1
                j0=floor(phip/dahphi+1)
                if (AH_Nphi-1<j0) j0=AH_Nphi-1
                ft=((thetap-(i0-1)*dahtheta))/dahtheta
                fp=((phip-(j0-1)*dahphi))/dahphi

                !bi-linear interpolated value of AH radius at (thetap,phip)
                AH_Rp=  AH_R(i0,j0)*(1-ft)*(1-fp)+
     &                  AH_R(i0+1,j0)*(ft)*(1-fp)+
     &                  AH_R(i0,j0+1)*(1-ft)*(fp)+
     &                  AH_R(i0+1,j0+1)*(ft)*(fp)

                if (AH_Rp.lt.1) then 
                    rhostop=AH_Rp
                else
                    rhostop=-1
                end if
              end if



            if ((rho0.gt.rhostop).and.
     &         (chr(i,j,k).ne.ex)) then

                en2dens_f0=sqrth20normdensity_f(i,j,k)**2
     &              +sqrth10normdensity_dfdt(i,j,k)**2
     &              +sqrth10normdensity_m1f(i,j,k)**2
     &              +sqrth10normdensity_m2f(i,j,k)**2
     &              +sqrth10normdensity_m3f(i,j,k)**2
     &              +sqrth0m2normdensity_d2fdtdt(i,j,k)**2


!Take square root, so we can use DV calculation of L2-norm, and then square it, to compute the H_AdS(1,s) norm.
                sqrten2density_f(i,j,k)=sqrt(en2dens_f0)

                 !the phi coordinate is not defined at y=z=0, i.e., theta=0,PI 
                !(not surprising since spherical coordinates are not valid everywhere on the sphere), 
                !hence we set the norm density to 0 at these points

                if ((abs(y0).lt.10.0d0**(-10)).and.
     &          (abs(z0).lt.10.0d0**(-10))) then
                    sqrten2density_f(i,j,k)=0.0d0
                end if

                if (ltrace) then
                 if ((abs(x0).gt.10.0d0**(-1.0d0)).and.
     &               (abs(y0).gt.10.0d0**(-1.0d0)).and.
     &               (abs(z0).gt.10.0d0**(-1.0d0)).and.
     &                (rho0.lt.0.9)  ) then
                   write(*,*) 'sqrten2density_func '
                   write(*,*) 'at i,j,k= ', i,j,k
                   write(*,*) 'i.e., x,y,z=', x(i),y(j),z(k)
                   write(*,*) ' sqrth20normdensity_f = ',
     &                  sqrth20normdensity_f(i,j,k)
                   write(*,*) ' sqrth10normdensity_dfdt = ',
     &                  sqrth10normdensity_dfdt(i,j,k)
                   write(*,*) ' sqrth10normdensity_m1f = ',
     &                  sqrth10normdensity_m1f(i,j,k)
                   write(*,*) ' sqrth10normdensity_m2f = ',
     &                  sqrth10normdensity_m2f(i,j,k)
                   write(*,*) ' sqrth10normdensity_m3f = ',
     &                  sqrth10normdensity_m3f(i,j,k)
                   write(*,*) ' sqrth0m2normdensity_d2fdtdt = ',
     &                  sqrth0m2normdensity_d2fdtdt(i,j,k)
                   write(*,*) ' sqrten2density_f = ',
     &                  sqrten2density_f(i,j,k)
                   !return
                   stop
                 end if
                end if


            else !i.e., the point is either excised or inside the apparent horizon
               sqrten2density_f(i,j,k)=0.0d0
            end if

           end do
          end do
        end do


        return
        end

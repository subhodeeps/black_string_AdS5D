c----------------------------------------------------------------------
c add a small Lorentz-boosted Gaussian perturbation to the initial scalar field profile
c and gbs
c centre of the perturbation is (boost_xu0,boost_yu0,boost_zu0), boost velocity is (boost_vx,boost_vy,boost_vz)
c----------------------------------------------------------------------

        subroutine boost_perturb(
     &                     phi1_np1,phi1_n,phi1_nm1,phi1_t_n,
     &					      gb_tt_np1,gb_tt_n,gb_tt_nm1,gb_tt_t_n,
     &                     gb_tx_np1,gb_tx_n,gb_tx_nm1,gb_tx_t_n,
     &                     gb_ty_np1,gb_ty_n,gb_ty_nm1,gb_ty_t_n,
     &                     gb_tz_np1,gb_tz_n,gb_tz_nm1,gb_tz_t_n,
     &                     gb_xx_np1,gb_xx_n,gb_xx_nm1,gb_xx_t_n,
     &                     gb_xy_np1,gb_xy_n,gb_xy_nm1,gb_xy_t_n,
     &                     gb_xz_np1,gb_xz_n,gb_xz_nm1,gb_xz_t_n,
     &                     gb_yy_np1,gb_yy_n,gb_yy_nm1,gb_yy_t_n,
     &                     gb_yz_np1,gb_yz_n,gb_yz_nm1,gb_yz_t_n,
     &                     gb_zz_np1,gb_zz_n,gb_zz_nm1,gb_zz_t_n,
     &                     boost_vx,boost_vy,boost_vz,
     &                     boost_amp,
     &                     boost_r0,boost_delta,
     &                     boost_xu0,boost_yu0,boost_zu0,
     &                     boost_ex,boost_ey,boost_ez,
     &                     L,x,y,z,dt,chr,exc,Nx,Ny,Nz)

        implicit none
        integer Nx,Ny,Nz
        real*8 f_n(Nx,Ny,Nz),f_t_n(Nx,Ny,Nz)
        real*8 f_np1(Nx,Ny,Nz),f_nm1(Nx,Ny,Nz)
        real*8 phi1_n(Nx,Ny,Nz),phi1_t_n(Nx,Ny,Nz)
        real*8 phi1_np1(Nx,Ny,Nz),phi1_nm1(Nx,Ny,Nz)
        real*8 x(Nx),y(Ny),z(Nz),dt
        real*8 boost_amp,boost_r0,boost_delta
        real*8 boost_ex,boost_ey,boost_ez
        real*8 boost_xu0,boost_yu0,boost_zu0
        real*8 L
        real*8 chr(Nx,Ny,Nz),exc

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
        real*8 gb_tt_t_n(Nx,Ny,Nz)
        real*8 gb_tx_t_n(Nx,Ny,Nz)
        real*8 gb_ty_t_n(Nx,Ny,Nz)
        real*8 gb_tz_t_n(Nx,Ny,Nz)
        real*8 gb_xx_t_n(Nx,Ny,Nz)
        real*8 gb_xy_t_n(Nx,Ny,Nz)
        real*8 gb_xz_t_n(Nx,Ny,Nz)
        real*8 gb_yy_t_n(Nx,Ny,Nz)
        real*8 gb_yz_t_n(Nx,Ny,Nz)
        real*8 gb_zz_t_n(Nx,Ny,Nz)

        logical is_nan

        integer i,j,k
        integer stype
        integer a,b
        real*8 r,x0,y0,z0,rho0,xi0,chi0,csr,xb,yb,zb
        real*8 boost_vx,boost_vy,boost_vz
        real*8 boost_vnorm
        real*8 gamma
        real*8 lambda_boost(4,4),lambdainv_boost(4,4)
        real*8 t0,t0_invboost
        real*8 x0_invboost,y0_invboost,z0_invboost
        real*8 xb_invboost,yb_invboost,zb_invboost
        real*8 r_invboost
        real*8 x0_invboost_t,y0_invboost_t,z0_invboost_t
        real*8 r_invboost_t
        real*8 derexp0

        real*8 rhoc,rhod
        real*8 f1,trans

        real*8 PI
        parameter (PI=3.141592653589793d0)

        ! initialize fixed-size variables
        data i,j,k/0,0,0/
        data r,x0,y0,z0,rho0,xi0,chi0,csr,xb,yb,zb
     &       /0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/
        data lambda_boost,lambdainv_boost/16*0.0,16*0.0/

        !--------------------------------------------------------------

       if ((abs(boost_vx).gt.10.0d0**(-10)).or.
     &     (abs(boost_vy).gt.10.0d0**(-10)).or.
     &     (abs(boost_vz).gt.10.0d0**(-10)) ) then
        do i=1,Nx
           do j=1,Ny
            do k=1,Nz
             if (chr(i,j,k).ne.exc) then
              f_n(i,j,k)=0
              f_np1(i,j,k)=0
              f_nm1(i,j,k)=0
              f_t_n(i,j,k)=0
              x0=x(i)
              y0=y(j)
              z0=z(k)
              rho0=sqrt(x0**2+y0**2+z0**2)
              if (rho0.ne.0.0d0) then
               xi0=acos(x0/rho0)
              end if
              if ((y0.ne.0.0d0).or.(z0.ne.0.0d0)) then
                chi0=atan2(z0,y0)
                if (chi0.lt.0) chi0=chi0+2*PI
              end if
!            write (*,*) 'DEBUG from gauss3d'
!            write(*,*) 'rho0=',rho0
!            write(*,*) 'chi0,xi0=',chi0,xi0


                boost_vnorm=sqrt(boost_vx**2+boost_vy**2+boost_vz**2)
                gamma=1/sqrt(1-boost_vnorm**2)

                lambda_boost(1,1)=gamma
                lambda_boost(1,2)=-gamma*(-boost_vx)
                lambda_boost(1,3)=-gamma*(-boost_vy)
                lambda_boost(1,4)=-gamma*(-boost_vz)
                lambda_boost(2,2)=1+(gamma-1)*((-boost_vx)**2)
     &           /boost_vnorm**2
                lambda_boost(2,3)=(gamma-1)*(-boost_vx)*(-boost_vy)
     &           /boost_vnorm**2
                lambda_boost(2,4)=(gamma-1)*(-boost_vx)*(-boost_vz)
     &           /boost_vnorm**2
                lambda_boost(3,3)=1+(gamma-1)*((-boost_vy)**2)
     &           /boost_vnorm**2
                lambda_boost(3,4)=(gamma-1)*(-boost_vy)*(-boost_vz)
     &           /boost_vnorm**2
                lambda_boost(4,4)=1+(gamma-1)*((-boost_vz)**2)
     &           /boost_vnorm**2
                !the matrix of Lorentz boosts is symmetric
                do a=1,4
                 do b=a+1,4
                  lambda_boost(b,a)=lambda_boost(a,b)
                 end do
                end do

                lambdainv_boost(1,1)=-1/((-1+boost_vnorm**2)*gamma)
                lambdainv_boost(1,2)=-(-boost_vx)
     &           /((-1+boost_vnorm**2)*gamma)
                lambdainv_boost(1,3)=-(-boost_vy)
     &           /((-1+boost_vnorm**2)*gamma)
                lambdainv_boost(1,4)=-(-boost_vz)
     &           /((-1+boost_vnorm**2)*gamma)
                lambdainv_boost(2,2)=(1/boost_vnorm**2)*
     &           ((-boost_vy)**2+(-boost_vz)**2-((-boost_vx)**2)
     &           /((-1+boost_vnorm**2)*gamma))
                lambdainv_boost(2,3)=-((-boost_vx)*(-boost_vy)*
     &           (1+(-1+boost_vnorm**2)*gamma))
     &           /((-1+boost_vnorm**2)*boost_vnorm**2*gamma)
                lambdainv_boost(2,4)=-((-boost_vx)*(-boost_vz)*
     &           (1+(-1+boost_vnorm**2)*gamma))
     &           /((-1+boost_vnorm**2)*boost_vnorm**2*gamma)
                lambdainv_boost(3,3)=(1/boost_vnorm**2)*
     &           ((-boost_vx)**2+(-boost_vz)**2-((-boost_vy)**2)
     &           /((-1+boost_vnorm**2)*gamma))
                lambdainv_boost(3,4)=-((-boost_vy)*(-boost_vz)*
     &           (1+(-1+boost_vnorm**2)*gamma))
     &           /((-1+boost_vnorm**2)*boost_vnorm**2*gamma)
                lambdainv_boost(4,4)=(1/boost_vnorm**2)*
     &           ((-boost_vx)**2+(-boost_vy)**2-((-boost_vz)**2)
     &           /((-1+boost_vnorm**2)*gamma))
                !the matrix of Lorentz boosts is symmetric
                do a=1,4
                 do b=a+1,4
                  lambdainv_boost(b,a)=lambdainv_boost(a,b)
                 end do
                end do

                !The inverse Lorentz boost is lambdainv_boost^(-1)^\mu_\nu x^\mu in UNCOMPACTIFIED CARTESIAN COORDINATES
                !We compute this in compactified Cartesian coords at t=0
                t0=0
                !computing t0_invboost is not necessary but we do it for completeness

                t0_invboost=
     -  lambdainv_boost(1,1)*t0 - 
     -  (2*(lambdainv_boost(1,2)*x0 + lambdainv_boost(1,3)*y0 + 
     -       lambdainv_boost(1,4)*z0))/
     -   (-1 + x0**2 + y0**2 + z0**2)
                x0_invboost=
     -         ((lambdainv_boost(1,2)*t0 - 
     -      (2*(lambdainv_boost(2,2)*x0 + lambdainv_boost(2,3)*y0 + 
     -           lambdainv_boost(2,4)*z0))/
     -       (-1 + x0**2 + y0**2 + z0**2))*
     -    (-1 + Sqrt(1 + 
     -        (lambdainv_boost(1,2)*t0 - 
     -           (2*(lambdainv_boost(2,2)*x0 + 
     -                lambdainv_boost(2,3)*y0 + 
     -                lambdainv_boost(2,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2
     -         + (lambdainv_boost(1,3)*t0 - 
     -           (2*(lambdainv_boost(2,3)*x0 + 
     -                lambdainv_boost(3,3)*y0 + 
     -                lambdainv_boost(3,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2
     -         + (lambdainv_boost(1,4)*t0 - 
     -           (2*(lambdainv_boost(2,4)*x0 + 
     -                lambdainv_boost(3,4)*y0 + 
     -                lambdainv_boost(4,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2)
     -      ))/
     -  ((lambdainv_boost(1,2)*t0 - 
     -       (2*(lambdainv_boost(2,2)*x0 + lambdainv_boost(2,3)*y0 + 
     -            lambdainv_boost(2,4)*z0))/
     -        (-1 + x0**2 + y0**2 + z0**2))**2 + 
     -    (lambdainv_boost(1,3)*t0 - 
     -       (2*(lambdainv_boost(2,3)*x0 + lambdainv_boost(3,3)*y0 + 
     -            lambdainv_boost(3,4)*z0))/
     -        (-1 + x0**2 + y0**2 + z0**2))**2 + 
     -    (lambdainv_boost(1,4)*t0 - 
     -       (2*(lambdainv_boost(2,4)*x0 + lambdainv_boost(3,4)*y0 + 
     -            lambdainv_boost(4,4)*z0))/
     -        (-1 + x0**2 + y0**2 + z0**2))**2)
                y0_invboost=
     -          ((-2*(lambdainv_boost(2,3)*x0 
     -      + lambdainv_boost(3,3)*y0 + 
     -         lambdainv_boost(3,4)*z0) + 
     -      lambdainv_boost(1,3)*t0*
     -       (-1 + x0**2 + y0**2 + z0**2))*
     -    Sqrt(1 - (lambdainv_boost(1,2)*t0 - 
     -          (2*(lambdainv_boost(2,2)*x0 + 
     -               lambdainv_boost(2,3)*y0 + 
     -               lambdainv_boost(2,4)*z0))/
     -           (-1 + x0**2 + y0**2 + z0**2))**2/
     -       ((lambdainv_boost(1,2)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,2)*x0 + 
     -                 lambdainv_boost(2,3)*y0 + 
     -                 lambdainv_boost(2,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -           + (lambdainv_boost(1,3)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,3)*x0 + 
     -                 lambdainv_boost(3,3)*y0 + 
     -                 lambdainv_boost(3,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -           + (lambdainv_boost(1,4)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,4)*x0 + 
     -                 lambdainv_boost(3,4)*y0 + 
     -                 lambdainv_boost(4,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -         ))*(-1 + 
     -      Sqrt(1 + 
     -        (lambdainv_boost(1,2)*t0 - 
     -           (2*(lambdainv_boost(2,2)*x0 + 
     -                lambdainv_boost(2,3)*y0 + 
     -                lambdainv_boost(2,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2
     -         + (lambdainv_boost(1,3)*t0 - 
     -           (2*(lambdainv_boost(2,3)*x0 + 
     -                lambdainv_boost(3,3)*y0 + 
     -                lambdainv_boost(3,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2
     -         + (lambdainv_boost(1,4)*t0 - 
     -           (2*(lambdainv_boost(2,4)*x0 + 
     -                lambdainv_boost(3,4)*y0 + 
     -                lambdainv_boost(4,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2)
     -      ))/
     -  ((-1 + x0**2 + y0**2 + z0**2)*
     -    Sqrt((lambdainv_boost(1,3)*t0 - 
     -         (2*(lambdainv_boost(2,3)*x0 + 
     -              lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -              ))/(-1 + x0**2 + y0**2 + z0**2)
     -         )**2 + 
     -      (lambdainv_boost(1,4)*t0 - 
     -         (2*(lambdainv_boost(2,4)*x0 + 
     -              lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -              ))/(-1 + x0**2 + y0**2 + z0**2)
     -         )**2)*
     -    Sqrt((lambdainv_boost(1,2)*t0 - 
     -         (2*(lambdainv_boost(2,2)*x0 + 
     -              lambdainv_boost(2,3)*y0 + lambdainv_boost(2,4)*z0
     -              ))/(-1 + x0**2 + y0**2 + z0**2)
     -         )**2 + 
     -      (lambdainv_boost(1,3)*t0 - 
     -         (2*(lambdainv_boost(2,3)*x0 + 
     -              lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -              ))/(-1 + x0**2 + y0**2 + z0**2)
     -         )**2 + 
     -      (lambdainv_boost(1,4)*t0 - 
     -         (2*(lambdainv_boost(2,4)*x0 + 
     -              lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -              ))/(-1 + x0**2 + y0**2 + z0**2)
     -         )**2))
                z0_invboost=
     -  ((-2*(lambdainv_boost(2,4)*x0 + lambdainv_boost(3,4)*y0 + 
     -         lambdainv_boost(4,4)*z0) + 
     -      lambdainv_boost(1,4)*t0*
     -       (-1 + x0**2 + y0**2 + z0**2))*
     -    Sqrt(1 - (lambdainv_boost(1,2)*t0 - 
     -          (2*(lambdainv_boost(2,2)*x0 + 
     -               lambdainv_boost(2,3)*y0 + 
     -               lambdainv_boost(2,4)*z0))/
     -           (-1 + x0**2 + y0**2 + z0**2))**2/
     -       ((lambdainv_boost(1,2)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,2)*x0 + 
     -                 lambdainv_boost(2,3)*y0 + 
     -                 lambdainv_boost(2,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -           + (lambdainv_boost(1,3)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,3)*x0 + 
     -                 lambdainv_boost(3,3)*y0 + 
     -                 lambdainv_boost(3,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -           + (lambdainv_boost(1,4)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,4)*x0 + 
     -                 lambdainv_boost(3,4)*y0 + 
     -                 lambdainv_boost(4,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -         ))*(-1 + 
     -      Sqrt(1 + 
     -        (lambdainv_boost(1,2)*t0 - 
     -           (2*(lambdainv_boost(2,2)*x0 + 
     -                lambdainv_boost(2,3)*y0 + 
     -                lambdainv_boost(2,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2
     -         + (lambdainv_boost(1,3)*t0 - 
     -           (2*(lambdainv_boost(2,3)*x0 + 
     -                lambdainv_boost(3,3)*y0 + 
     -                lambdainv_boost(3,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2
     -         + (lambdainv_boost(1,4)*t0 - 
     -           (2*(lambdainv_boost(2,4)*x0 + 
     -                lambdainv_boost(3,4)*y0 + 
     -                lambdainv_boost(4,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2)
     -      ))/
     -  ((-1 + x0**2 + y0**2 + z0**2)*
     -    Sqrt((lambdainv_boost(1,3)*t0 - 
     -         (2*(lambdainv_boost(2,3)*x0 + 
     -              lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -              ))/(-1 + x0**2 + y0**2 + z0**2)
     -         )**2 + 
     -      (lambdainv_boost(1,4)*t0 - 
     -         (2*(lambdainv_boost(2,4)*x0 + 
     -              lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -              ))/(-1 + x0**2 + y0**2 + z0**2)
     -         )**2)*
     -    Sqrt((lambdainv_boost(1,2)*t0 - 
     -         (2*(lambdainv_boost(2,2)*x0 + 
     -              lambdainv_boost(2,3)*y0 + lambdainv_boost(2,4)*z0
     -              ))/(-1 + x0**2 + y0**2 + z0**2)
     -         )**2 + 
     -      (lambdainv_boost(1,3)*t0 - 
     -         (2*(lambdainv_boost(2,3)*x0 + 
     -              lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -              ))/(-1 + x0**2 + y0**2 + z0**2)
     -         )**2 + 
     -      (lambdainv_boost(1,4)*t0 - 
     -         (2*(lambdainv_boost(2,4)*x0 + 
     -              lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -              ))/(-1 + x0**2 + y0**2 + z0**2)
     -         )**2))

              if ((abs(y0).lt.10.0d0**(-10))
     -         .and.(abs(z0).lt.10.0d0**(-10)) ) then
                y0_invboost=0
                z0_invboost=0
              end if

              if ((abs(x0).lt.10.0d0**(-10))
     -         .and.(abs(z0).lt.10.0d0**(-10)) ) then
                x0_invboost=0
                z0_invboost=0
              end if

              if ((abs(x0).lt.10.0d0**(-10))
     -         .and.(abs(y0).lt.10.0d0**(-10)) ) then
                x0_invboost=0
                y0_invboost=0
              end if

              if ((abs(x0).lt.10.0d0**(-10))
     -         .and.(abs(y0).lt.10.0d0**(-10))
     -          .and.(abs(z0).lt.10.0d0**(-10)) ) then
                x0_invboost=0
                y0_invboost=0
                z0_invboost=0
              end if




               xb_invboost=x0_invboost-boost_xu0
               yb_invboost=y0_invboost-boost_yu0
               zb_invboost=z0_invboost-boost_zu0

               r_invboost=sqrt(xb_invboost**2+yb_invboost**2
     &          +zb_invboost**2
     &           -boost_ex**2*xb_invboost**2-boost_ey**2*yb_invboost**2
     &           -boost_ez**2*zb_invboost**2)

     			f_n(i,j,k)=boost_amp*
     &             exp(-((r_invboost-boost_r0)/boost_delta)**2)

               phi1_n(i,j,k)=phi1_n(i,j,k)+f_n(i,j,k)
              
               gb_tt_n(i,j,k)=gb_tt_n(i,j,k)+f_n(i,j,k)
               gb_tx_n(i,j,k)=gb_tx_n(i,j,k)+f_n(i,j,k)
               gb_ty_n(i,j,k)=gb_ty_n(i,j,k)+f_n(i,j,k)
               gb_tz_n(i,j,k)=gb_tz_n(i,j,k)+f_n(i,j,k)
               gb_xx_n(i,j,k)=gb_xx_n(i,j,k)+f_n(i,j,k)
               gb_xy_n(i,j,k)=gb_xy_n(i,j,k)+f_n(i,j,k)
               gb_xz_n(i,j,k)=gb_xz_n(i,j,k)+f_n(i,j,k)
               gb_yy_n(i,j,k)=gb_yy_n(i,j,k)+f_n(i,j,k)
               gb_yz_n(i,j,k)=gb_yz_n(i,j,k)+f_n(i,j,k)
               gb_zz_n(i,j,k)=gb_zz_n(i,j,k)+f_n(i,j,k)















               !set derivative of f(i,j,k) w.r.t. t at t=0

               x0_invboost_t=
     -         ((lambdainv_boost(1,2)*t0 - 
     -       (2*(lambdainv_boost(2,2)*x0 + lambdainv_boost(2,3)*y0 + 
     -            lambdainv_boost(2,4)*z0))/
     -        (-1 + x0**2 + y0**2 + z0**2))*
     -     (2*lambdainv_boost(1,2)*
     -        (lambdainv_boost(1,2)*t0 - 
     -          (2*(lambdainv_boost(2,2)*x0 + 
     -               lambdainv_boost(2,3)*y0 + lambdainv_boost(2,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          ) + 2*lambdainv_boost(1,3)*
     -        (lambdainv_boost(1,3)*t0 - 
     -          (2*(lambdainv_boost(2,3)*x0 + 
     -               lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          ) + 2*lambdainv_boost(1,4)*
     -        (lambdainv_boost(1,4)*t0 - 
     -          (2*(lambdainv_boost(2,4)*x0 + 
     -               lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )))/
     -   (2.*((lambdainv_boost(1,2)*t0 - 
     -          (2*(lambdainv_boost(2,2)*x0 + 
     -               lambdainv_boost(2,3)*y0 + lambdainv_boost(2,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2 + 
     -       (lambdainv_boost(1,3)*t0 - 
     -          (2*(lambdainv_boost(2,3)*x0 + 
     -               lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2 + 
     -       (lambdainv_boost(1,4)*t0 - 
     -          (2*(lambdainv_boost(2,4)*x0 + 
     -               lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2)*
     -     Sqrt(1 + (lambdainv_boost(1,2)*t0 - 
     -          (2*(lambdainv_boost(2,2)*x0 + 
     -               lambdainv_boost(2,3)*y0 + lambdainv_boost(2,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2 + 
     -       (lambdainv_boost(1,3)*t0 - 
     -          (2*(lambdainv_boost(2,3)*x0 + 
     -               lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2 + 
     -       (lambdainv_boost(1,4)*t0 - 
     -          (2*(lambdainv_boost(2,4)*x0 + 
     -               lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2)) - 
     -  ((lambdainv_boost(1,2)*t0 - 
     -       (2*(lambdainv_boost(2,2)*x0 + lambdainv_boost(2,3)*y0 + 
     -            lambdainv_boost(2,4)*z0))/
     -        (-1 + x0**2 + y0**2 + z0**2))*
     -     (2*lambdainv_boost(1,2)*
     -        (lambdainv_boost(1,2)*t0 - 
     -          (2*(lambdainv_boost(2,2)*x0 + 
     -               lambdainv_boost(2,3)*y0 + lambdainv_boost(2,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          ) + 2*lambdainv_boost(1,3)*
     -        (lambdainv_boost(1,3)*t0 - 
     -          (2*(lambdainv_boost(2,3)*x0 + 
     -               lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          ) + 2*lambdainv_boost(1,4)*
     -        (lambdainv_boost(1,4)*t0 - 
     -          (2*(lambdainv_boost(2,4)*x0 + 
     -               lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          ))*(-1 + 
     -       Sqrt(1 + 
     -         (lambdainv_boost(1,2)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,2)*x0 + 
     -                 lambdainv_boost(2,3)*y0 + 
     -                 lambdainv_boost(2,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -          + (lambdainv_boost(1,3)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,3)*x0 + 
     -                 lambdainv_boost(3,3)*y0 + 
     -                 lambdainv_boost(3,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -          + (lambdainv_boost(1,4)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,4)*x0 + 
     -                 lambdainv_boost(3,4)*y0 + 
     -                 lambdainv_boost(4,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2)
     -       ))/
     -   ((lambdainv_boost(1,2)*t0 - 
     -         (2*(lambdainv_boost(2,2)*x0 + 
     -              lambdainv_boost(2,3)*y0 + lambdainv_boost(2,4)*z0)
     -            )/(-1 + x0**2 + y0**2 + z0**2))**2
     -        + (lambdainv_boost(1,3)*t0 - 
     -         (2*(lambdainv_boost(2,3)*x0 + 
     -              lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0)
     -            )/(-1 + x0**2 + y0**2 + z0**2))**2
     -        + (lambdainv_boost(1,4)*t0 - 
     -         (2*(lambdainv_boost(2,4)*x0 + 
     -              lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0)
     -            )/(-1 + x0**2 + y0**2 + z0**2))**2
     -      )**2 + (lambdainv_boost(1,2)*
     -     (-1 + Sqrt(1 + 
     -         (lambdainv_boost(1,2)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,2)*x0 + 
     -                 lambdainv_boost(2,3)*y0 + 
     -                 lambdainv_boost(2,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -          + (lambdainv_boost(1,3)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,3)*x0 + 
     -                 lambdainv_boost(3,3)*y0 + 
     -                 lambdainv_boost(3,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -          + (lambdainv_boost(1,4)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,4)*x0 + 
     -                 lambdainv_boost(3,4)*y0 + 
     -                 lambdainv_boost(4,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2)
     -       ))/
     -   ((lambdainv_boost(1,2)*t0 - 
     -        (2*(lambdainv_boost(2,2)*x0 + lambdainv_boost(2,3)*y0 + 
     -             lambdainv_boost(2,4)*z0))/
     -         (-1 + x0**2 + y0**2 + z0**2))**2 + 
     -     (lambdainv_boost(1,3)*t0 - 
     -        (2*(lambdainv_boost(2,3)*x0 + lambdainv_boost(3,3)*y0 + 
     -             lambdainv_boost(3,4)*z0))/
     -         (-1 + x0**2 + y0**2 + z0**2))**2 + 
     -     (lambdainv_boost(1,4)*t0 - 
     -        (2*(lambdainv_boost(2,4)*x0 + lambdainv_boost(3,4)*y0 + 
     -             lambdainv_boost(4,4)*z0))/
     -         (-1 + x0**2 + y0**2 + z0**2))**2)

               y0_invboost_t=
     -         ((-2*(lambdainv_boost(2,3)*x0 + lambdainv_boost(3,3)*y0 + 
     -          lambdainv_boost(3,4)*z0) + 
     -       lambdainv_boost(1,3)*t0*
     -        (-1 + x0**2 + y0**2 + z0**2))*
     -     (2*lambdainv_boost(1,2)*
     -        (lambdainv_boost(1,2)*t0 - 
     -          (2*(lambdainv_boost(2,2)*x0 + 
     -               lambdainv_boost(2,3)*y0 + lambdainv_boost(2,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          ) + 2*lambdainv_boost(1,3)*
     -        (lambdainv_boost(1,3)*t0 - 
     -          (2*(lambdainv_boost(2,3)*x0 + 
     -               lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          ) + 2*lambdainv_boost(1,4)*
     -        (lambdainv_boost(1,4)*t0 - 
     -          (2*(lambdainv_boost(2,4)*x0 + 
     -               lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          ))*Sqrt(1 - 
     -       (lambdainv_boost(1,2)*t0 - 
     -           (2*(lambdainv_boost(2,2)*x0 + 
     -                lambdainv_boost(2,3)*y0 + 
     -                lambdainv_boost(2,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2/
     -        ((lambdainv_boost(1,2)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,2)*x0 + 
     -                  lambdainv_boost(2,3)*y0 + 
     -                  lambdainv_boost(2,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -            + (lambdainv_boost(1,3)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,3)*x0 + 
     -                  lambdainv_boost(3,3)*y0 + 
     -                  lambdainv_boost(3,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -            + (lambdainv_boost(1,4)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,4)*x0 + 
     -                  lambdainv_boost(3,4)*y0 + 
     -                  lambdainv_boost(4,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -          )))/
     -   (2.*(-1 + x0**2 + y0**2 + z0**2)*
     -     Sqrt((lambdainv_boost(1,3)*t0 - 
     -          (2*(lambdainv_boost(2,3)*x0 + 
     -               lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2 + 
     -       (lambdainv_boost(1,4)*t0 - 
     -          (2*(lambdainv_boost(2,4)*x0 + 
     -               lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2)*
     -     Sqrt((lambdainv_boost(1,2)*t0 - 
     -          (2*(lambdainv_boost(2,2)*x0 + 
     -               lambdainv_boost(2,3)*y0 + lambdainv_boost(2,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2 + 
     -       (lambdainv_boost(1,3)*t0 - 
     -          (2*(lambdainv_boost(2,3)*x0 + 
     -               lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2 + 
     -       (lambdainv_boost(1,4)*t0 - 
     -          (2*(lambdainv_boost(2,4)*x0 + 
     -               lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2)*
     -     Sqrt(1 + (lambdainv_boost(1,2)*t0 - 
     -          (2*(lambdainv_boost(2,2)*x0 + 
     -               lambdainv_boost(2,3)*y0 + lambdainv_boost(2,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2 + 
     -       (lambdainv_boost(1,3)*t0 - 
     -          (2*(lambdainv_boost(2,3)*x0 + 
     -               lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2 + 
     -       (lambdainv_boost(1,4)*t0 - 
     -          (2*(lambdainv_boost(2,4)*x0 + 
     -               lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2)) + 
     -  ((-2*(lambdainv_boost(2,3)*x0 + lambdainv_boost(3,3)*y0 + 
     -          lambdainv_boost(3,4)*z0) + 
     -       lambdainv_boost(1,3)*t0*
     -        (-1 + x0**2 + y0**2 + z0**2))*
     -     (((lambdainv_boost(1,2)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,2)*x0 + 
     -                  lambdainv_boost(2,3)*y0 + 
     -                  lambdainv_boost(2,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -           *(2*lambdainv_boost(1,2)*
     -             (lambdainv_boost(1,2)*t0 - 
     -               (2*
     -                  (lambdainv_boost(2,2)*x0 + 
     -                  lambdainv_boost(2,3)*y0 + 
     -                  lambdainv_boost(2,4)*z0))/
     -                (-1 + x0**2 + y0**2 + z0**2))
     -             + 2*lambdainv_boost(1,3)*
     -             (lambdainv_boost(1,3)*t0 - 
     -               (2*
     -                  (lambdainv_boost(2,3)*x0 + 
     -                  lambdainv_boost(3,3)*y0 + 
     -                  lambdainv_boost(3,4)*z0))/
     -                (-1 + x0**2 + y0**2 + z0**2))
     -             + 2*lambdainv_boost(1,4)*
     -             (lambdainv_boost(1,4)*t0 - 
     -               (2*
     -                  (lambdainv_boost(2,4)*x0 + 
     -                  lambdainv_boost(3,4)*y0 + 
     -                  lambdainv_boost(4,4)*z0))/
     -                (-1 + x0**2 + y0**2 + z0**2)))
     -          )/
     -        ((lambdainv_boost(1,2)*t0 - 
     -              (2*
     -                 (lambdainv_boost(2,2)*x0 + 
     -                  lambdainv_boost(2,3)*y0 + 
     -                  lambdainv_boost(2,4)*z0))/
     -               (-1 + x0**2 + y0**2 + z0**2))**
     -            2 + 
     -           (lambdainv_boost(1,3)*t0 - 
     -              (2*
     -                 (lambdainv_boost(2,3)*x0 + 
     -                  lambdainv_boost(3,3)*y0 + 
     -                  lambdainv_boost(3,4)*z0))/
     -               (-1 + x0**2 + y0**2 + z0**2))**
     -            2 + 
     -           (lambdainv_boost(1,4)*t0 - 
     -              (2*
     -                 (lambdainv_boost(2,4)*x0 + 
     -                  lambdainv_boost(3,4)*y0 + 
     -                  lambdainv_boost(4,4)*z0))/
     -               (-1 + x0**2 + y0**2 + z0**2))**
     -            2)**2 - 
     -       (2*lambdainv_boost(1,2)*
     -          (lambdainv_boost(1,2)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,2)*x0 + 
     -                 lambdainv_boost(2,3)*y0 + 
     -                 lambdainv_boost(2,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2)))/
     -        ((lambdainv_boost(1,2)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,2)*x0 + 
     -                  lambdainv_boost(2,3)*y0 + 
     -                  lambdainv_boost(2,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -            + (lambdainv_boost(1,3)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,3)*x0 + 
     -                  lambdainv_boost(3,3)*y0 + 
     -                  lambdainv_boost(3,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -            + (lambdainv_boost(1,4)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,4)*x0 + 
     -                  lambdainv_boost(3,4)*y0 + 
     -                  lambdainv_boost(4,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -          ))*(-1 + 
     -       Sqrt(1 + 
     -         (lambdainv_boost(1,2)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,2)*x0 + 
     -                 lambdainv_boost(2,3)*y0 + 
     -                 lambdainv_boost(2,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -          + (lambdainv_boost(1,3)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,3)*x0 + 
     -                 lambdainv_boost(3,3)*y0 + 
     -                 lambdainv_boost(3,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -          + (lambdainv_boost(1,4)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,4)*x0 + 
     -                 lambdainv_boost(3,4)*y0 + 
     -                 lambdainv_boost(4,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2)
     -       ))/
     -   (2.*(-1 + x0**2 + y0**2 + z0**2)*
     -     Sqrt((lambdainv_boost(1,3)*t0 - 
     -          (2*(lambdainv_boost(2,3)*x0 + 
     -               lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2 + 
     -       (lambdainv_boost(1,4)*t0 - 
     -          (2*(lambdainv_boost(2,4)*x0 + 
     -               lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2)*
     -     Sqrt((lambdainv_boost(1,2)*t0 - 
     -          (2*(lambdainv_boost(2,2)*x0 + 
     -               lambdainv_boost(2,3)*y0 + lambdainv_boost(2,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2 + 
     -       (lambdainv_boost(1,3)*t0 - 
     -          (2*(lambdainv_boost(2,3)*x0 + 
     -               lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2 + 
     -       (lambdainv_boost(1,4)*t0 - 
     -          (2*(lambdainv_boost(2,4)*x0 + 
     -               lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2)*
     -     Sqrt(1 - (lambdainv_boost(1,2)*t0 - 
     -           (2*(lambdainv_boost(2,2)*x0 + 
     -                lambdainv_boost(2,3)*y0 + 
     -                lambdainv_boost(2,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2/
     -        ((lambdainv_boost(1,2)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,2)*x0 + 
     -                  lambdainv_boost(2,3)*y0 + 
     -                  lambdainv_boost(2,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -            + (lambdainv_boost(1,3)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,3)*x0 + 
     -                  lambdainv_boost(3,3)*y0 + 
     -                  lambdainv_boost(3,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -            + (lambdainv_boost(1,4)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,4)*x0 + 
     -                  lambdainv_boost(3,4)*y0 + 
     -                  lambdainv_boost(4,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -          ))) - 
     -  ((-2*(lambdainv_boost(2,3)*x0 + lambdainv_boost(3,3)*y0 + 
     -          lambdainv_boost(3,4)*z0) + 
     -       lambdainv_boost(1,3)*t0*
     -        (-1 + x0**2 + y0**2 + z0**2))*
     -     (2*lambdainv_boost(1,2)*
     -        (lambdainv_boost(1,2)*t0 - 
     -          (2*(lambdainv_boost(2,2)*x0 + 
     -               lambdainv_boost(2,3)*y0 + lambdainv_boost(2,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          ) + 2*lambdainv_boost(1,3)*
     -        (lambdainv_boost(1,3)*t0 - 
     -          (2*(lambdainv_boost(2,3)*x0 + 
     -               lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          ) + 2*lambdainv_boost(1,4)*
     -        (lambdainv_boost(1,4)*t0 - 
     -          (2*(lambdainv_boost(2,4)*x0 + 
     -               lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          ))*Sqrt(1 - 
     -       (lambdainv_boost(1,2)*t0 - 
     -           (2*(lambdainv_boost(2,2)*x0 + 
     -                lambdainv_boost(2,3)*y0 + 
     -                lambdainv_boost(2,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2/
     -        ((lambdainv_boost(1,2)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,2)*x0 + 
     -                  lambdainv_boost(2,3)*y0 + 
     -                  lambdainv_boost(2,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -            + (lambdainv_boost(1,3)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,3)*x0 + 
     -                  lambdainv_boost(3,3)*y0 + 
     -                  lambdainv_boost(3,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -            + (lambdainv_boost(1,4)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,4)*x0 + 
     -                  lambdainv_boost(3,4)*y0 + 
     -                  lambdainv_boost(4,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -          ))*(-1 + 
     -       Sqrt(1 + 
     -         (lambdainv_boost(1,2)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,2)*x0 + 
     -                 lambdainv_boost(2,3)*y0 + 
     -                 lambdainv_boost(2,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -          + (lambdainv_boost(1,3)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,3)*x0 + 
     -                 lambdainv_boost(3,3)*y0 + 
     -                 lambdainv_boost(3,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -          + (lambdainv_boost(1,4)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,4)*x0 + 
     -                 lambdainv_boost(3,4)*y0 + 
     -                 lambdainv_boost(4,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2)
     -       ))/
     -   (2.*(-1 + x0**2 + y0**2 + z0**2)*
     -     Sqrt((lambdainv_boost(1,3)*t0 - 
     -          (2*(lambdainv_boost(2,3)*x0 + 
     -               lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2 + 
     -       (lambdainv_boost(1,4)*t0 - 
     -          (2*(lambdainv_boost(2,4)*x0 + 
     -               lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2)*
     -     ((lambdainv_boost(1,2)*t0 - 
     -           (2*(lambdainv_boost(2,2)*x0 + 
     -                lambdainv_boost(2,3)*y0 + 
     -                lambdainv_boost(2,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2
     -         + (lambdainv_boost(1,3)*t0 - 
     -           (2*(lambdainv_boost(2,3)*x0 + 
     -                lambdainv_boost(3,3)*y0 + 
     -                lambdainv_boost(3,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2
     -         + (lambdainv_boost(1,4)*t0 - 
     -           (2*(lambdainv_boost(2,4)*x0 + 
     -                lambdainv_boost(3,4)*y0 + 
     -                lambdainv_boost(4,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2)
     -       **1.5) - 
     -  ((-2*(lambdainv_boost(2,3)*x0 + lambdainv_boost(3,3)*y0 + 
     -          lambdainv_boost(3,4)*z0) + 
     -       lambdainv_boost(1,3)*t0*
     -        (-1 + x0**2 + y0**2 + z0**2))*
     -     (2*lambdainv_boost(1,3)*
     -        (lambdainv_boost(1,3)*t0 - 
     -          (2*(lambdainv_boost(2,3)*x0 + 
     -               lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          ) + 2*lambdainv_boost(1,4)*
     -        (lambdainv_boost(1,4)*t0 - 
     -          (2*(lambdainv_boost(2,4)*x0 + 
     -               lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          ))*Sqrt(1 - 
     -       (lambdainv_boost(1,2)*t0 - 
     -           (2*(lambdainv_boost(2,2)*x0 + 
     -                lambdainv_boost(2,3)*y0 + 
     -                lambdainv_boost(2,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2/
     -        ((lambdainv_boost(1,2)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,2)*x0 + 
     -                  lambdainv_boost(2,3)*y0 + 
     -                  lambdainv_boost(2,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -            + (lambdainv_boost(1,3)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,3)*x0 + 
     -                  lambdainv_boost(3,3)*y0 + 
     -                  lambdainv_boost(3,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -            + (lambdainv_boost(1,4)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,4)*x0 + 
     -                  lambdainv_boost(3,4)*y0 + 
     -                  lambdainv_boost(4,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -          ))*(-1 + 
     -       Sqrt(1 + 
     -         (lambdainv_boost(1,2)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,2)*x0 + 
     -                 lambdainv_boost(2,3)*y0 + 
     -                 lambdainv_boost(2,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -          + (lambdainv_boost(1,3)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,3)*x0 + 
     -                 lambdainv_boost(3,3)*y0 + 
     -                 lambdainv_boost(3,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -          + (lambdainv_boost(1,4)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,4)*x0 + 
     -                 lambdainv_boost(3,4)*y0 + 
     -                 lambdainv_boost(4,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2)
     -       ))/
     -   (2.*(-1 + x0**2 + y0**2 + z0**2)*
     -     ((lambdainv_boost(1,3)*t0 - 
     -           (2*(lambdainv_boost(2,3)*x0 + 
     -                lambdainv_boost(3,3)*y0 + 
     -                lambdainv_boost(3,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2
     -         + (lambdainv_boost(1,4)*t0 - 
     -           (2*(lambdainv_boost(2,4)*x0 + 
     -                lambdainv_boost(3,4)*y0 + 
     -                lambdainv_boost(4,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2)
     -       **1.5*Sqrt((lambdainv_boost(1,2)*t0 - 
     -          (2*(lambdainv_boost(2,2)*x0 + 
     -               lambdainv_boost(2,3)*y0 + lambdainv_boost(2,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2 + 
     -       (lambdainv_boost(1,3)*t0 - 
     -          (2*(lambdainv_boost(2,3)*x0 + 
     -               lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2 + 
     -       (lambdainv_boost(1,4)*t0 - 
     -          (2*(lambdainv_boost(2,4)*x0 + 
     -               lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2)) + 
     -  (lambdainv_boost(1,3)*
     -     Sqrt(1 - (lambdainv_boost(1,2)*t0 - 
     -           (2*(lambdainv_boost(2,2)*x0 + 
     -                lambdainv_boost(2,3)*y0 + 
     -                lambdainv_boost(2,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2/
     -        ((lambdainv_boost(1,2)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,2)*x0 + 
     -                  lambdainv_boost(2,3)*y0 + 
     -                  lambdainv_boost(2,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -            + (lambdainv_boost(1,3)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,3)*x0 + 
     -                  lambdainv_boost(3,3)*y0 + 
     -                  lambdainv_boost(3,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -            + (lambdainv_boost(1,4)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,4)*x0 + 
     -                  lambdainv_boost(3,4)*y0 + 
     -                  lambdainv_boost(4,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -          ))*(-1 + 
     -       Sqrt(1 + 
     -         (lambdainv_boost(1,2)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,2)*x0 + 
     -                 lambdainv_boost(2,3)*y0 + 
     -                 lambdainv_boost(2,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -          + (lambdainv_boost(1,3)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,3)*x0 + 
     -                 lambdainv_boost(3,3)*y0 + 
     -                 lambdainv_boost(3,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -          + (lambdainv_boost(1,4)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,4)*x0 + 
     -                 lambdainv_boost(3,4)*y0 + 
     -                 lambdainv_boost(4,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2)
     -       ))/
     -   (Sqrt((lambdainv_boost(1,3)*t0 - 
     -          (2*(lambdainv_boost(2,3)*x0 + 
     -               lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2 + 
     -       (lambdainv_boost(1,4)*t0 - 
     -          (2*(lambdainv_boost(2,4)*x0 + 
     -               lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2)*
     -     Sqrt((lambdainv_boost(1,2)*t0 - 
     -          (2*(lambdainv_boost(2,2)*x0 + 
     -               lambdainv_boost(2,3)*y0 + lambdainv_boost(2,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2 + 
     -       (lambdainv_boost(1,3)*t0 - 
     -          (2*(lambdainv_boost(2,3)*x0 + 
     -               lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2 + 
     -       (lambdainv_boost(1,4)*t0 - 
     -          (2*(lambdainv_boost(2,4)*x0 + 
     -               lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2))

               z0_invboost_t=
     -         ((-2*(lambdainv_boost(2,4)*x0 + lambdainv_boost(3,4)*y0 + 
     -          lambdainv_boost(4,4)*z0) + 
     -       lambdainv_boost(1,4)*t0*
     -        (-1 + x0**2 + y0**2 + z0**2))*
     -     (2*lambdainv_boost(1,2)*
     -        (lambdainv_boost(1,2)*t0 - 
     -          (2*(lambdainv_boost(2,2)*x0 + 
     -               lambdainv_boost(2,3)*y0 + lambdainv_boost(2,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          ) + 2*lambdainv_boost(1,3)*
     -        (lambdainv_boost(1,3)*t0 - 
     -          (2*(lambdainv_boost(2,3)*x0 + 
     -               lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          ) + 2*lambdainv_boost(1,4)*
     -        (lambdainv_boost(1,4)*t0 - 
     -          (2*(lambdainv_boost(2,4)*x0 + 
     -               lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          ))*Sqrt(1 - 
     -       (lambdainv_boost(1,2)*t0 - 
     -           (2*(lambdainv_boost(2,2)*x0 + 
     -                lambdainv_boost(2,3)*y0 + 
     -                lambdainv_boost(2,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2/
     -        ((lambdainv_boost(1,2)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,2)*x0 + 
     -                  lambdainv_boost(2,3)*y0 + 
     -                  lambdainv_boost(2,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -            + (lambdainv_boost(1,3)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,3)*x0 + 
     -                  lambdainv_boost(3,3)*y0 + 
     -                  lambdainv_boost(3,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -            + (lambdainv_boost(1,4)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,4)*x0 + 
     -                  lambdainv_boost(3,4)*y0 + 
     -                  lambdainv_boost(4,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -          )))/
     -   (2.*(-1 + x0**2 + y0**2 + z0**2)*
     -     Sqrt((lambdainv_boost(1,3)*t0 - 
     -          (2*(lambdainv_boost(2,3)*x0 + 
     -               lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2 + 
     -       (lambdainv_boost(1,4)*t0 - 
     -          (2*(lambdainv_boost(2,4)*x0 + 
     -               lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2)*
     -     Sqrt((lambdainv_boost(1,2)*t0 - 
     -          (2*(lambdainv_boost(2,2)*x0 + 
     -               lambdainv_boost(2,3)*y0 + lambdainv_boost(2,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2 + 
     -       (lambdainv_boost(1,3)*t0 - 
     -          (2*(lambdainv_boost(2,3)*x0 + 
     -               lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2 + 
     -       (lambdainv_boost(1,4)*t0 - 
     -          (2*(lambdainv_boost(2,4)*x0 + 
     -               lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2)*
     -     Sqrt(1 + (lambdainv_boost(1,2)*t0 - 
     -          (2*(lambdainv_boost(2,2)*x0 + 
     -               lambdainv_boost(2,3)*y0 + lambdainv_boost(2,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2 + 
     -       (lambdainv_boost(1,3)*t0 - 
     -          (2*(lambdainv_boost(2,3)*x0 + 
     -               lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2 + 
     -       (lambdainv_boost(1,4)*t0 - 
     -          (2*(lambdainv_boost(2,4)*x0 + 
     -               lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2)) + 
     -  ((-2*(lambdainv_boost(2,4)*x0 + lambdainv_boost(3,4)*y0 + 
     -          lambdainv_boost(4,4)*z0) + 
     -       lambdainv_boost(1,4)*t0*
     -        (-1 + x0**2 + y0**2 + z0**2))*
     -     (((lambdainv_boost(1,2)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,2)*x0 + 
     -                  lambdainv_boost(2,3)*y0 + 
     -                  lambdainv_boost(2,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -           *(2*lambdainv_boost(1,2)*
     -             (lambdainv_boost(1,2)*t0 - 
     -               (2*
     -                  (lambdainv_boost(2,2)*x0 + 
     -                  lambdainv_boost(2,3)*y0 + 
     -                  lambdainv_boost(2,4)*z0))/
     -                (-1 + x0**2 + y0**2 + z0**2))
     -             + 2*lambdainv_boost(1,3)*
     -             (lambdainv_boost(1,3)*t0 - 
     -               (2*
     -                  (lambdainv_boost(2,3)*x0 + 
     -                  lambdainv_boost(3,3)*y0 + 
     -                  lambdainv_boost(3,4)*z0))/
     -                (-1 + x0**2 + y0**2 + z0**2))
     -             + 2*lambdainv_boost(1,4)*
     -             (lambdainv_boost(1,4)*t0 - 
     -               (2*
     -                  (lambdainv_boost(2,4)*x0 + 
     -                  lambdainv_boost(3,4)*y0 + 
     -                  lambdainv_boost(4,4)*z0))/
     -                (-1 + x0**2 + y0**2 + z0**2)))
     -          )/
     -        ((lambdainv_boost(1,2)*t0 - 
     -              (2*
     -                 (lambdainv_boost(2,2)*x0 + 
     -                  lambdainv_boost(2,3)*y0 + 
     -                  lambdainv_boost(2,4)*z0))/
     -               (-1 + x0**2 + y0**2 + z0**2))**
     -            2 + 
     -           (lambdainv_boost(1,3)*t0 - 
     -              (2*
     -                 (lambdainv_boost(2,3)*x0 + 
     -                  lambdainv_boost(3,3)*y0 + 
     -                  lambdainv_boost(3,4)*z0))/
     -               (-1 + x0**2 + y0**2 + z0**2))**
     -            2 + 
     -           (lambdainv_boost(1,4)*t0 - 
     -              (2*
     -                 (lambdainv_boost(2,4)*x0 + 
     -                  lambdainv_boost(3,4)*y0 + 
     -                  lambdainv_boost(4,4)*z0))/
     -               (-1 + x0**2 + y0**2 + z0**2))**
     -            2)**2 - 
     -       (2*lambdainv_boost(1,2)*
     -          (lambdainv_boost(1,2)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,2)*x0 + 
     -                 lambdainv_boost(2,3)*y0 + 
     -                 lambdainv_boost(2,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2)))/
     -        ((lambdainv_boost(1,2)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,2)*x0 + 
     -                  lambdainv_boost(2,3)*y0 + 
     -                  lambdainv_boost(2,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -            + (lambdainv_boost(1,3)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,3)*x0 + 
     -                  lambdainv_boost(3,3)*y0 + 
     -                  lambdainv_boost(3,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -            + (lambdainv_boost(1,4)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,4)*x0 + 
     -                  lambdainv_boost(3,4)*y0 + 
     -                  lambdainv_boost(4,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -          ))*(-1 + 
     -       Sqrt(1 + 
     -         (lambdainv_boost(1,2)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,2)*x0 + 
     -                 lambdainv_boost(2,3)*y0 + 
     -                 lambdainv_boost(2,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -          + (lambdainv_boost(1,3)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,3)*x0 + 
     -                 lambdainv_boost(3,3)*y0 + 
     -                 lambdainv_boost(3,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -          + (lambdainv_boost(1,4)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,4)*x0 + 
     -                 lambdainv_boost(3,4)*y0 + 
     -                 lambdainv_boost(4,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2)
     -       ))/
     -   (2.*(-1 + x0**2 + y0**2 + z0**2)*
     -     Sqrt((lambdainv_boost(1,3)*t0 - 
     -          (2*(lambdainv_boost(2,3)*x0 + 
     -               lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2 + 
     -       (lambdainv_boost(1,4)*t0 - 
     -          (2*(lambdainv_boost(2,4)*x0 + 
     -               lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2)*
     -     Sqrt((lambdainv_boost(1,2)*t0 - 
     -          (2*(lambdainv_boost(2,2)*x0 + 
     -               lambdainv_boost(2,3)*y0 + lambdainv_boost(2,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2 + 
     -       (lambdainv_boost(1,3)*t0 - 
     -          (2*(lambdainv_boost(2,3)*x0 + 
     -               lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2 + 
     -       (lambdainv_boost(1,4)*t0 - 
     -          (2*(lambdainv_boost(2,4)*x0 + 
     -               lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2)*
     -     Sqrt(1 - (lambdainv_boost(1,2)*t0 - 
     -           (2*(lambdainv_boost(2,2)*x0 + 
     -                lambdainv_boost(2,3)*y0 + 
     -                lambdainv_boost(2,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2/
     -        ((lambdainv_boost(1,2)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,2)*x0 + 
     -                  lambdainv_boost(2,3)*y0 + 
     -                  lambdainv_boost(2,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -            + (lambdainv_boost(1,3)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,3)*x0 + 
     -                  lambdainv_boost(3,3)*y0 + 
     -                  lambdainv_boost(3,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -            + (lambdainv_boost(1,4)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,4)*x0 + 
     -                  lambdainv_boost(3,4)*y0 + 
     -                  lambdainv_boost(4,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -          ))) - 
     -  ((-2*(lambdainv_boost(2,4)*x0 + lambdainv_boost(3,4)*y0 + 
     -          lambdainv_boost(4,4)*z0) + 
     -       lambdainv_boost(1,4)*t0*
     -        (-1 + x0**2 + y0**2 + z0**2))*
     -     (2*lambdainv_boost(1,2)*
     -        (lambdainv_boost(1,2)*t0 - 
     -          (2*(lambdainv_boost(2,2)*x0 + 
     -               lambdainv_boost(2,3)*y0 + lambdainv_boost(2,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          ) + 2*lambdainv_boost(1,3)*
     -        (lambdainv_boost(1,3)*t0 - 
     -          (2*(lambdainv_boost(2,3)*x0 + 
     -               lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          ) + 2*lambdainv_boost(1,4)*
     -        (lambdainv_boost(1,4)*t0 - 
     -          (2*(lambdainv_boost(2,4)*x0 + 
     -               lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          ))*Sqrt(1 - 
     -       (lambdainv_boost(1,2)*t0 - 
     -           (2*(lambdainv_boost(2,2)*x0 + 
     -                lambdainv_boost(2,3)*y0 + 
     -                lambdainv_boost(2,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2/
     -        ((lambdainv_boost(1,2)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,2)*x0 + 
     -                  lambdainv_boost(2,3)*y0 + 
     -                  lambdainv_boost(2,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -            + (lambdainv_boost(1,3)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,3)*x0 + 
     -                  lambdainv_boost(3,3)*y0 + 
     -                  lambdainv_boost(3,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -            + (lambdainv_boost(1,4)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,4)*x0 + 
     -                  lambdainv_boost(3,4)*y0 + 
     -                  lambdainv_boost(4,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -          ))*(-1 + 
     -       Sqrt(1 + 
     -         (lambdainv_boost(1,2)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,2)*x0 + 
     -                 lambdainv_boost(2,3)*y0 + 
     -                 lambdainv_boost(2,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -          + (lambdainv_boost(1,3)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,3)*x0 + 
     -                 lambdainv_boost(3,3)*y0 + 
     -                 lambdainv_boost(3,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -          + (lambdainv_boost(1,4)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,4)*x0 + 
     -                 lambdainv_boost(3,4)*y0 + 
     -                 lambdainv_boost(4,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2)
     -       ))/
     -   (2.*(-1 + x0**2 + y0**2 + z0**2)*
     -     Sqrt((lambdainv_boost(1,3)*t0 - 
     -          (2*(lambdainv_boost(2,3)*x0 + 
     -               lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2 + 
     -       (lambdainv_boost(1,4)*t0 - 
     -          (2*(lambdainv_boost(2,4)*x0 + 
     -               lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2)*
     -     ((lambdainv_boost(1,2)*t0 - 
     -           (2*(lambdainv_boost(2,2)*x0 + 
     -                lambdainv_boost(2,3)*y0 + 
     -                lambdainv_boost(2,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2
     -         + (lambdainv_boost(1,3)*t0 - 
     -           (2*(lambdainv_boost(2,3)*x0 + 
     -                lambdainv_boost(3,3)*y0 + 
     -                lambdainv_boost(3,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2
     -         + (lambdainv_boost(1,4)*t0 - 
     -           (2*(lambdainv_boost(2,4)*x0 + 
     -                lambdainv_boost(3,4)*y0 + 
     -                lambdainv_boost(4,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2)
     -       **1.5) - 
     -  ((-2*(lambdainv_boost(2,4)*x0 + lambdainv_boost(3,4)*y0 + 
     -          lambdainv_boost(4,4)*z0) + 
     -       lambdainv_boost(1,4)*t0*
     -        (-1 + x0**2 + y0**2 + z0**2))*
     -     (2*lambdainv_boost(1,3)*
     -        (lambdainv_boost(1,3)*t0 - 
     -          (2*(lambdainv_boost(2,3)*x0 + 
     -               lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          ) + 2*lambdainv_boost(1,4)*
     -        (lambdainv_boost(1,4)*t0 - 
     -          (2*(lambdainv_boost(2,4)*x0 + 
     -               lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          ))*Sqrt(1 - 
     -       (lambdainv_boost(1,2)*t0 - 
     -           (2*(lambdainv_boost(2,2)*x0 + 
     -                lambdainv_boost(2,3)*y0 + 
     -                lambdainv_boost(2,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2/
     -        ((lambdainv_boost(1,2)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,2)*x0 + 
     -                  lambdainv_boost(2,3)*y0 + 
     -                  lambdainv_boost(2,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -            + (lambdainv_boost(1,3)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,3)*x0 + 
     -                  lambdainv_boost(3,3)*y0 + 
     -                  lambdainv_boost(3,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -            + (lambdainv_boost(1,4)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,4)*x0 + 
     -                  lambdainv_boost(3,4)*y0 + 
     -                  lambdainv_boost(4,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -          ))*(-1 + 
     -       Sqrt(1 + 
     -         (lambdainv_boost(1,2)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,2)*x0 + 
     -                 lambdainv_boost(2,3)*y0 + 
     -                 lambdainv_boost(2,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -          + (lambdainv_boost(1,3)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,3)*x0 + 
     -                 lambdainv_boost(3,3)*y0 + 
     -                 lambdainv_boost(3,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -          + (lambdainv_boost(1,4)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,4)*x0 + 
     -                 lambdainv_boost(3,4)*y0 + 
     -                 lambdainv_boost(4,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2)
     -       ))/
     -   (2.*(-1 + x0**2 + y0**2 + z0**2)*
     -     ((lambdainv_boost(1,3)*t0 - 
     -           (2*(lambdainv_boost(2,3)*x0 + 
     -                lambdainv_boost(3,3)*y0 + 
     -                lambdainv_boost(3,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2
     -         + (lambdainv_boost(1,4)*t0 - 
     -           (2*(lambdainv_boost(2,4)*x0 + 
     -                lambdainv_boost(3,4)*y0 + 
     -                lambdainv_boost(4,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2)
     -       **1.5*Sqrt((lambdainv_boost(1,2)*t0 - 
     -          (2*(lambdainv_boost(2,2)*x0 + 
     -               lambdainv_boost(2,3)*y0 + lambdainv_boost(2,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2 + 
     -       (lambdainv_boost(1,3)*t0 - 
     -          (2*(lambdainv_boost(2,3)*x0 + 
     -               lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2 + 
     -       (lambdainv_boost(1,4)*t0 - 
     -          (2*(lambdainv_boost(2,4)*x0 + 
     -               lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2)) + 
     -  (lambdainv_boost(1,4)*
     -     Sqrt(1 - (lambdainv_boost(1,2)*t0 - 
     -           (2*(lambdainv_boost(2,2)*x0 + 
     -                lambdainv_boost(2,3)*y0 + 
     -                lambdainv_boost(2,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2/
     -        ((lambdainv_boost(1,2)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,2)*x0 + 
     -                  lambdainv_boost(2,3)*y0 + 
     -                  lambdainv_boost(2,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -            + (lambdainv_boost(1,3)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,3)*x0 + 
     -                  lambdainv_boost(3,3)*y0 + 
     -                  lambdainv_boost(3,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -            + (lambdainv_boost(1,4)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,4)*x0 + 
     -                  lambdainv_boost(3,4)*y0 + 
     -                  lambdainv_boost(4,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -          ))*(-1 + 
     -       Sqrt(1 + 
     -         (lambdainv_boost(1,2)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,2)*x0 + 
     -                 lambdainv_boost(2,3)*y0 + 
     -                 lambdainv_boost(2,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -          + (lambdainv_boost(1,3)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,3)*x0 + 
     -                 lambdainv_boost(3,3)*y0 + 
     -                 lambdainv_boost(3,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -          + (lambdainv_boost(1,4)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,4)*x0 + 
     -                 lambdainv_boost(3,4)*y0 + 
     -                 lambdainv_boost(4,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2)
     -       ))/
     -   (Sqrt((lambdainv_boost(1,3)*t0 - 
     -          (2*(lambdainv_boost(2,3)*x0 + 
     -               lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2 + 
     -       (lambdainv_boost(1,4)*t0 - 
     -          (2*(lambdainv_boost(2,4)*x0 + 
     -               lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2)*
     -     Sqrt((lambdainv_boost(1,2)*t0 - 
     -          (2*(lambdainv_boost(2,2)*x0 + 
     -               lambdainv_boost(2,3)*y0 + lambdainv_boost(2,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2 + 
     -       (lambdainv_boost(1,3)*t0 - 
     -          (2*(lambdainv_boost(2,3)*x0 + 
     -               lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2 + 
     -       (lambdainv_boost(1,4)*t0 - 
     -          (2*(lambdainv_boost(2,4)*x0 + 
     -               lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2))


              if ((abs(y0).lt.10.0d0**(-10))
     -         .and.(abs(z0).lt.10.0d0**(-10)) ) then
                y0_invboost_t=0
                z0_invboost_t=0
              end if

              if ((abs(x0).lt.10.0d0**(-10))
     -         .and.(abs(z0).lt.10.0d0**(-10)) ) then
                x0_invboost_t=0
                z0_invboost_t=0
              end if

              if ((abs(x0).lt.10.0d0**(-10))
     -         .and.(abs(y0).lt.10.0d0**(-10)) ) then
                x0_invboost_t=0
                y0_invboost_t=0
              end if

              if ((abs(x0).lt.10.0d0**(-10))
     -         .and.(abs(y0).lt.10.0d0**(-10))
     -          .and.(abs(z0).lt.10.0d0**(-10)) ) then
                x0_invboost_t=0
                y0_invboost_t=0
                z0_invboost_t=0
              end if


              r_invboost_t=(1/(2*r_invboost))*
     &             ((1-boost_ex**2)*2*xb_invboost*x0_invboost_t+
     &              (1-boost_ey**2)*2*yb_invboost*y0_invboost_t+
     &              (1-boost_ez**2)*2*zb_invboost*z0_invboost_t)

              derexp0=(-1/boost_delta**2)*2*
     &         (r_invboost-boost_r0)*r_invboost_t

               f_t_n(i,j,k)=f_n(i,j,k)*derexp0

               if (rho0.lt.10.0d0**(-10)) then
                f_t_n(i,j,k)=0
               end if

               phi1_t_n(i,j,k)=phi1_t_n(i,j,k)+f_t_n(i,j,k)
               gb_tt_t_n(i,j,k)=gb_tt_t_n(i,j,k)+f_t_n(i,j,k)
               gb_tx_t_n(i,j,k)=gb_tx_t_n(i,j,k)+f_t_n(i,j,k)
               gb_ty_t_n(i,j,k)=gb_ty_t_n(i,j,k)+f_t_n(i,j,k)
               gb_tz_t_n(i,j,k)=gb_tz_t_n(i,j,k)+f_t_n(i,j,k)
               gb_xx_t_n(i,j,k)=gb_xx_t_n(i,j,k)+f_t_n(i,j,k)
               gb_xy_t_n(i,j,k)=gb_xy_t_n(i,j,k)+f_t_n(i,j,k)
               gb_xz_t_n(i,j,k)=gb_xz_t_n(i,j,k)+f_t_n(i,j,k)
               gb_yy_t_n(i,j,k)=gb_yy_t_n(i,j,k)+f_t_n(i,j,k)
               gb_yz_t_n(i,j,k)=gb_yz_t_n(i,j,k)+f_t_n(i,j,k)
               gb_zz_t_n(i,j,k)=gb_zz_t_n(i,j,k)+f_t_n(i,j,k)


               phi1_nm1(i,j,k)=phi1_n(i,j,k) - phi1_t_n(i,j,k)*dt
               phi1_np1(i,j,k)=phi1_n(i,j,k) + phi1_t_n(i,j,k)*dt
               gb_tt_nm1(i,j,k)=gb_tt_n(i,j,k) - gb_tt_t_n(i,j,k)*dt
               gb_tt_np1(i,j,k)=gb_tt_n(i,j,k) + gb_tt_t_n(i,j,k)*dt
               gb_tx_nm1(i,j,k)=gb_tx_n(i,j,k) - gb_tx_t_n(i,j,k)*dt
               gb_tx_np1(i,j,k)=gb_tx_n(i,j,k) + gb_tx_t_n(i,j,k)*dt
               gb_ty_nm1(i,j,k)=gb_ty_n(i,j,k) - gb_ty_t_n(i,j,k)*dt
               gb_ty_np1(i,j,k)=gb_ty_n(i,j,k) + gb_ty_t_n(i,j,k)*dt
               gb_tz_nm1(i,j,k)=gb_tz_n(i,j,k) - gb_tz_t_n(i,j,k)*dt
               gb_tz_np1(i,j,k)=gb_tz_n(i,j,k) + gb_tz_t_n(i,j,k)*dt
               gb_xx_nm1(i,j,k)=gb_xx_n(i,j,k) - gb_xx_t_n(i,j,k)*dt
               gb_xx_np1(i,j,k)=gb_xx_n(i,j,k) + gb_xx_t_n(i,j,k)*dt
               gb_xy_nm1(i,j,k)=gb_xy_n(i,j,k) - gb_xy_t_n(i,j,k)*dt
               gb_xy_np1(i,j,k)=gb_xy_n(i,j,k) + gb_xy_t_n(i,j,k)*dt
               gb_xz_nm1(i,j,k)=gb_xz_n(i,j,k) - gb_xz_t_n(i,j,k)*dt
               gb_xz_np1(i,j,k)=gb_xz_n(i,j,k) + gb_xz_t_n(i,j,k)*dt
               gb_yy_nm1(i,j,k)=gb_yy_n(i,j,k) - gb_yy_t_n(i,j,k)*dt
               gb_yy_np1(i,j,k)=gb_yy_n(i,j,k) + gb_yy_t_n(i,j,k)*dt
               gb_yz_nm1(i,j,k)=gb_yz_n(i,j,k) - gb_yz_t_n(i,j,k)*dt
               gb_yz_np1(i,j,k)=gb_yz_n(i,j,k) + gb_yz_t_n(i,j,k)*dt
               gb_zz_nm1(i,j,k)=gb_zz_n(i,j,k) - gb_zz_t_n(i,j,k)*dt
               gb_zz_np1(i,j,k)=gb_zz_n(i,j,k) + gb_zz_t_n(i,j,k)*dt



!!DEBUG
!              if ((is_nan(f_n(i,j,k)))
!     &         .or.(is_nan(f_t_n(i,j,k)))) then
!                write (*,*) "x0,y0,z0=",x0,y0,z0
!                write (*,*) "boost_xu0,boost_yu0,boost_zu0=",boost_xu0,boost_yu0,boost_zu0
!                write (*,*) "A0,boost_r0,boost_delta=",A0,boost_r0,boost_delta
!               do a=1,4
!                do b=1,4
!                  write (*,*) "a,b,lambda_boost(a,b)="
!     &             ,a,b,lambda_boost(a,b)
!                  write (*,*) "a,b,lambdainv_boost(a,b)="
!     &             ,a,b,lambdainv_boost(a,b)
!                end do
!               end do
!                  write (*,*) "x0_invboost=",x0_invboost
!                  write (*,*) "y0_invboost=",y0_invboost
!                  write (*,*) "z0_invboost=",z0_invboost
!                  write (*,*) "r_invboost=",r_invboost
!                  write (*,*) "x0_invboost_t=",x0_invboost_t
!                  write (*,*) "y0_invboost_t=",y0_invboost_t
!                  write (*,*) "z0_invboost_t=",z0_invboost_t
!                  write (*,*) "r_invboost_t=",r_invboost_t
!                  write (*,*) "f_n(i,j,k)=",f_n(i,j,k)
!                  write (*,*) "f_t_n(i,j,k)=",f_t_n(i,j,k)
!              end if
!
!
!              if ((abs(x0-0.375d0).lt.10.0d0**(-10)).and.
!     &            (abs(y0-0.125d0).lt.10.0d0**(-10)).and.
!     &            (abs(z0-0.0d0).lt.10.0d0**(-10))) then
!                write (*,*) "x0,y0,z0=",x0,y0,z0
!                write (*,*) "boost_xu0,boost_yu0,boost_zu0=",boost_xu0,boost_yu0,boost_zu0
!                write (*,*) "A0,boost_r0,boost_delta=",A0,boost_r0,boost_delta
!               do a=1,4
!                do b=1,4
!                  write (*,*) "a,b,lambda_boost(a,b)="
!     &             ,a,b,lambda_boost(a,b)
!                  write (*,*) "a,b,lambdainv_boost(a,b)="
!     &             ,a,b,lambdainv_boost(a,b)
!                end do
!               end do
!                  write (*,*) "xb_invboost=",xb_invboost
!                  write (*,*) "yb_invboost=",yb_invboost
!                  write (*,*) "zb_invboost=",zb_invboost
!                  write (*,*) "r_invboost=",r_invboost
!                  write (*,*) "f(i,j,k)=",f(i,j,k)
!                  write (*,*) "f_t(i,j,k)=",f_t(i,j,k)
!              end if

             end if
            end do
           end do
        end do
       end if

        return
        end
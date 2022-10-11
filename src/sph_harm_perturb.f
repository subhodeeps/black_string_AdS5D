c----------------------------------------------------------------------
c add l=2,m=2 scalar and vector perturbations to metric components in Cartesian (code) coords
c----------------------------------------------------------------------

        subroutine sph_harm_perturb(
     &                     phi1,
     &                     gb_tt,
     &                     gb_tx,
     &                     gb_ty,
     &                     gb_tz,
     &                     gb_xx,
     &                     gb_xy,
     &                     gb_xz,
     &                     gb_yy,
     &                     gb_yz,
     &                     gb_zz,
     &                     amp_Y,amp_V,
     &                     L,x,y,z,dt,chr,exc,Nx,Ny,Nz,
     &                     rhoamp_a,rhoamp_b)


        implicit none

        integer Nx,Ny,Nz
        real*8 phi1(Nx,Ny,Nz)
        real*8 x(Nx),y(Ny),z(Nz),dt
        real*8 amp_Y,amp_V
        real*8 L
        real*8 chr(Nx,Ny,Nz),exc
        real*8 rhoamp_a,rhoamp_b
        real*8 f0,trans

        real*8 gb_tt(Nx,Ny,Nz)
        real*8 gb_tx(Nx,Ny,Nz)
        real*8 gb_ty(Nx,Ny,Nz)
        real*8 gb_tz(Nx,Ny,Nz)
        real*8 gb_xx(Nx,Ny,Nz)
        real*8 gb_xy(Nx,Ny,Nz)
        real*8 gb_xz(Nx,Ny,Nz)
        real*8 gb_yy(Nx,Ny,Nz)
        real*8 gb_yz(Nx,Ny,Nz)
        real*8 gb_zz(Nx,Ny,Nz)

        real*8 pert_gbcart_ll(4,4),pert_gbqssph_ll(4,4)
        real*8 dxqssph_dxcar(4,4)

        real*8 Re_Y_22
        real*8 Re_V_122_chi,Re_V_122_xi
        real*8 Re_V_222_chi,Re_V_222_xi

        logical is_nan

        integer i,j,k
        integer stype
        integer a,b,c,d
        real*8 r,x0,y0,z0,rho0,xi0,chi0,csr,xb,yb,zb


        real*8 PI
        parameter (PI=3.141592653589793d0)

        ! initialize fixed-size variables
        data i,j,k/0,0,0/
        data r,x0,y0,z0,rho0,xi0,chi0,csr,xb,yb,zb
     &       /0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/

     	data pert_gbcart_ll/16*0.0/
     	data pert_gbqssph_ll/16*0.0/
     	data dxqssph_dxcar/16*0.0/

        !--------------------------------------------------------------

       if ((abs(amp_Y).gt.10.0d0**(-10))
     &  .or.(abs(amp_V).gt.10.0d0**(-10))) then

        do i=1,Nx
         do j=1,Ny
          do k=1,Nz
           if (chr(i,j,k).ne.exc) then
              x0=x(i)
              y0=y(j)
              z0=z(k)
              rho0=sqrt(x0**2+y0**2+z0**2)
              if (rho0.ne.0.0d0) then
               chi0=(1/PI)*acos(x0/rho0)
              end if
              if ((y0.ne.0.0d0).or.(z0.ne.0.0d0)) then
                xi0=(1/(2*PI))*atan2(z0,y0)
                if (xi0.lt.0) xi0=xi0+1
              end if

              f0=trans(rho0,rhoamp_a,rhoamp_b)

             !real part of scalar spherical harmonic with l=2,m=2 
              Re_Y_22=(1/4.0d0)*sqrt(15/(2*PI))
     &		         *cos(2*2*PI*xi0)*(sin(PI*chi0)**2)

             !real part of even (s=1) vector spherical harmonic with l=2,m=2
             !in (chi,xi) coords, where chi=theta/PI,xi=phi/(2*PI) and (theta,phi) are the usual 
             !coordinates on the 2-sphere
              Re_V_122_chi=PI*(1/8.0d0)*sqrt(5/PI)
     &		         *cos(2*2*PI*xi0)*sin(2*PI*chi0)
              Re_V_122_xi=2*PI*(-1/4.0d0)*sqrt(5/PI)
     &		         *(sin(PI*chi0)**2)*sin(2*2*PI*xi0)

             !real part of odd (s=2) vector spherical harmonic with l=2,m=2
             !in (chi,xi) coords, where chi=theta/PI,xi=phi/(2*PI) and (theta,phi) are the usual 
             !coordinates on the 2-sphere
              Re_V_222_chi=PI*(-1/4.0d0)*sqrt(5/PI)
     &		         *sin(PI*chi0)*sin(2*2*PI*xi0)
              Re_V_222_xi=2*PI*(-1/4.0d0)*sqrt(5/PI)*
     &        		 cos(PI*chi0)*cos(2*2*PI*xi0)*(sin(PI*chi0)**2)

     		pert_gbqssph_ll(1,1)=0
     		pert_gbqssph_ll(1,2)=0
     		pert_gbqssph_ll(1,3)=f0*16*PI/3*(1-rho0)*amp_V
     &		         *(Re_V_122_chi+Re_V_222_chi)
     		pert_gbqssph_ll(1,4)=f0*16*PI/3*(1-rho0)*amp_V
     &		         *(Re_V_122_xi+Re_V_222_xi)
     		pert_gbqssph_ll(2,2)=f0*(1.0d0/3.0d0)*(64*PI**3)
     &             /(8*PI**2)*(1-rho0)*amp_Y*Re_Y_22
     		pert_gbqssph_ll(2,3)=0
     		pert_gbqssph_ll(2,4)=0
     		pert_gbqssph_ll(3,3)=f0*(1.0d0/3.0d0)*(64*PI**3)
     &             /12*(1-rho0)*amp_Y*Re_Y_22
     		pert_gbqssph_ll(3,4)=0
     		pert_gbqssph_ll(4,4)=f0*(1.0d0/3.0d0)*64*PI**3*(sin(PI*chi0)**2)
     &		         /3*(1-rho0)*amp_Y*Re_Y_22

     	      do a=1,3
           	 do b=a+1,4
              pert_gbqssph_ll(b,a)=pert_gbqssph_ll(a,b)
             end do
            end do

!DEBUG!!!!!
!        if ((abs(x0-(0.875d0)).lt.10.0d0**(-10))
!     -  .and.(abs(y0-(0.25d0)).lt.10.0d0**(-10))
!     -  .and.(abs(z0-(0.125d0)).lt.10.0d0**(-10))) then
!
!         write (*,*) "x0,y0,z0=",x0,y0,z0
!        write (*,*) "rho0,chi0,xi0=",rho0,chi0,xi0
!        write (*,*) "amp_Y,amp_V=",amp_Y,amp_V
!        write (*,*) "Re_Y_22=",Re_Y_22
!        write (*,*) "Re_V_122_chi,Re_V_122_xi=",
!     -   Re_V_122_chi,Re_V_122_xi
!        write (*,*) "Re_V_222_chi,Re_V_222_xi=",
!     -   Re_V_222_chi,Re_V_222_xi
!        do a=1,4
!          do b=1,4
!            write (*,*) "a,b,pert_gbqssph_ll(a,b)=",
!     -       a,b,pert_gbqssph_ll(a,b)
!          end do
!        end do
!
!        end if
!!!!!!!!!!!!!!!!!


			dxqssph_dxcar(1,1)=1
			dxqssph_dxcar(1,2)=0
			dxqssph_dxcar(1,3)=0
			dxqssph_dxcar(1,4)=0

			dxqssph_dxcar(2,1)=0
			dxqssph_dxcar(2,2)=x0/rho0
			dxqssph_dxcar(2,3)=y0/rho0
			dxqssph_dxcar(2,4)=z0/rho0

	        dxqssph_dxcar(3,1)=0
	        dxqssph_dxcar(3,2)=-(Sqrt(rho0**2 - x0**2)
     &		    /(PI*rho0**2))
	        dxqssph_dxcar(3,3)=(x0*y0)
     &		    /(PI*rho0**2*Sqrt(rho0**2 - x0**2))
	        dxqssph_dxcar(3,4)=(x0*z0)
     &		    /(PI*rho0**2*Sqrt(rho0**2 - x0**2))

			dxqssph_dxcar(4,1)=0
			dxqssph_dxcar(4,2)=0
			dxqssph_dxcar(4,3)=-(z0
     &		     /(2*PI*(y0**2 + z0**2)))
			dxqssph_dxcar(4,4)=y0
     &		     /(2*PI*(y0**2 + z0**2))

	        !compute Cartesian quantities in terms of quasi-spherical ones
	        do a=1,4
	         do b=1,4
	          pert_gbcart_ll(a,b)=0
	          do c=1,4
	           do d=1,4
	            pert_gbcart_ll(a,b)=
     &		     pert_gbcart_ll(a,b)
     &	   		  +dxqssph_dxcar(c,a)
     &			  *dxqssph_dxcar(d,b) 
     &            	*pert_gbqssph_ll(c,d)
	     	   end do
	     	  end do
	     	 end do
	     	end do

	       !some of the dxqssph_dxcar diverge at y=z=0 so we need to consider this case separately
    	    if ((abs(y0).lt.10.0d0**(-10)).and.
     &	     (abs(z0).lt.10.0d0**(-10))) then
     			do a=1,4
	         	 do b=1,4
	         	 	pert_gbcart_ll(a,b)=0
	         	 end do
	         	end do
     	    end if


!DEBUG!!!!!!!!!!!
!         if ((is_nan(pert_gbcart_ll(1,1))).or.
!     -      (is_nan(pert_gbcart_ll(1,2))).or.
!     -      (is_nan(pert_gbcart_ll(1,3))).or.
!     -      (is_nan(pert_gbcart_ll(1,4))).or.
!     -      (is_nan(pert_gbcart_ll(2,2))).or.
!     -      (is_nan(pert_gbcart_ll(2,3))).or.
!     -      (is_nan(pert_gbcart_ll(2,4))).or.
!     -      (is_nan(pert_gbcart_ll(3,3))).or.
!     -      (is_nan(pert_gbcart_ll(3,4))).or.
!     -      (is_nan(pert_gbcart_ll(4,4))) ) then
!
!
!        write (*,*) "x0,y0,z0=",x0,y0,z0
!        write (*,*) "rho0,chi0,xi0=",rho0,chi0,xi0
!        write (*,*) "amp_Y,amp_V=",amp_Y,amp_V
!        do a=1,4
!          do b=1,4
!            write (*,*) "a,b,pert_gbcart_ll(a,b)=",
!     -       a,b,pert_gbcart_ll(a,b)
!          end do
!        end do
!
!
!
!        end if
!!!!!!!!!!!!!!!!!!!!!
!DEBUG!!!!!
!        if ((abs(x0-(0.875d0)).lt.10.0d0**(-10))
!     -  .and.(abs(y0-(0.25d0)).lt.10.0d0**(-10))
!     -  .and.(abs(z0-(0.125d0)).lt.10.0d0**(-10))) then
!
!         write (*,*) "x0,y0,z0=",x0,y0,z0
!        write (*,*) "rho0,chi0,xi0=",rho0,chi0,xi0
!        write (*,*) "amp_Y,amp_V=",amp_Y,amp_V
!        write (*,*) "Re_Y_22=",Re_Y_22
!        write (*,*) "Re_V_122_chi,Re_V_122_xi=",
!     -   Re_V_122_chi,Re_V_122_xi
!        write (*,*) "Re_V_222_chi,Re_V_222_xi=",
!     -   Re_V_222_chi,Re_V_222_xi
!        do a=1,4
!          do b=1,4
!            write (*,*) "a,b,pert_gbcart_ll(a,b)=",
!     -       a,b,pert_gbcart_ll(a,b)
!          end do
!        end do
!
!        end if
!!!!!!!!!!!!!!!!!

     	         !phi1(i,j,k)=phi1(i,j,k)+0.0d0
               gb_tt(i,j,k)=gb_tt(i,j,k)+pert_gbcart_ll(1,1)
               gb_tx(i,j,k)=gb_tx(i,j,k)+pert_gbcart_ll(1,2)
               gb_ty(i,j,k)=gb_ty(i,j,k)+pert_gbcart_ll(1,3)
               gb_tz(i,j,k)=gb_tz(i,j,k)+pert_gbcart_ll(1,4)
               gb_xx(i,j,k)=gb_xx(i,j,k)+pert_gbcart_ll(2,2)
               gb_xy(i,j,k)=gb_xy(i,j,k)+pert_gbcart_ll(2,3)
               gb_xz(i,j,k)=gb_xz(i,j,k)+pert_gbcart_ll(2,4)
               gb_yy(i,j,k)=gb_yy(i,j,k)+pert_gbcart_ll(3,3)
               gb_yz(i,j,k)=gb_yz(i,j,k)+pert_gbcart_ll(3,4)
               gb_zz(i,j,k)=gb_zz(i,j,k)+pert_gbcart_ll(4,4)




           end if
          end do
         end do
        end do
       end if



        return
        end

c----------------------------------------------------------------------
c Relaxation, Lop, Residual routines for MG solver for zeta
c satisfying L.zeta=0, with induced MG right-hand side R giving L.zeta=R
c----------------------------------------------------------------------

c-----------------------------------------------------------------------
c if 
c 
c action = 1 : performs 1 RB relaxation sweep of zeta
c              (_lop, _res vars unused) ... sets norm
c action = 2 : calculates MG residual L.zeta-R --> _res vars ... sets norm
c              (_lop unused)
c action = 3 : calculates normal residual L.zeta --> _lop
c              (_res, _rhs unused)

c-----------------------------------------------------------------------      
        subroutine mg_sup(action,zeta,zeta_rhs,zeta_lop,zeta_res,phi1,
     &                    L,cmask,phys_bdy,chr,ex,x,y,z,norm,Nx,Ny,Nz)
        implicit none
        integer Nx,Ny,Nz,action
        integer phys_bdy(6)
        real*8 zeta(Nx,Ny,Nz),zeta_rhs(Nx,Ny,Nz),zeta_lop(Nx,Ny,Nz)
        real*8 zeta_res(Nx,Ny,Nz)
        real*8 phi1(Nx,Ny,Nz)
        real*8 cmask(Nx,Ny,Nz),chr(Nx,Ny,Nz)
        real*8 x(Nx),y(Ny),z(Nz),norm,ex,L
        real*8 trhoE_grad,trhoE_ptl
        real*8 zeta_x(4),ddzeta,ddzeta_Jac,grad_zeta_sq
        real*8 phi1_x(4),ddphi1,ddphi1_Jac,grad_phi1_sq

        real*8 phi10(Nx,Ny,Nz)

        integer relax,lop,residual,cdiff_method
        parameter (relax=1,lop=3,residual=2)
        real*8 csr,snr,x0,y0,z0,rho0,x2,y2,x3,y3,x4,y4
        real*8 h,dx,dy,dz

        real*8 zeta0
        real*8 phi10_0

        real*8 Jac,res,rhs,new_rhs
        real*8 lambda4
        real*8 PI
        parameter (PI=0.3141592653589793D1)

        real*8 ddzeta_J_xx,ddzeta_J_yy

        integer i,j,k,pass,sum,ii,jj
        integer is

        logical first,do_red_black,extrap
        parameter (extrap=.true.) !matches extrap=.true. in df2_int
        data first/.true./
        data do_red_black/.true./
        save first

        ! initialize fixed-size variables
        data i,j,k,pass,sum,ii,jj/0,0,0,0,0,0,0/
        data is/0/

        data csr,snr,x0,y0,rho0/0.0,0.0,0.0,0.0,0.0/
        data x2,y2,x3,y3,x4,y4/0.0,0.0,0.0,0.0,0.0,0.0/
        data h,dx,dy,dz/0.0,0.0,0.0,0.0/

        data zeta0/0.0/
        data phi10_0/0.0/

        data Jac,res,rhs,new_rhs/0.0,0.0,0.0,0.0/
        data lambda4/0.0/

        data ddzeta_J_xx,ddzeta_J_yy/0.0,0.0/

        !--------------------------------------------------------------

        first=.false.

        norm=0
        sum=0

        dx=(x(2)-x(1))
        dy=(y(2)-y(1))
        dz=(z(2)-z(1))

        h=x(2)-x(1)

        if ((abs(1-(y(2)-y(1))/h).gt.1.0d-5)
     &     .or.(abs(1-(z(2)-z(1))/h).gt.1.0d-5)) then
           write(*,*) 'error ... relax() expects dx=dy and dx=dz'
           stop
        end if

        rhs=0

        ! sets AdS4D cosmological constant
        lambda4=-3/L/L

        ! manually reconstruct phi10=phi1*(1-rho^2)^3 
        do i=1,Nx
          do j=1,Ny
           do k=1,Nz
            x0=x(i)
            y0=y(j)
            z0=z(k)
            rho0=sqrt(x0**2+y0**2+z0**2)
            if (phi1(i,j,k).ne.0) phi10(i,j,k)=phi1(i,j,k)
     &                                         *(1-rho0**2)**2
            if (phi1(i,j,k).eq.0) phi10(i,j,k)=0
           end do
          end do
        end do

        ! (REGION) interior, solve L.zeta=0 Hamiltonian constraint 
        do pass=0,1
          do i=2,Nx-1
            do j=2,Ny-1
             do k=2,Nz-1
              x0=x(i)
              y0=y(j)
              z0=z(k)
              x2=x0*x0
              y2=y0*y0
              x3=x0*x0*x0
              y3=y0*y0*y0
              x4=x0*x0*x0*x0
              y4=y0*y0*y0*y0
              rho0=sqrt(x0**2+y0**2+z0**2)

              if (
     &            ((do_red_black.and.mod(i+j+k+pass,2).eq.0)
     &             .or.(.not.do_red_black.and.pass.eq.0)) 
     &             .and.(cmask(i,j,k).eq.1)
     &             .and.(chr(i,j,k).ne.ex)      
     &           ) then

                ! fill in zeta 
                zeta0=zeta(i,j,k)

                ! fill in phi10_0
                phi10_0=phi10(i,j,k)

                ! computes initial energy density at i,j, with initial data
                ! time-symmetric so phi_t=0, and free scalar so V(phi)=0
                call df_int(phi10,phi1_x,ddphi1,ddphi1_Jac,
     &                      grad_phi1_sq,
     &                      x,y,z,i,j,k,chr,L,ex,Nx,Ny,Nz)
                trhoE_grad=grad_phi1_sq/2
                trhoE_ptl=0

                ! computes normal residual L.zeta
                !(NOTE: the physical energy density rhoE is such that
                ! rhoE_grad=trhoE_grad*zeta^(-2), rhoE_ptl=trhoE_ptl)
                call df_int(zeta,zeta_x,ddzeta,ddzeta_Jac,
     &                      grad_zeta_sq,
     &                      x,y,z,i,j,k,chr,L,ex,Nx,Ny,Nz)
                res=ddzeta
     &              -lambda4*zeta0/4
     &              +(lambda4+8*PI*trhoE_ptl)*(zeta0**5)/4
     &              +8*PI*trhoE_grad*zeta0/4
            
                ! computes MG residual L.zeta-R
                rhs=res-zeta_rhs(i,j,k)

                Jac=ddzeta_Jac
     &              -( lambda4/4 )
     &              +(lambda4+8*PI*trhoE_ptl)*(zeta0**4)*5/4
     &              +8*PI*trhoE_grad/4

                ! performs action
                if (action.eq.residual) then
                  zeta_res(i,j,k)=rhs
                else if (action.eq.lop) then
                  zeta_lop(i,j,k)=res
                else if (action.eq.relax) then
                  zeta(i,j,k)=zeta(i,j,k)-rhs/Jac
                end if

                norm=norm+rhs**2
                sum=sum+1
                
              end if
             end do
            end do
          end do
        end do

!        ! (REGION) y=0 axis, solve d/dy(zeta)=0  
!        if (phys_bdy(3).ne.0) then
!          do i=2,Nx-1
!            do k=2,Nz-1
!
!            ! computes normal residual d/dy(zeta)
!            res=zeta(i,1,k)-(4*zeta(i,2,k)-zeta(i,3,k))/3
!
!            ! computes MG residual d/dy(zeta)-R
!            rhs=res-zeta_rhs(i,1,k)
!
!            ! computes diag. Jacobian of zeta->L.zeta transformation
!            ! by differentiating L.zeta wrt. z(i,1,k) diag. entries
!            Jac=1
!
!            if (action.eq.residual) then
!              zeta_res(i,1,k)=rhs
!            else if (action.eq.lop) then
!              zeta_lop(i,1,k)=res
!            else if (action.eq.relax) then
!              zeta(i,1,k)=zeta(i,1,k)-rhs/Jac                  
!            end if
!
!            norm=norm+rhs**2
!            sum=sum+1
!           end do
!          end do
!        end if

        norm=sqrt(norm/sum)

        return
        end

c-----------------------------------------------------------------------
c The following initializes the rest of the metric and Hb,
c given zeta
c
c NOTE: if we ever add gb_xy,gb_xz,gb_yz, must define them
c       as AMRD MG_cnst vars
c-----------------------------------------------------------------------
        subroutine init_ghb(zeta,gb_tt,gb_tx,gb_ty,
     &                      gb_tz,
     &                      gb_xx,gb_xy,
     &                      gb_xz,
     &                      gb_yy,
     &                      gb_yz,
     &                      gb_zz,Hb_t,Hb_x,Hb_y,
     &                      Hb_z,
     &                      L,cmask,phys_bdy,
     &                      x,y,z,chr,ex,Nx,Ny,Nz,regtype,rhoa,rhob)
        implicit none
        integer Nx,Ny,Nz
        real*8 zeta(Nx,Ny,Nz)
        real*8 gb_xx(Nx,Ny,Nz),gb_tt(Nx,Ny,Nz)
        real*8 gb_yy(Nx,Ny,Nz),gb_tx(Nx,Ny,Nz),gb_zz(Nx,Ny,Nz)
        real*8 gb_ty(Nx,Ny,Nz)
        real*8 gb_tz(Nx,Ny,Nz)
        real*8 gb_xy(Nx,Ny,Nz)
        real*8 gb_xz(Nx,Ny,Nz)
        real*8 gb_yz(Nx,Ny,Nz)
        real*8 Hb_t(Nx,Ny,Nz)
        real*8 Hb_x(Nx,Ny,Nz)
        real*8 Hb_y(Nx,Ny,Nz)
        real*8 Hb_z(Nx,Ny,Nz)
        real*8 cmask(Nx,Ny,Nz),chr(Nx,Ny,Nz),ex,L
        real*8 x(Nx),y(Ny),z(Nz)
        real*8 rhoa,rhob
        integer phys_bdy(6),regtype

        real*8 zeros(Nx,Ny,Nz)

        real*8 zeta0
        
        real*8 PI
        parameter (PI=0.3141592653589793D1)
        integer i,j,k
        real*8 x0,y0,z0
        real*8 rho0
        real*8 f0

        real*8 tfunction(Nx,Ny,Nz)

        real*8 f1,trans

        real*8 g0_tt_ads0,g0_xx_ads0
        real*8 g0_tx_ads0,g0_ty_ads0,g0_tz_ads0
        real*8 g0_xy_ads0,g0_yy_ads0,g0_zz_ads0
        real*8 g0_xz_ads0,g0_yz_ads0

        !--------------------------------------------------------------

        ! initialize metric given zeta, using full metric expression
        ! g0_ij=g0_ij_ads*zeta^2, and 
        ! g-bar expression g0_ij=g0_ij_ads+gb_ij
        do i=2,Nx-1
          do j=2,Ny-1
           do k=2,Nz-1
 
            zeros(i,j,k)=0

            x0=x(i)
            y0=y(j)
            z0=z(k)
            rho0=sqrt(x0**2+y0**2+z0**2)


!!!!!!!!!2+1 VERSION!!!!!!!!
!            f0=(1-rho0**2)**2+4*rho0**2/L**2
!
!            ! set gads values
!            g0_tt_ads0 =-f0
!     &                 /(1-rho0**2)**2
!            g0_tx_ads0 =0
!            g0_ty_ads0 =0
!            g0_tz_ads0 =0
!            g0_xx_ads0 =(x0**2*(1+rho0**2)**2/f0+y0**2)
!     &                 /(1-rho0**2)**2
!     &                 /rho0**2
!     &                 *4
!            g0_xy_ads0 =((1+rho0**2)**2/f0-1)
!     &                 /(1-rho0**2)**2
!     &                 /rho0**2
!     &                 *x0*y0
!     &                 *4
!            g0_xz_ads0 =0
!            g0_yy_ads0 =(y0**2*(1+rho0**2)**2/f0+x0**2)
!     &                 /(1-rho0**2)**2
!     &                 /rho0**2
!     &                 *4
!            g0_yz_ads0 =0
!            g0_psi_ads0=(y0**2)
!     &                 /(1-rho0**2)**2
!     &                 *4
!
!!!!!!!!!!!!!!!!!!!!!!!!

!3+1 version
            ! set gads values
            g0_tt_ads0 =-(4*rho0**2+L**2*(-1+rho0**2)**2)
     &               /L**2/(-1+rho0**2)**2
            g0_tx_ads0 =0
            g0_ty_ads0 =0
            g0_tz_ads0 =0
            g0_xx_ads0 =(8*(-1+L**2)*(x0**2-y0**2-z0**2)
     &              +8*rho0**2+4*L**2*(1+rho0**4))
     &              /(-1+rho0**2)**2/(4*rho0**2+L**2*(-1+rho0**2)**2)
            g0_xy_ads0 =(16 *(-1 + L**2) *x0* y0)
     &              /((-1 + rho0**2)**2 
     &               *(4 *rho0**2 +L**2 *(-1 +rho0**2)**2))
            g0_xz_ads0 =(16 *(-1 + L**2) *x0* z0)
     &              /((-1 + rho0**2)**2
     &               *(4 *rho0**2 +L**2 *(-1 +rho0**2)**2))
            g0_yy_ads0 =(4*(4*(x0**2+z0**2)+L**2*(x0**4+(1+y0**2)**2
     &              +2*(-1+y0**2)*z0**2+z0**4
     &              +2*x0**2*(-1+y0**2+z0**2))))
     &              /(L**2*(-1+rho0**2)**4+4*(-1+rho0**2)**2*(rho0**2))
            g0_yz_ads0 =(16 *(-1 + L**2) *y0* z0)
     &              /((-1 + rho0**2)**2
     &               *(4 *rho0**2 +L**2 *(-1 +rho0**2)**2))
            g0_zz_ads0=(4*(4*(x0**2+y0**2)+L**2*((-1+x0**2+y0**2)**2
     &              +2*(1+x0**2+y0**2)*z0**2+z0**4)))
     &              /(L**2*(-1+rho0**2)**4
     &              +4*(-1+rho0**2)**2*(rho0**2))

            if (chr(i,j,k).ne.ex) then 

              zeta0=zeta(i,j,k)

              gb_tx(i,j,k)=0
              gb_ty(i,j,k)=0
              gb_tz(i,j,k)=0
              gb_xx(i,j,k)=g0_xx_ads0*(zeta0**4-1)
              gb_xy(i,j,k)=g0_xy_ads0*(zeta0**4-1)
              gb_xz(i,j,k)=g0_xz_ads0*(zeta0**4-1)
              gb_yy(i,j,k)=g0_yy_ads0*(zeta0**4-1)
              gb_yz(i,j,k)=g0_yz_ads0*(zeta0**4-1)
              gb_zz(i,j,k)=g0_zz_ads0*(zeta0**4-1)

              f1=trans(rho0,rhoa,rhob)

              !consistent with target gauge -gbtt+gbxx+gbyy+gbgb_zz=0
              gb_tt(i,j,k)=(gb_xx(i,j,k)+gb_yy(i,j,k)+gb_zz(i,j,k))*f1

            endif
           end do
          end do
        end do

!        ! y=0 axis regularization
!        call axi_reg_g(gb_tt,gb_tx,gb_ty,
!     &                 gb_xx,gb_xy,gb_yy,psi,tfunction,chr,ex,
!     &                 L,x,y,z,Nx,Ny,Nz,regtype)

        return
        end

c----------------------------------------------------------------------
c the following computes ddf, the background Laplacian on f, and 
c grad_f_sq, the background squared gradient of f; 
c evaluated at point i,j and at the initial time.
c----------------------------------------------------------------------
        subroutine df_int(f0,f0_x,ddf,ddf_Jac,grad_f_sq,
     &                    x,y,z,i,j,k,chr,L,ex,Nx,Ny,Nz)
        implicit none
        integer Nx,Ny,Nz
        real*8 chr(Nx,Ny,Nz),ex
        real*8 x(Nx),y(Ny),z(Nz),L

        real*8 f0(Nx,Ny,Nz),f0_x(4),f0_xx(4,4),ddf,ddf_Jac,grad_f_sq

        real*8 PI
        parameter (PI=0.3141592653589793D1)

        integer i,j,k,i1,j1,a,b,c,d

        real*8 dx,dy,dz,dt
        real*8 x0,y0,z0,rho0
        real*8 x2,y2,x3,y3,x4,y4

        dx=(x(2)-x(1))
        dy=(y(2)-y(1))
        dz=(z(2)-z(1))
        dt=1.0d0 ! only a placeholder for df2_int (no time derivatives needed in MG ID-solver)

        ! sets x0,y0,rho0 to values at i,j
        x0=x(i)
        y0=y(j)
        z0=z(k)
        rho0=sqrt(x0**2+y0**2+z0**2)

        x2=x0*x0
        y2=y0*y0
        x3=x0*x0*x0
        y3=y0*y0*y0
        x4=x0*x0*x0*x0
        y4=y0*y0*y0*y0

!debugging!
!        do a=1,Nx
!         do b=1,Ny
!          do c=1,Nz
!           f0(a,b,c)=x(a)**2+y(b)**3+x(a)*y(b)**4*z(c)**7
!          end do
!         end do
!        end do
!
!       write(*,*) "L,x0,y0,z0=",L,x0,y0,z0,dx
!       write(*,*) "f0(i,j,k)=",f0(i,j,k)
!       write(*,*) "f0(i+1,j,k)=",f0(i+1,j,k)
!       write(*,*) "f0(i-1,j,k)=",f0(i-1,j,k)
!       write(*,*) "(f0(i+1,j,k)-f0(i-1,j,k))/2/dx="
!     &             ,(f0(i+1,j,k)-f0(i-1,j,k))/2/dx

        ! set first and second derivatives
        !(time symmetric initial data for now, so time derivatives vanish, and
        ! S2 symmetric initial data for now, so phi, theta derivatives vanish) 
        f0_x(1)=0
        f0_x(2)=(f0(i+1,j,k)-f0(i-1,j,k))/2/dx!f_x
        f0_x(3)=(f0(i,j+1,k)-f0(i,j-1,k))/2/dy!f_y
        f0_x(4)=(f0(i,j,k+1)-f0(i,j,k-1))/2/dz!f_z

        f0_xx(1,1)=0
        f0_xx(1,2)=0
        f0_xx(1,3)=0
        f0_xx(1,4)=0
        f0_xx(2,2)=(f0(i+1,j,k)-2*f0(i,j,k)+f0(i-1,j,k))/dx/dx
        f0_xx(2,3)=( (f0(i+1,j+1,k)-f0(i+1,j-1,k))/2/dy
     &              -(f0(i-1,j+1,k)-f0(i-1,j-1,k))/2/dy )/2/dx
        f0_xx(2,4)=( (f0(i+1,j,k+1)-f0(i+1,j,k-1))/2/dz
     &              -(f0(i-1,j,k+1)-f0(i-1,j,k-1))/2/dz )/2/dx
        f0_xx(3,3)=(f0(i,j+1,k)-2*f0(i,j,k)+f0(i,j-1,k))/dy/dy 
        f0_xx(3,4)=( (f0(i,j+1,k+1)-f0(i,j+1,k-1))/2/dz
     &              -(f0(i,j-1,k+1)-f0(i,j-1,k-1))/2/dz )/2/dy
        f0_xx(4,4)=(f0(i,j,k+1)-2*f0(i,j,k)+f0(i,j,k-1))/dz/dz

!       write(*,*) "L,x0,y0,z0=",L,x0,y0,z0
!       write(*,*) "f0_x(2),f0_x(3),f0_x(4)=",f0_x(2),f0_x(3),f0_x(4)
!       write(*,*) "f0_xx(2,2),f0_xx(2,3),f0_xx(2,4)="
!     &             ,f0_xx(2,2),f0_xx(2,3),f0_xx(2,4)
!       write(*,*) "f0_xx(3,3),f0_xx(3,4)="
!     &             ,f0_xx(3,3),f0_xx(3,4)
!       write(*,*) "f0_xx(4,4)="
!     &             ,f0_xx(4,4)

        do a=1,3
          do b=a+1,4
            f0_xx(b,a)=f0_xx(a,b)
          end do
        end do

!!!!2+1 version!!!!!!!!
!        ! calculate ddf, background Laplacian acting on f
!        ! and ddf_Jac, Jacobian of f->ddf transformation 
!        !(DDf = g^ab D_a D_b f)
!        ddf= 
!     &        f0_x(3)*
!     &        ((-1 + x0**2 - y0**2)*(-1 + x0**2 + y0**2))/(4*y0)
!     &       +f0_xx(3,3)*
!     &        (-1 + x0**2 + y0**2)**2/4
!     &       +f0_x(2)*
!     &        (-(x0*(-1 + x0**2 + y0**2)))/2
!     &       +f0_xx(2,3)*
!     &        0.0d0
!     &       +f0_xx(2,2)*
!     &        (-1 + x0**2 + y0**2)**2/4
!
!        ddf_Jac= 
!     &            (-2/dy/dy)*
!     &            (-1 + x0**2 + y0**2)**2/4
!     &           +(-2/dx/dx)*
!     &            (-1 + x0**2 + y0**2)**2/4
!
!        ! calculate grad_f_sq, squared gradient of f 
!        !(Df^2 = g^ab D_a f D_b f)
!        grad_f_sq= 
!     &              f0_x(3)*f0_x(3)*
!     &              (-1 + x0**2 + y0**2)**2/4
!     &             +2*f0_x(3)*f0_x(2)*
!     &              0.0d0
!     &             +f0_x(2)*f0_x(2)*
!     &              (-1 + x0**2 + y0**2)**2/4
!!!!!!!!!!!!!!!!!!!!!!!!

!3+1version

        ddf= 
     &      (((-1+rho0**2)**2*(4*x0**2+L**2*(1+(-2+x0)*x0+y0**2+z0**2)
     &       *(1+x0*(2+x0)+y0**2+z0**2)))
     &      /(4*L**2*(1+rho0**2)**2))
     &      *f0_xx(2,2)
     &     +2*(-(((-1+L**2)*x0*y0*(-1+rho0**2)**2)
     &      /(L**2*(1+rho0**2)**2)) )
     &      *f0_xx(2,3)
     &     +2*(-(((-1+L**2)*x0*z0*(-1+rho0**2)**2)
     &      /(L**2*(1+rho0**2)**2)) )
     &      *f0_xx(2,4)
     &     +(((-1+rho0**2)**2*(4*y0**2+L**2*(x0**2+(-1+y0)**2+z0**2)
     &      *(x0**2+(1+y0)**2+z0**2)))
     &      /(4*L**2*(1+rho0**2)**2))
     &      *f0_xx(3,3)
     &     +2*(-(((-1+L**2)*y0*z0*(-1+rho0**2)**2)
     &      /(L**2*(1+rho0**2)**2)) )
     &      *f0_xx(3,4)
     &     +(((-1+rho0**2)**2*(4*z0**2+L**2*((1+x0**2+y0**2)**2
     &      +2*(-1+x0**2+y0**2)*z0**2+z0**4)))
     &      /(4*L**2*(1+rho0**2)**2))
     &      *f0_xx(4,4)
     &   
     &     +(-(1/(2*L**2*(1+rho0**2)**3))*x0
     &      *(-1+rho0**2)*(6+2*(rho0**2)**2
     &      +L**2*(-1+rho0**2)*(5+x0**4+2*y0**2+2*z0**2
     &      +(y0**2+z0**2)**2+2*x0**2*(1+y0**2+z0**2))))
     &      *f0_x(2)
     &     +(-(1/(2*L**2*(1+rho0**2)**3))*y0
     &      *(-1+rho0**2)*(6+2*(rho0**2)**2
     &      +L**2*(-1+rho0**2)*(5+x0**4+2*y0**2+2*z0**2
     &      +(y0**2+z0**2)**2+2*x0**2*(1+y0**2+z0**2))))
     &      *f0_x(3)
     &     +(-(1/(2*L**2*(1+rho0**2)**3))*z0
     &      *(-1+rho0**2)*(6+2*(rho0**2)**2
     &      +L**2*(-1+rho0**2)*(5+x0**4+2*y0**2+2*z0**2
     &      +(y0**2+z0**2)**2+2*x0**2*(1+y0**2+z0**2))))
     &      *f0_x(4)

        ddf_Jac= 
     &      (((-1+rho0**2)**2*(4*x0**2+L**2*(1+(-2+x0)*x0+y0**2+z0**2)
     &       *(1+x0*(2+x0)+y0**2+z0**2)))
     &      /(4*L**2*(1+rho0**2)**2))
     &      *(-2/dx/dx)
     &     +(((-1+rho0**2)**2*(4*y0**2+L**2*(x0**2+(-1+y0)**2+z0**2)
     &      *(x0**2+(1+y0)**2+z0**2)))
     &      /(4*L**2*(1+rho0**2)**2))
     &      *(-2/dy/dy)
     &     +(((-1+rho0**2)**2*(4*z0**2+L**2*((1+x0**2+y0**2)**2
     &      +2*(-1+x0**2+y0**2)*z0**2+z0**4)))
     &      /(4*L**2*(1+rho0**2)**2))
     &      *(-2/dz/dz)
            
        grad_f_sq= 
     &       (((-1+rho0**2)**2*(4*x0**2+L**2*(1+(-2+x0)*x0+y0**2+z0**2)
     &       *(1+x0*(2+x0)+y0**2+z0**2)))
     &      /(4*L**2*(1+rho0**2)**2))
     &      *f0_x(2)*f0_x(2)
     &     +2*(-(((-1+L**2)*x0*y0*(-1+rho0**2)**2)
     &      /(L**2*(1+rho0**2)**2)) )
     &      *f0_x(2)*f0_x(3)
     &     +2*(-(((-1+L**2)*x0*z0*(-1+rho0**2)**2)
     &      /(L**2*(1+rho0**2)**2)) )
     &      *f0_x(2)*f0_x(4)
     &     +(((-1+rho0**2)**2*(4*y0**2+L**2*(x0**2+(-1+y0)**2+z0**2)
     &      *(x0**2+(1+y0)**2+z0**2)))
     &      /(4*L**2*(1+rho0**2)**2))
     &      *f0_x(3)*f0_x(3)
     &     +2*(-(((-1+L**2)*y0*z0*(-1+rho0**2)**2)
     &      /(L**2*(1+rho0**2)**2)) )
     &      *f0_x(3)*f0_x(4)
     &     +(((-1+rho0**2)**2*(4*z0**2+L**2*((1+x0**2+y0**2)**2
     &      +2*(-1+x0**2+y0**2)*z0**2+z0**4)))
     &      /(4*L**2*(1+rho0**2)**2))
     &      *f0_x(4)*f0_x(4)

!       write(*,*) "L,x0,y0,z0=",L,x0,y0,z0
!       write(*,*) "ddf,grad_f_sq=",ddf,grad_f_sq

        return
        
        end

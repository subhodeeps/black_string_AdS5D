c----------------------------------------------------------------------
c NOT NEEDED FOR GENERAL STUDIES SUCH AS 3+1 Cartesian in 4 dimensions
c in cartesian coordinates t,x,y for x in [-1,1], y in [0,1]
c
c routines for computing y=0 axis regularity conditions for
c metric and scalar fied
c----------------------------------------------------------------------
        subroutine axi_reg_g(gb_tt,gb_tx,gb_ty,
     &                       gb_xx,gb_xy,gb_yy,psi,tfunction,chr,ex,
     &                       L,x,y,z,Nx,Ny,Nz,regtype)
        implicit none
        integer Nx,Ny,Nz
        integer regtype
        real*8 gb_tt(Nx,Ny,Nz),gb_tx(Nx,Ny,Nz),gb_ty(Nx,Ny,Nz)
        real*8 gb_xx(Nx,Ny,Nz),gb_xy(Nx,Ny,Nz),gb_yy(Nx,Ny,Nz)
        real*8 psi(Nx,Ny,Nz),tfunction(Nx,Ny,Nz),chr(Nx,Ny,Nz),ex,L
        real*8 x(Nx),y(Ny),z(Nz)

        integer i,k
        real*8 dx,dy,dz
        real*8 PI
        parameter (PI=3.141592653589793d0)

        dx=x(2)-x(1)
        dy=y(2)-y(1)
        dz=z(2)-z(1)

        if (abs(y(1)).gt.dy/2) return
 
        do i=1,Nx
         do k=1,Nz
           if (chr(i,1,k).ne.ex) then

              ! 1-pt axis bcs, using 3-pt zero derivative to get axis pt; 
              ! set psi=gb_yy at y=0
              if (regtype.eq.1) then 
                gb_tt(i,1,k)=(4*gb_tt(i,2,k)-gb_tt(i,3,k))/3
                gb_tx(i,1,k)=(4*gb_tx(i,2,k)-gb_tx(i,3,k))/3
                gb_xx(i,1,k)=(4*gb_xx(i,2,k)-gb_xx(i,3,k))/3
                gb_yy(i,1,k)=(4*gb_yy(i,2,k)-gb_yy(i,3,k))/3
                psi(i,1,k)=gb_yy(i,1,k)
                gb_xy(i,1,k)=0
                gb_ty(i,1,k)=0
              ! 1-pt axis bcs, using 3-pt zero derivative to get axis pt; 
              ! set gb_yy=psi at y=0
              else if (regtype.eq.2) then 
                gb_tt(i,1,k)=(4*gb_tt(i,2,k)-gb_tt(i,3,k))/3
                gb_tx(i,1,k)=(4*gb_tx(i,2,k)-gb_tx(i,3,k))/3
                gb_xx(i,1,k)=(4*gb_xx(i,2,k)-gb_xx(i,3,k))/3
                psi(i,1,k)=(4*psi(i,2,k)-psi(i,3,k))/3
                gb_yy(i,1,k)=psi(i,1,k)
                gb_xy(i,1,k)=0
                gb_ty(i,1,k)=0
              !(experimental 2-pt axis bcs)!
              else if (regtype.eq.3) then  
                gb_tt(i,1,k)=(4*gb_tt(i,3,k)-gb_tt(i,5,k))/3
                gb_tx(i,1,k)=(4*gb_tx(i,3,k)-gb_tx(i,5,k))/3
                gb_xx(i,1,k)=(4*gb_xx(i,3,k)-gb_xx(i,5,k))/3
                gb_yy(i,1,k)=(4*gb_yy(i,3,k)-gb_yy(i,5,k))/3
                psi(i,1,k)=gb_yy(i,1,k)
                gb_xy(i,1,k)=0
                gb_ty(i,1,k)=0
                gb_tt(i,2,k)= gb_tt(i,1,k)/4 + 3*gb_tt(i,3,k)/2
     &                     -gb_tt(i,4,k)   +   gb_tt(i,5,k)/4
                gb_tx(i,2,k)= gb_tx(i,1,k)/4 + 3*gb_tx(i,3,k)/2
     &                     -gb_tx(i,4,k)   +   gb_tx(i,5,k)/4
                gb_xx(i,2,k)= gb_xx(i,1,k)/4 + 3*gb_xx(i,3,k)/2
     &                     -gb_xx(i,4,k)   +   gb_xx(i,5,k)/4
                gb_yy(i,2,k)= gb_yy(i,1,k)/4 + 3*gb_yy(i,3,k)/2
     &                     -gb_yy(i,4,k)   +   gb_yy(i,5,k)/4
                psi(i,2,k)  = psi(i,1,k)/4   + 3*psi(i,3,k)/2
     &                     -psi(i,4,k)     +   psi(i,5,k)/4
                gb_xy(i,2,k)=0.5*gb_xy(i,3,k)
                gb_ty(i,2,k)=0.5*gb_ty(i,3,k)
              ! 2-pt axis bcs, using 5-pt zero derivative to get non-axis pt; 
              ! set psi=gb_yy at y=0
              else if (regtype.eq.4) then  
                gb_tt(i,1,k)=(4*gb_tt(i,3,k)-gb_tt(i,5,k))/3
                gb_tx(i,1,k)=(4*gb_tx(i,3,k)-gb_tx(i,5,k))/3
                gb_xx(i,1,k)=(4*gb_xx(i,3,k)-gb_xx(i,5,k))/3
                gb_yy(i,1,k)=(4*gb_yy(i,3,k)-gb_yy(i,5,k))/3
                psi(i,1,k)=gb_yy(i,1,k)
                gb_xy(i,1,k)=0
                gb_ty(i,1,k)=0
                gb_tt(i,2,k)=(
     &            25*gb_tt(i,1,k)+36*gb_tt(i,3,k)
     &            -16*gb_tt(i,4,k)+3*gb_tt(i,5,k)
     &            )/48
                gb_tx(i,2,k)=(
     &            25*gb_tx(i,1,k)+36*gb_tx(i,3,k)
     &            -16*gb_tx(i,4,k)+3*gb_tx(i,5,k)
     &            )/48
                gb_xx(i,2,k)=(
     &            25*gb_xx(i,1,k)+36*gb_xx(i,3,k)
     &            -16*gb_xx(i,4,k)+3*gb_xx(i,5,k)
     &            )/48
                gb_yy(i,2,k)=(
     &            25*gb_yy(i,1,k)+36*gb_yy(i,3,k)
     &            -16*gb_yy(i,4,k)+3*gb_yy(i,5,k)
     &            )/48
                psi(i,2,k)=(
     &            25*psi(i,1,k)+36*psi(i,3,k)
     &            -16*psi(i,4,k)+3*psi(i,5,k)
     &            )/48
                gb_xy(i,2,k)=0.5*gb_xy(i,3,k)
                gb_ty(i,2,k)=0.5*gb_ty(i,3,k)
              ! 2-pt axis bcs, using 5-pt zero derivative to get non-axis pt; 
              ! set gb_yy=psi at y=0
              else if (regtype.eq.5) then  
                gb_tt(i,1,k)=(4*gb_tt(i,3,k)-gb_tt(i,5,k))/3
                gb_tx(i,1,k)=(4*gb_tx(i,3,k)-gb_tx(i,5,k))/3
                gb_xx(i,1,k)=(4*gb_xx(i,3,k)-gb_xx(i,5,k))/3
                gb_yy(i,1,k)=(4*gb_yy(i,3,k)-gb_yy(i,5,k))/3
                gb_yy(i,1,k)=psi(i,1,k)
                gb_xy(i,1,k)=0
                gb_ty(i,1,k)=0
                gb_tt(i,2,k)=(
     &            25*gb_tt(i,1,k)+36*gb_tt(i,3,k)
     &            -16*gb_tt(i,4,k)+3*gb_tt(i,5,k)
     &            )/48
                gb_tx(i,2,k)=(
     &            25*gb_tx(i,1,k)+36*gb_tx(i,3,k)
     &            -16*gb_tx(i,4,k)+3*gb_tx(i,5,k)
     &            )/48
                gb_xx(i,2,k)=(
     &            25*gb_xx(i,1,k)+36*gb_xx(i,3,k)
     &            -16*gb_xx(i,4,k)+3*gb_xx(i,5,k)
     &            )/48
                gb_yy(i,2,k)=(
     &            25*gb_yy(i,1,k)+36*gb_yy(i,3,k)
     &            -16*gb_yy(i,4,k)+3*gb_yy(i,5,k)
     &            )/48
                psi(i,2,k)=(
     &            25*psi(i,1,k)+36*psi(i,3,k)
     &            -16*psi(i,4,k)+3*psi(i,5,k)
     &            )/48
                gb_xy(i,2,k)=0.5*gb_xy(i,3,k)
                gb_ty(i,2,k)=0.5*gb_ty(i,3,k)
              ! 3-pt axis bcs, using 5-pt zero derivative to get non-axis pts; 
              ! set psi=gb_yy at y=0
              else if (regtype.eq.6) then  
                gb_tt(i,1,k)=(4*gb_tt(i,4,k)-gb_tt(i,7,k))/3
                gb_tx(i,1,k)=(4*gb_tx(i,4,k)-gb_tx(i,7,k))/3
                gb_xx(i,1,k)=(4*gb_xx(i,4,k)-gb_xx(i,7,k))/3
                gb_yy(i,1,k)=(4*gb_yy(i,4,k)-gb_yy(i,7,k))/3
                psi(i,1,k)=gb_yy(i,1,k)
                gb_xy(i,1,k)=0
                gb_ty(i,1,k)=0
                gb_tt(i,2,k)=( 107*gb_tt(i,1,k) + 100*gb_tt(i,4,k)
     &                       -75*gb_tt(i,5,k) + 18*gb_tt(i,6,k) )/150
                gb_tx(i,2,k)=( 107*gb_tx(i,1,k) + 100*gb_tx(i,4,k)
     &                       -75*gb_tx(i,5,k) + 18*gb_tx(i,6,k) )/150
                gb_xx(i,2,k)=( 107*gb_xx(i,1,k) + 100*gb_xx(i,4,k)
     &                       -75*gb_xx(i,5,k) + 18*gb_xx(i,6,k) )/150
                gb_yy(i,2,k)=( 107*gb_yy(i,1,k) + 100*gb_yy(i,4,k)
     &                       -75*gb_yy(i,5,k) + 18*gb_yy(i,6,k) )/150
                psi(i,2,k)  =( 107*psi(i,1,k)   + 100*psi(i,4,k)
     &                       -75*psi(i,5,k)   + 18*psi(i,6,k) )/150
                gb_xy(i,2,k)=gb_xy(i,4,k)/3
                gb_ty(i,2,k)=gb_ty(i,4,k)/3
                gb_tt(i,3,k)=( 77*gb_tt(i,1,k)   + 400*gb_tt(i,4,k)
     &                       -225*gb_tt(i,5,k) + 48*gb_tt(i,6,k) )/300
                gb_tx(i,3,k)=( 77*gb_tx(i,1,k)   + 400*gb_tx(i,4,k)
     &                       -225*gb_tx(i,5,k) + 48*gb_tx(i,6,k) )/300
                gb_xx(i,3,k)=( 77*gb_xx(i,1,k)   + 400*gb_xx(i,4,k)
     &                       -225*gb_xx(i,5,k) + 48*gb_xx(i,6,k) )/300
                gb_yy(i,3,k)=( 77*gb_yy(i,1,k)   + 400*gb_yy(i,4,k)
     &                       -225*gb_yy(i,5,k) + 48*gb_yy(i,6,k) )/300
                psi(i,3,k)  =( 77*psi(i,1,k)     + 400*psi(i,4,k)
     &                       -225*psi(i,5,k)   + 48*psi(i,6,k) )/300
                gb_xy(i,3,k)=2*gb_xy(i,4,k)/3
                gb_ty(i,3,k)=2*gb_ty(i,4,k)/3
              ! 3-pt axis bcs, using 5-pt zero derivative to get non-axis pts; 
              ! set gb_yy=psi at y=0
              else if (regtype.eq.7) then  
                gb_tt(i,1,k)=(4*gb_tt(i,4,k)-gb_tt(i,7,k))/3
                gb_tx(i,1,k)=(4*gb_tx(i,4,k)-gb_tx(i,7,k))/3
                gb_xx(i,1,k)=(4*gb_xx(i,4,k)-gb_xx(i,7,k))/3
                psi(i,1,k)=(4*psi(i,4,k)-psi(i,7,k))/3
                gb_yy(i,1,k)=psi(i,1,k)
                gb_xy(i,1,k)=0
                gb_ty(i,1,k)=0
                gb_tt(i,2,k)=( 107*gb_tt(i,1,k) + 100*gb_tt(i,4,k)
     &                       -75*gb_tt(i,5,k) + 18*gb_tt(i,6,k) )/150
                gb_tx(i,2,k)=( 107*gb_tx(i,1,k) + 100*gb_tx(i,4,k)
     &                       -75*gb_tx(i,5,k) + 18*gb_tx(i,6,k) )/150
                gb_xx(i,2,k)=( 107*gb_xx(i,1,k) + 100*gb_xx(i,4,k)
     &                       -75*gb_xx(i,5,k) + 18*gb_xx(i,6,k) )/150
                gb_yy(i,2,k)=( 107*gb_yy(i,1,k) + 100*gb_yy(i,4,k)
     &                       -75*gb_yy(i,5,k) + 18*gb_yy(i,6,k) )/150
                psi(i,2,k)  =( 107*psi(i,1,k)   + 100*psi(i,4,k)
     &                       -75*psi(i,5,k)   + 18*psi(i,6,k) )/150
                gb_xy(i,2,k)=gb_xy(i,4,k)/3
                gb_ty(i,2,k)=gb_ty(i,4,k)/3
                gb_tt(i,3,k)=( 77*gb_tt(i,1,k)   + 400*gb_tt(i,4,k)
     &                       -225*gb_tt(i,5,k) + 48*gb_tt(i,6,k) )/300
                gb_tx(i,3,k)=( 77*gb_tx(i,1,k)   + 400*gb_tx(i,4,k)
     &                       -225*gb_tx(i,5,k) + 48*gb_tx(i,6,k) )/300
                gb_xx(i,3,k)=( 77*gb_xx(i,1,k)   + 400*gb_xx(i,4,k)
     &                       -225*gb_xx(i,5,k) + 48*gb_xx(i,6,k) )/300
                gb_yy(i,3,k)=( 77*gb_yy(i,1,k)   + 400*gb_yy(i,4,k)
     &                       -225*gb_yy(i,5,k) + 48*gb_yy(i,6,k) )/300
                psi(i,3,k)  =( 77*psi(i,1,k)     + 400*psi(i,4,k)
     &                       -225*psi(i,5,k)   + 48*psi(i,6,k) )/300
                gb_xy(i,3,k)=2*gb_xy(i,4,k)/3
                gb_ty(i,3,k)=2*gb_ty(i,4,k)/3
              else
                write(*,*) 'axi_reg_g error: invalid regtype'
                write(*,*) 'regtype=',regtype
                stop
              endif

           else

              gb_tt(i,1,k)=0 
              gb_tx(i,1,k)=0
              gb_xx(i,1,k)=0
              gb_yy(i,1,k)=0
              gb_xy(i,1,k)=0
              gb_ty(i,1,k)=0
              psi(i,1,k)=0

           end if
         end do
        end do

        return
        end

c----------------------------------------------------------------------
        subroutine axi_reg_phi(phi,chr,ex,L,x,y,z,Nx,Ny,Nz,regtype)
        implicit none
        integer Nx,Ny,Nz
        real*8 phi(Nx,Ny,Nz),chr(Nx,Ny,Nz),ex,L
        real*8 x(Nx),y(Ny),z(Nz)

        real*8 PI,dx,dy,dz
        parameter (PI=3.141592653589793d0)
        integer i,k

        integer regtype

        dx=x(2)-x(1)
        dy=y(2)-y(1)
        dz=z(2)-z(1)

        if (abs(y(1)).gt.dy/2) return

        do i=1,Nx
         do k=1,Nz
           if (chr(i,1,k).ne.ex) then

              if (regtype.eq.1) then
                phi(i,1,k)=(4*phi(i,2,k)-phi(i,3,k))/3
              else if (regtype.eq.2) then
                phi(i,1,k)=(4*phi(i,2,k)-phi(i,3,k))/3
              else if (regtype.eq.3) then
                phi(i,1,k)=(4*phi(i,3,k)-phi(i,5,k))/3
                phi(i,2,k)  = phi(i,1,k)/4   + 3*phi(i,3,k)/2
     &                     -phi(i,4,k)     +   phi(i,5,k)/4
              else if (regtype.eq.4) then
                phi(i,1,k)=(4*phi(i,3,k)-phi(i,5,k))/3
                phi(i,2,k)=(
     &             phi(i,1,k)+6*phi(i,3,k)-4*phi(i,4,k)+phi(i,5,k)
     &             )/4
              else if (regtype.eq.5) then
                phi(i,1,k)=(4*phi(i,3,k)-phi(i,5,k))/3
                phi(i,2,k)=(
     &            25*phi(i,1,k)+36*phi(i,3,k)-16*phi(i,4,k)+3*phi(i,5,k)
     &            )/48
              else if (regtype.eq.6) then
                phi(i,1,k)=(4*phi(i,4,k)-phi(i,7,k))/3
                phi(i,2,k)  =( 107*phi(i,1,k)   + 100*phi(i,4,k)
     &                       -75*phi(i,5,k)   + 18*phi(i,6,k) )/150
                phi(i,3,k)  =( 77*phi(i,1,k)     + 400*phi(i,4,k)
     &                       -225*phi(i,5,k)   + 48*phi(i,6,k) )/300
              else if (regtype.eq.7) then
                phi(i,1,k)=(4*phi(i,4,k)-phi(i,7,k))/3
                phi(i,2,k)  =( 107*phi(i,1,k)   + 100*phi(i,4,k)
     &                       -75*phi(i,5,k)   + 18*phi(i,6,k) )/150
                phi(i,3,k)  =( 77*phi(i,1,k)     + 400*phi(i,4,k)
     &                       -225*phi(i,5,k)   + 48*phi(i,6,k) )/300
              else
                write(*,*) 'axi_reg_phi error: invalid regtype'
                write(*,*) 'regtype=',regtype
                stop
              endif
           
           else 

              phi(i,1,k)=0

           end if
         end do
        end do
           
        return
        end

c----------------------------------------------------------------------
        subroutine axi_reg_Hb(Hb_t,Hb_x,Hb_y,chr,ex,
     &                        L,x,y,z,Nx,Ny,Nz,regtype)
        implicit none
        integer Nx,Ny,Nz
        real*8 Hb_t(Nx,Ny,Nz),Hb_x(Nx,Ny,Nz),Hb_y(Nx,Ny,Nz)
        real*8 gb_xx(Nx,Ny,Nz),gb_xy(Nx,Ny,Nz),gb_yy(Nx,Ny,Nz)
        real*8 psi(Nx,Ny,Nz),chr(Nx,Ny,Nz),ex,L
        real*8 x(Nx),y(Ny),z(Nz)

        real*8 PI,dx,dy,dz
        parameter (PI=3.141592653589793d0)
        integer i,k
 
        integer regtype

        dx=x(2)-x(1)
        dy=y(2)-y(1)
        dz=z(2)-z(1)

        if (abs(y(1)).gt.dy/2) return

        do i=1,Nx
         do k=1,Nz
           if (chr(i,1,k).ne.ex) then

              if (regtype.eq.1) then
                Hb_t(i,1,k)=(4*Hb_t(i,2,k)-Hb_t(i,3,k))/3
                Hb_x(i,1,k)=(4*Hb_x(i,2,k)-Hb_x(i,3,k))/3
                Hb_y(i,1,k)=0 
              else if (regtype.eq.2) then
                Hb_t(i,1,k)=(4*Hb_t(i,2,k)-Hb_t(i,3,k))/3
                Hb_x(i,1,k)=(4*Hb_x(i,2,k)-Hb_x(i,3,k))/3
                Hb_y(i,1,k)=0 
              else if (regtype.eq.3) then
                Hb_t(i,1,k)=(4*Hb_t(i,3,k)-Hb_t(i,5,k))/3
                Hb_x(i,1,k)=(4*Hb_x(i,3,k)-Hb_x(i,5,k))/3
                Hb_y(i,1,k)=0
                Hb_t(i,2,k)  = Hb_t(i,1,k)/4   + 3*Hb_t(i,3,k)/2
     &                     -Hb_t(i,4,k)     +   Hb_t(i,5,k)/4
                Hb_x(i,2,k)  = Hb_x(i,1,k)/4   + 3*Hb_x(i,3,k)/2
     &                     -Hb_x(i,4,k)     +   Hb_x(i,5,k)/4
                Hb_y(i,2,k)=0.5*Hb_y(i,3,k)
              else if (regtype.eq.4) then
                Hb_t(i,1,k)=(4*Hb_t(i,3,k)-Hb_t(i,5,k))/3
                Hb_x(i,1,k)=(4*Hb_x(i,3,k)-Hb_x(i,5,k))/3
                Hb_y(i,1,k)=0
                Hb_t(i,2,k)=(
     &             Hb_t(i,1,k)+6*Hb_t(i,3,k)-4*Hb_t(i,4,k)+Hb_t(i,5,k)
     &             )/4
                Hb_x(i,2,k)=(
     &             Hb_x(i,1,k)+6*Hb_x(i,3,k)-4*Hb_x(i,4,k)+Hb_x(i,5,k)
     &             )/4
                Hb_y(i,2,k)=0.5*Hb_y(i,3,k)
              else if (regtype.eq.5) then
                Hb_t(i,1,k)=(4*Hb_t(i,3,k)-Hb_t(i,5,k))/3
                Hb_x(i,1,k)=(4*Hb_x(i,3,k)-Hb_x(i,5,k))/3
                Hb_y(i,1,k)=0
                Hb_t(i,2,k)=(
     &            25*Hb_t(i,1,k)+36*Hb_t(i,3,k)
     &            -16*Hb_t(i,4,k)+3*Hb_t(i,5,k)
     &            )/48
                Hb_x(i,2,k)=(
     &            25*Hb_x(i,1,k)+36*Hb_x(i,3,k)
     &            -16*Hb_x(i,4,k)+3*Hb_x(i,5,k)
     &            )/48
                Hb_y(i,2,k)=0.5*Hb_y(i,3,k)
              else if (regtype.eq.6) then
                Hb_t(i,1,k)=(4*Hb_t(i,4,k)-Hb_t(i,7,k))/3
                Hb_x(i,1,k)=(4*Hb_x(i,4,k)-Hb_x(i,7,k))/3
                Hb_y(i,1,k)=0
                Hb_t(i,2,k)  =( 107*Hb_t(i,1,k)   + 100*Hb_t(i,4,k)
     &                       -75*Hb_t(i,5,k)   + 18*Hb_t(i,6,k) )/150
                Hb_x(i,2,k)  =( 107*Hb_x(i,1,k)   + 100*Hb_x(i,4,k)
     &                       -75*Hb_x(i,5,k)   + 18*Hb_x(i,6,k) )/150
                Hb_y(i,2,k)=Hb_y(i,4,k)/3
                Hb_t(i,3,k)  =( 77*Hb_t(i,1,k)     + 400*Hb_t(i,4,k)
     &                       -225*Hb_t(i,5,k)   + 48*Hb_t(i,6,k) )/300
                Hb_x(i,3,k)  =( 77*Hb_x(i,1,k)     + 400*Hb_x(i,4,k)
     &                       -225*Hb_x(i,5,k)   + 48*Hb_x(i,6,k) )/300
                Hb_y(i,3,k)=2*Hb_y(i,4,k)/3
              else if (regtype.eq.7) then
                Hb_t(i,1,k)=(4*Hb_t(i,4,k)-Hb_t(i,7,k))/3
                Hb_x(i,1,k)=(4*Hb_x(i,4,k)-Hb_x(i,7,k))/3
                Hb_y(i,1,k)=0
                Hb_t(i,2,k)  =( 107*Hb_t(i,1,k)   + 100*Hb_t(i,4,k)
     &                       -75*Hb_t(i,5,k)   + 18*Hb_t(i,6,k) )/150
                Hb_x(i,2,k)  =( 107*Hb_x(i,1,k)   + 100*Hb_x(i,4,k)
     &                       -75*Hb_x(i,5,k)   + 18*Hb_x(i,6,k) )/150
                Hb_y(i,2,k)=Hb_y(i,4,k)/3
                Hb_t(i,3,k)  =( 77*Hb_t(i,1,k)     + 400*Hb_t(i,4,k)
     &                       -225*Hb_t(i,5,k)   + 48*Hb_t(i,6,k) )/300
                Hb_x(i,3,k)  =( 77*Hb_x(i,1,k)     + 400*Hb_x(i,4,k)
     &                       -225*Hb_x(i,5,k)   + 48*Hb_x(i,6,k) )/300
                Hb_y(i,3,k)=2*Hb_y(i,4,k)/3
              else
                write(*,*) 'axi_reg_Hb error: invalid regtype'
                write(*,*) 'regtype=',regtype
                stop
              endif

           else 

              Hb_t(i,1,k)=0
              Hb_x(i,1,k)=0
              Hb_y(i,1,k)=0

           end if
         end do
        end do
           
        return
        end

c----------------------------------------------------------------------
cTHE FOLLOWING ROUTINES ARE STILL IN THEIR 2+1 VERSION
!        subroutine ref_sym_g(gb_tt,gb_tx,gb_ty,
!     &                       gb_xx,gb_xy,gb_yy,psi,tfunction,chr,ex,
!     &                       L,x,y,z,Nx,Ny,Nz,regtype)
!        implicit none
!        integer Nx,Ny,Nz
!        integer regtype
!        real*8 gb_tt(Nx,Ny),gb_tx(Nx,Ny),gb_ty(Nx,Ny)
!        real*8 gb_xx(Nx,Ny),gb_xy(Nx,Ny),gb_yy(Nx,Ny)
!        real*8 psi(Nx,Ny),tfunction(Nx,Ny),chr(Nx,Ny),ex,L
!        real*8 x(Nx),y(Ny),z(Nz)
!
!        integer j
!        real*8 dx,dy,dz
!        real*8 PI
!        parameter (PI=3.141592653589793d0)
!
!        dx=x(2)-x(1)
!        dy=y(2)-y(1)
!        dz=z(2)-z(1)
!
!        if (abs(x(1)).gt.dx/2) return
! 
!        do j=1,Ny
!           if (chr(1,j).ne.ex) then
!
!              gb_tt(1,j)=(4*gb_tt(2,j)-gb_tt(3,j))/3
!              gb_tx(1,j)=(4*gb_tx(2,j)-gb_tx(3,j))/3
!              gb_ty(1,j)=(4*gb_ty(2,j)-gb_ty(3,j))/3
!              gb_xx(1,j)=(4*gb_xx(2,j)-gb_xx(3,j))/3
!              gb_xy(1,j)=(4*gb_xy(2,j)-gb_xy(3,j))/3
!              gb_yy(1,j)=(4*gb_yy(2,j)-gb_yy(3,j))/3
!              psi(1,j)=(4*psi(2,j)-psi(3,j))/3
!
!           end if
!        end do
!
!        return
!        end
!
!c----------------------------------------------------------------------
!        subroutine ref_sym_phi(phi,chr,ex,L,x,y,z,Nx,Ny,Nz,regtype)
!        implicit none
!        integer Nx,Ny,Nz
!        real*8 phi(Nx,Ny),chr(Nx,Ny),ex,L
!        real*8 x(Nx),y(Ny),z(Nz)
!
!        real*8 PI,dx,dy,dz
!        parameter (PI=3.141592653589793d0)
!        integer j
!
!        integer regtype
!
!        dx=x(2)-x(1)
!        dy=y(2)-y(1)
!        dz=z(2)-z(1)
!
!        if (abs(x(1)).gt.dx/2) return
!
!        do j=1,Ny
!           if (chr(1,j).ne.ex) then
!
!              phi(1,j)=(4*phi(2,j)-phi(3,j))/3
!
!           end if
!        end do
!
!        return
!        end
!
!c----------------------------------------------------------------------
!        subroutine ref_sym_Hb(Hb_t,Hb_x,Hb_y,chr,ex,
!     &                        L,x,y,z,Nx,Ny,Nz,regtype)
!        implicit none
!        integer Nx,Ny,Nz
!        real*8 Hb_t(Nx,Ny),Hb_x(Nx,Ny),Hb_y(Nx,Ny)
!        real*8 gb_xx(Nx,Ny),gb_xy(Nx,Ny),gb_yy(Nx,Ny)
!        real*8 psi(Nx,Ny),chr(Nx,Ny),ex,L
!        real*8 x(Nx),y(Ny),z(Nz)
!
!        real*8 PI,dx,dy,dz
!        parameter (PI=3.141592653589793d0)
!        integer j
!
!        integer regtype
!
!        dx=x(2)-x(1)
!        dy=y(2)-y(1)
!        dz=z(2)-z(1)
!
!        if (abs(x(1)).gt.dx/2) return
!
!        do j=1,Ny
!           if (chr(1,j).ne.ex) then
!
!              Hb_t(1,j)=(4*Hb_t(2,j)-Hb_t(3,j))/3
!              Hb_x(1,j)=(4*Hb_x(2,j)-Hb_x(3,j))/3
!              Hb_y(1,j)=(4*Hb_y(2,j)-Hb_y(3,j))/3
!
!           end if
!        end do
!
!        return
!        end

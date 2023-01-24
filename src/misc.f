c----------------------------------------------------------------------
c Initialise evolution variables gb and Hb to zero
c----------------------------------------------------------------------
        subroutine init_gbhb(gb_tt,gb_tx,gb_ty,
     &                      gb_xx,gb_xy,gb_yy,
     &                      gb_zz,
     &                      Hb_t,Hb_x,Hb_y,
     &                      L,x,y,chr,ex,Nx,Ny)
        implicit none
        integer Nx,Ny
        real*8 L,x(Nx),y(Ny)
        integer i,j

        real*8 gb_tt(Nx,Ny),gb_tx(Nx,Ny),gb_ty(Nx,Ny)
        real*8 gb_xx(Nx,Ny),gb_xy(Nx,Ny),gb_yy(Nx,Ny)
        real*8 gb_zz(Nx,Ny)
        real*8 Hb_t(Nx,Ny),Hb_x(Nx,Ny),Hb_y(Nx,Ny)
        real*8 chr(Nx,Ny),ex
  
        !--------------------------------------------------------------
 
        do i=1,Nx
            do j=1,Ny
            gb_tt(i,j)=0
            gb_tx(i,j)=0
            gb_ty(i,j)=0
            gb_xx(i,j)=0
            gb_xy(i,j)=0
            gb_yy(i,j)=0
            gb_zz(i,j)=0

            Hb_t(i,j)=0
            Hb_x(i,j)=0
            Hb_y(i,j)=0
           end do
        end do

!        call axi_reg_g(gb_tt,gb_tx,gb_ty,
!     &                 gb_xx,gb_xy,gb_yy,psi,tfunction,chr,ex,
!     &                 L,x,y,z,Nx,Ny,Nz,regtype)

!        call axi_reg_Hb(Hb_t,Hb_x,Hb_y,chr,ex,L,x,y,z,Nx,Ny,Nz,regtype)

        return
        end

c-----------------------------------------------------------------------
c for variables that go to zero linearly at the physical boundary, 
c set residual to zero there
c-----------------------------------------------------------------------
        subroutine lin_zero_bnd_res(f,phys_bdy,all,Nx,Ny,app_dim)
        
        implicit none
        integer Nx,Ny,app_dim,all
        real*8 f(Nx,Ny)
        integer phys_bdy(2*app_dim)

        integer i,j,is,ie,js,je

        ! initialize fixed-size variables
        data i,j,is,ie,js,je/0,0,0,0,0,0/
  
        !--------------------------------------------------------------
 
        if (phys_bdy(1).eq.1.or.all.eq.1) then
            do j=2,Ny-1
                f(2,j)=0
            end do
        end if

        if (phys_bdy(2).eq.1.or.all.eq.1) then
            do j=2,Ny-1
                f(Nx-1,j)=0
            end do
        end if

        if (phys_bdy(3).eq.1.or.all.eq.1) then
            do i=3,Nx-2
                f(i,2)=0
            end do
        end if

        if (phys_bdy(4).eq.1.or.all.eq.1) then
            do i=3,Nx-2
                f(i,Ny-1)=0
            end do
        end if

        return
        end
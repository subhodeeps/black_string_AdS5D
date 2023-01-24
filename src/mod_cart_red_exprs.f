c--------------------------------------------------------------------------------
c Compute the non-vanishing components of the necessary tensorial quantities  
c according to modified Cartoon reduction expressions.
c We reduce the solution to the hypersurface S at omega_1=z=0, y>=0. 
c The other coordinates are time t, and the radial coordinate x.
c We compute the quantities at x(i),y(j) on S.
c
c We use first and second derivatives of the metric as input.
c--------------------------------------------------------------------------------

c--------------------------------------------------------------------------------
c Metric components and their first and second derivatives, 
c as well as Christoffel symbols and their first derivatives 
c--------------------------------------------------------------------------------

        subroutine modcartred_inders_metr(
     &                  g_ll,g_lzlz,
     &					g_ll_x,g_lzlz_x,
     &					g_ll_xx,g_lzlz_xx,
     &					g_uu,g_uzuz,
     &					g_llz_z,
     &					g_uu_x,g_uzuz_x,
     &					g_uuz_z,
     &					gamma_ull,
     &                  gamma_uzllz,gamma_ulzlz,
     &					g_uu_xx,g_uzuz_xx,
     &                  g_llz_xz,g_ll_zz,
     &                  g1_lzlz_zz,g2_lzlz_zz,
     &                  gamma_ull_x,
     &                  gamma_uzllz_x,gamma_ulzlz_x,
     &                  gamma_ullz_z,gamma_uzll_z,
     &                  gamma1_uzlzlz_z,gamma2_uzlzlz_z,
     &                  nsym,x,y,dt,L,Nx,Ny,i,j)


        implicit none

        integer Nx,Ny
        real*8 x(Nx),y(Ny)
        integer i,j
        kron(3,3)
        real*8 x0,y0
        integer a,b,c,d,e,f,g

        integer nsym

        !zeroth derivatives
        real*8 g_ll(3,3),g_lzlz
        real*8 detg
        real*8 g_uu(3,3),g_uzuz
        !first derivatives
        real*8 g_ll_x(3,3,3),g_lzlz_x(3)
        real*8 g_llz_z(3)
        real*8 g_uu(3,3),g_uzuz
        real*8 g_uu_x(3,3,3),g_uzuz_x(3)
        real*8 g_uuz_z(3)
        real*8 gamma_ull(3,3,3)
        real*8 gamma_uzllz(3),gamma_ulzlz(3)
        !second derivatives
        real*8 g_ll_xx(3,3,3,3),g_lzlz_xx(3,3)
        real*8 g_uu_xx(3,3,3,3),g_uzuz_xx(3,3)
        real*8 g_llz_xz(3,3),g_ll_zz(3,3)
        real*8 g1_lzlz_zz,g2_lzlz_zz
        real*8 gamma_ull_x(3,3,3,3)
        real*8 gamma_uzllz_x(3,3),gamma_ulzlz_x(3,3)
        real*8 gamma_ullz_z(3,3),gamma_uzll_z(3,3)
        real*8 gamma1_uzlzlz_z,gamma2_uzlzlz_z

        ! INITIALISE
        !zeroth derivatives
        data detg/0.0/
        data g_uu,g_uzuz/9*0.0,0.0/
        !first derivatives
        data g_ll_x,g_lzlz_x/27*0.0,3*0.0/
        data g_llz_z/3*0.0/
        data g_uu_x,g_uzuz_x/27*0.0,3*0.0/
        data g_uuz_z/3*0.0/
        data gamma_ull/27*0.0/
        data gamma_uzllz,gamma_ulzlz/3*0.0,3*0.0/
        !second derivatives
        data g_ll_xx,g_lzlz_xx/81*0.0,9*0.0/
        data g_uu_xx,g_uzuz_xx/81*0.0,9*0.0/
        data g_llz_xz,g_ll_zz/9*0.0,9*0.0/
        data g1_lzlz_zz,g2_lzlz_zz/0.0,0..0/
        data gamma_ull_x(3,3,3,3)/81*0.0/
        data gamma_uzllz_x,gamma_ulzlz_x/9*0.0,9*0.0/
        data gamma_ullz_z,gamma_uzll_z/9*0.0,9*0.0/
        data gamma1_uzlzlz_z,gamma2_uzlzlz_z/0.0,0.0/


!----------------------------------------------------------------------

        x0=x(i)
        y0=y(j)

		!define Kronecker delta
		do a=1,3
		 do b=1,3
		 	if (a.eq.b) then
		 		kron(a,b)=1
		 	else
		 		kron(a,b)=0
		 	end if
		 end do
		end do

		!zeroth derivatives
		call calc_inv(g_ll,detg,g_uu)
		g_uzuz=1/g_lzlz

		!first derivatives
		do a=1,3
			g_llz_z(a)=(g_ll(a,3)-kron(a,3)*g_lzlz)/y0
		end do


        do a=1,3
          do b=1,3
            do c=1,3
               g_uu_x(a,b,c)=0
               	do d=1,3
               	 do e=1,3
               	  g_uu_x(a,b,c)=g_uu_x(a,b,c)+
     &               (-g_uu(a,d)
     &                   *g_uu(b,e)*g_ll_x(d,e,c))
     			 end do
     			end do
            end do
          end do
        end do

        do a=1,3
        	g_uzuz_x(a)=-(1/g_lzlz**2)*g_lzlz_x(a)
        end do	

        do a=1,3
            g_uuz_z(a)=(g_uu(a,3)-kron(a,3)*g_uzuz)/y0
        end do

        do a=1,3
         do b=1,3
          do c=1,3
            gamma_ull(a,b,c)=0
            do d=1,3
                gamma_ull(a,b,c)=gamma_ull(a,b,c)
     &                          +0.5d0*g_uu(a,d)
     &                                *(g_ll_x(c,d,b)
     &                                 -g_ll_x(b,c,d)
     &                                 +g_ll_x(d,b,c))
          end do
         end do
        end do

        do a=1,3
            gamma_uzllz(a)=0.5d0*g_uzuz*g_lzlz_x(a)
        end do

        do a=1,3
            gamma_ulzlz(a)=0
            do b=1,3
                gamma_ulzlz(a)=gamma_ulzlz(a)
     &              -0.5d0*g_uu(a,b)*g_lzlz_x(b)
     &              +g_uu(a,b)*g_llz_z(b)
            end do
        end do

		!second derivatives
        do a=1,3
         do b=1,3
            g_llz_xz(a,b)=(1/y0)*(g_llz_x(a,b)-kron(a,3)*g_lzlz_x(b))
     &          -kron(b,3)*(1/y0**2)*(g_llz(a)-kron(a,3)*g_lzlz)       
         end do
        end do

        do a=1,3
         do b=1,3
            g_ll_zz(a,b)=(1/y0)*g_ll_z(a,b)
     &          -(1/y0**2)*(kron(a,3)*g_ll(b,3)+kron(b,3)*g_ll(a,3)
     &                  -2*kron(a,3)*kron(b,3)*g_lzlz)       
         end do
        end do

        g1_lzlz_zz=(g_ll(3,3)-g_lzlz)/y0*2

        g2_lzlz_zz=g_lzlz_x(3)/y0

        do a=1,3
         do b=1,3
          do c=1,3
           do d=1,3
            gamma_ull_x(a,b,c,d)=0
            do e=1,3
                gamma_ull_x(a,b,c,d)=gamma_ull_x(a,b,c,d)
     &           +0.5d0*g_uu_x(a,e,d)
     &              *(g_ll_x(b,e,c)-g_ll_x(b,c,e)+g_ll_x(c,e,b))
     &           +0.5d0*g_uu(a,e)
     &              *(g_ll_xx(b,e,c,d)-g_ll_xx(b,c,d,e)+g_ll_xx(c,e,b,d))
            end do   
           end do
          end do
         end do
        end do


        do a=1,3
         do b=1,3
            gamma_uzllz_x(a,b)=0.5d0*g_uzuz_x(b)*g_lzlz(a)
     &          +0.5d0*g_uzuz*g_lzlz_xx(a,b)
         end do
        end do

        do a=1,3
         do b=1,3
            gamma_ulzlz_x(a,b)=0
            do c=1,3
                gamma_ulzlz_x(a,b)=gamma_ulzlz_x(a,b)
     &              +g_uu_x(a,c,b)*(-0.5d0*g_lzlz_x(c)+g_llz_z(c))
     &              +g_uu(a,c)*(-0.5d0*g_lzlz_xx(c,b)+g_llz_xz(c,b))
            end do
         end do
        end do

        do a=1,3
         do b=1,3
            gamma_ullz_z(a,b)=0.5d0*g_uuz_z(a)*g_lzlz_x(b)
            do c=1,3
                gamma_ullz_z(a,b)=gamma_ullz_z(a,b)
     &              +0.5d0*g_uu(a,c)*(g_ll_zz(b,c)-g_llz_xz(b,c)+g_llz_xz(c,b))
            end do
         end do
        end do

        do a=1,3
         do b=1,3
            gamma_uzll_z(a,b)=0.5d0*g_uzuz
     &             *(g_llz_xz(a,b)-g_ll_zz(a,b)+g_llz_xz(b,a))
            do c=1,3
                gamma_uzll_z(a,b)=gamma_uzll_z(a,b)
     &              +0.5d0*g_uuz_z(c)*(g_ll_x(a,c,b)-g_ll_x(a,b,c)+g_ll_x(b,c,a))
            end do
         end do
        end do

        gamma1_uzlzlz_z=0.5d0*g_uzuz*g2_lzlz_zz

        gamma2_uzlzlz_z=g_uzuz*g1_lzlz_zz-0.5d0*g_uzuz*g2_lzlz_zz
        do a=1,3
            gamma2_uzlzlz_z=gamma2_uzlzlz_z
     &          -0.5d0*g_uuz_z(a)*g_lzlz_x(a)+g_uuz_z(a)*g_llz_z(a)       
        end do


        return
        end
!--------------------------------------------------------------------------------


c--------------------------------------------------------------------------------
c Source functions and their first derivatives
c--------------------------------------------------------------------------------


        subroutine modcartred_inders_sourcefuncs(
     &                  g_ll,g_lzlz,
     &                  g_ll_x,g_lzlz_x,
     &                  g_ll_xx,g_lzlz_xx,
     &                  H_u,H_l,
     &                  H_l_x,H_lz_z,
     &                  nsym,x,y,dt,L,Nx,Ny,i,j)


        implicit none

        integer Nx,Ny
        real*8 x(Nx),y(Ny)
        integer i,j
        integer a,b,c,d,e,f,g

        integer nsym

        !zeroth derivatives
        real*8 g_ll(3,3),g_lzlz
        real*8 detg
        real*8 g_uu(3,3),g_uzuz
        !first derivatives
        real*8 g_ll_x(3,3,3),g_lzlz_x(3)
        real*8 g_llz_z(3)
        real*8 g_uu(3,3),g_uzuz
        real*8 g_uu_x(3,3,3),g_uzuz_x(3)
        real*8 g_uuz_z(3)
        real*8 gamma_ull(3,3,3)
        real*8 gamma_uzllz(3),gamma_ulzlz(3)
        real*8 H_u(3),H_l(3)
        !second derivatives
        real*8 g_ll_xx(3,3,3,3),g_lzlz_xx(3,3)
        real*8 g_uu_xx(3,3,3,3),g_uzuz_xx(3,3)
        real*8 g_llz_xz(3,3),g_ll_zz(3,3)
        real*8 g1_lzlz_zz,g2_lzlz_zz
        real*8 gamma_ull_x(3,3,3,3)
        real*8 gamma_uzllz_x(3,3),gamma_ulzlz_x(3,3)
        real*8 gamma_ullz_z(3,3),gamma_uzll_z(3,3)
        real*8 gamma1_uzlzlz_z,gamma2_uzlzlz_z
        real*8 H_l_x(3,3),H_lz_z

        ! INITIALISE
        !zeroth derivatives
        data detg/0.0/
        data g_uu,g_uzuz/9*0.0,0.0/
        !first derivatives
        data g_ll_x,g_lzlz_x/27*0.0,3*0.0/
        data g_llz_z/3*0.0/
        data g_uu_x,g_uzuz_x/27*0.0,3*0.0/
        data g_uuz_z/3*0.0/
        data gamma_ull/27*0.0/
        data gamma_uzllz,gamma_ulzlz/3*0.0,3*0.0/
        data H_u,H_l/3*0.0,3*0.0/
        !second derivatives
        data g_ll_xx,g_lzlz_xx/81*0.0,9*0.0/
        data g_uu_xx,g_uzuz_xx/81*0.0,9*0.0/
        data g_llz_xz,g_ll_zz/9*0.0,9*0.0/
        data g1_lzlz_zz,g2_lzlz_zz/0.0,0..0/
        data gamma_ull_x(3,3,3,3)/81*0.0/
        data gamma_uzllz_x,gamma_ulzlz_x/9*0.0,9*0.0/
        data gamma_ullz_z,gamma_uzll_z/9*0.0,9*0.0/
        data gamma1_uzlzlz_z,gamma2_uzlzlz_z/0.0,0.0/
        data H_l_x,H_lz_z/9*0.0,0.0/


!----------------------------------------------------------------------

        call modcartred_inders_metr(
     &                  g_ll,g_lzlz,
     &                  g_ll_x,g_lzlz_x,
     &                  g_ll_xx,g_lzlz_xx,
     &                  g_uu,g_uzuz,
     &                  g_llz_z,
     &                  g_uu_x,g_uzuz_x,
     &                  g_uuz_z,
     &                  gamma_ull,
     &                  gamma_uzllz,gamma_ulzlz,
     &                  g_uu_xx,g_uzuz_xx,
     &                  g_llz_xz,g_ll_zz,
     &                  g1_lzlz_zz,g2_lzlz_zz,
     &                  gamma_ull_x,
     &                  gamma_uzllz_x,gamma_ulzlz_x,
     &                  gamma_ullz_z,gamma_uzll_z,
     &                  gamma1_uzlzlz_z,gamma2_uzlzlz_z,
     &                  nsym,x,y,dt,L,Nx,Ny,i,j)

        do a=1,3
            H_u(a)=0
            do b=1,3
             do c=1,3
                H_u(a)=H_u(a)
     &              -g_uu(b,c)*gamma_ull(a,b,c)
             end do
            end do
        end do

        do a=1,3
            H_l(a)=0
            do b=1,3
                H_l(a)=H_l(a)
     &              +g_ll(a,b)*H_u(b)
            end do
        end do

        do a=1,3
         do b=1,3
            H_l_x(a,b)=0
            do c=1,3
                H_l_x(a,b)=H_l_x(a,b)
     &              +g_ll_x(a,c,b)*H_u(c)
                do d=1,3
                 do e=1,3
                    H_l_x(a,b)=H_l_x(a,b)
     &                  +g_ll(a,c)*(-g_uu_x(d,e,b)*gamma_ull(c,d,e)
     &                              -g_uu(d,e)*gamma_ull_x(c,d,e,b))
                 end do
                end do
            end do
         end do
        end do

        H_lz_z=-2.0d0*gamma1_uzlzlz_z-nsym*gamma2_uzlzlz_z
        do a=1,3
            H_lz_z=H_lz_z
     &          -g_llz_z(a)*(-H_u(a)+nsym*g_uzuz*gamma_ulzlz(a))
     &          -2.0d0*g_lzlz*g_uuz_z(a)*gamma_uz_llz(a)
            do b=1,3
                H_lz_z=H_lz_z
     &          -g_lzlz*g_uu(a,b)*gamma_uzll_z(a,b)
            end do
        end do

        return
        end
!--------------------------------------------------------------------------------


c--------------------------------------------------------------------------------
c Ricci tensor in generalised harmonic form
c--------------------------------------------------------------------------------



        subroutine modcartred_inders_ricci(
     &                  g_ll,g_lzlz,
     &                  g_ll_x,g_lzlz_x,
     &                  g_ll_xx,g_lzlz_xx,
     &                  ricci_ll,ricci_lzlz,
     &                  nsym,x,y,dt,L,Nx,Ny,i,j)


        implicit none

        integer Nx,Ny
        real*8 x(Nx),y(Ny)
        integer i,j
        integer a,b,c,d,e,f,g

        integer nsym

        !zeroth derivatives
        real*8 g_ll(3,3),g_lzlz
        real*8 detg
        real*8 g_uu(3,3),g_uzuz
        !first derivatives
        real*8 g_ll_x(3,3,3),g_lzlz_x(3)
        real*8 g_llz_z(3)
        real*8 g_uu(3,3),g_uzuz
        real*8 g_uu_x(3,3,3),g_uzuz_x(3)
        real*8 g_uuz_z(3)
        real*8 gamma_ull(3,3,3)
        real*8 gamma_uzllz(3),gamma_ulzlz(3)
        real*8 H_u(3),H_l(3)
        !second derivatives
        real*8 g_ll_xx(3,3,3,3),g_lzlz_xx(3,3)
        real*8 g_uu_xx(3,3,3,3),g_uzuz_xx(3,3)
        real*8 g_llz_xz(3,3),g_ll_zz(3,3)
        real*8 g1_lzlz_zz,g2_lzlz_zz
        real*8 gamma_ull_x(3,3,3,3)
        real*8 gamma_uzllz_x(3,3),gamma_ulzlz_x(3,3)
        real*8 gamma_ullz_z(3,3),gamma_uzll_z(3,3)
        real*8 gamma1_uzlzlz_z,gamma2_uzlzlz_z
        real*8 H_l_x(3,3),H_lz_z
        real*8 ricci_ll(3,3),ricci_lzlz

        ! INITIALISE
        !zeroth derivatives
        data detg/0.0/
        data g_uu,g_uzuz/9*0.0,0.0/
        !first derivatives
        data g_ll_x,g_lzlz_x/27*0.0,3*0.0/
        data g_llz_z/3*0.0/
        data g_uu_x,g_uzuz_x/27*0.0,3*0.0/
        data g_uuz_z/3*0.0/
        data gamma_ull/27*0.0/
        data gamma_uzllz,gamma_ulzlz/3*0.0,3*0.0/
        data H_u,H_l/3*0.0,3*0.0/
        !second derivatives
        data g_ll_xx,g_lzlz_xx/81*0.0,9*0.0/
        data g_uu_xx,g_uzuz_xx/81*0.0,9*0.0/
        data g_llz_xz,g_ll_zz/9*0.0,9*0.0/
        data g1_lzlz_zz,g2_lzlz_zz/0.0,0..0/
        data gamma_ull_x(3,3,3,3)/81*0.0/
        data gamma_uzllz_x,gamma_ulzlz_x/9*0.0,9*0.0/
        data gamma_ullz_z,gamma_uzll_z/9*0.0,9*0.0/
        data gamma1_uzlzlz_z,gamma2_uzlzlz_z/0.0,0.0/
        data H_l_x,H_lz_z/9*0.0,0.0/
        data ricci_ll,ricci_lzlz/9*0.0,0.0/


!----------------------------------------------------------------------


        call modcartred_inders_metr(
     &                  g_ll,g_lzlz,
     &                  g_ll_x,g_lzlz_x,
     &                  g_ll_xx,g_lzlz_xx,
     &                  g_uu,g_uzuz,
     &                  g_llz_z,
     &                  g_uu_x,g_uzuz_x,
     &                  g_uuz_z,
     &                  gamma_ull,
     &                  gamma_uzllz,gamma_ulzlz,
     &                  g_uu_xx,g_uzuz_xx,
     &                  g_llz_xz,g_ll_zz,
     &                  g1_lzlz_zz,g2_lzlz_zz,
     &                  gamma_ull_x,
     &                  gamma_uzllz_x,gamma_ulzlz_x,
     &                  gamma_ullz_z,gamma_uzll_z,
     &                  gamma1_uzlzlz_z,gamma2_uzlzlz_z,
     &                  nsym,x,y,dt,L,Nx,Ny,i,j)

        call modcartred_inders_sourcefuncs(
     &                  g_ll,g_lzlz,
     &                  g_ll_x,g_lzlz_x,
     &                  g_ll_xx,g_lzlz_xx,
     &                  H_u,H_l,
     &                  H_l_x,H_lz_z,
     &                  nsym,x,y,dt,L,Nx,Ny,i,j)

        do a=1,3
         do b=1,3
            ricci_ll(a,b)=-0.5d0*(H_l_x(a,b)+H_l_x(b,a))
            do c=1,3
                ricci_ll(a,b)=ricci_ll(a,b)
     &              +H_l(c)*gamma_ull(c,a,b)
                do d=1,3
                    ricci_ll(a,b)=ricci_ll(a,b)
     &                  -0.5d0*g_uu(c,d)*g_ll_xx(a,b,c,d)
     &                  -0.5d0*(g_uu_x(c,d,a)*g_ll_x(b,c,d)
     &                      +g_uu_x(c,d,b)*g_ll_x(a,c,d))
     &                  -gamma_ull(c,d,a)*gamma_ull(d,c,b)
                end do
            end do
         end do
        end do

        ricci_lzlz=-g_uzuz*g1_lzlz_zz-(nsym/2)*g_uzuz*g2_lzlz_zz
     &              -H_lz_z
        do a=1,3
            ricci_lzlz=ricci_lzlz
     &          -g_uuz_z(a)*g_llz_z(a)-g_uuz_z(a)*g_lzlz_x(a)
     &      +(H_l(a)+(nsym/2)*g_uzuz*g_lzlz_x(a)
     &          -nsym*g_uzuz*g_llz_z(a))*gamma_ulzlz(a)
     &      -2.0d0*gamma_ulzlz(a)*gamma_uzllz(a)
            do b=1,3
                ricci_lzlz=ricci_lzlz
     &              -0.5d0*g_uu(a,b)*g_lzlz_xx(a,b)
            end do
        end do


        return
        end
!--------------------------------------------------------------------------------

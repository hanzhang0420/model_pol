      SUBROUTINE cal_integrate(ARGC, ARGV) 
C     external function that is called in pm_integrate.pro 
      implicit none
      INTEGER*8 ARGC
      INTEGER*8 ARGV(*)                   !Argc and argv are integers copied from the example on website
C      j = LOC(argc)
           
      call integrate( %VAL(ARGV(1)),%VAL(ARGV(2)),%VAL(ARGV(3)),
     &    %VAL(ARGV(4)),%VAL(ARGV(5)),%VAL(ARGV(6)),%VAL(ARGV(7)),
     &     %VAL(ARGV(8)),%VAL(ARGV(9)),%VAL(ARGV(10)),%VAL(ARGV(11)),
     &     %VAL(ARGV(12)),%VAL(ARGV(13)))

      return
      end 
     

      subroutine integrate(kp,kl,wave,dens,temp,R,
     &     I_i,Q_i,U_i,ncube,xgrid,gam,xi)
      IMPLICIT NONE
c     return I_i,Q_i,U_i
     
      real*8,parameter :: pi = 4*atan(1.0D0)
      real*8,parameter :: deg2rad = pi / 180.d0
      
      integer*4 :: nx,ny,nz,ncube
c     parameter (px=202,py=202,pz=202)
      
      REAL*8,dimension(1:ncube,1:ncube,1:ncube) :: dens,temp
      REAL*8,dimension(1:ncube,1:ncube) :: I_i,Q_i,U_i
      REAL*8:: xgrid(ncube+1)
      REAL*8:: del_z,dtau,dtaup,zeta,I_p_i,beta
      REAL*8:: wave,R,kp,kl,bbflux,I,U,Q,em,incl
      REAL*8,dimension(1:ncube,1:ncube,1:ncube)::gam,xi
c     REAL*8,dimension(1:ncube,1:ncube,1:ncube)::bx,by,bz
      REAL*8:: cos_beta,sin_beta

      print*,'dimension',ncube    
      print*,'external fortran code'
      
      do 111 nx=1,ncube            
         do 112 ny=1,ncube 
            
             I=0.0d0 
             U=0.0d0 
             Q=0.0d0
             do 113 nz=1,ncube
                
                del_z=(xgrid(nz+1)-xgrid(nz))*1.49598d13 !not linear spaced
  
                dtau=((kp+kl)*dens(nx,ny,nz)*del_z)/2.0d0
       
                dtaup=((kp-kl)*dens(nx,ny,nz)*del_z
     &*((COS(gam(nx,ny,nz)))**2))*R/2.0d0

c print*,(COS(gam(nx,ny,nz))),gam(nx,ny,nz)
                
                beta=2.0d0*xi(nx,ny,nz)
                zeta=-Q*COS(beta)-U*SIN(beta)
                
                I_p_i=I*dsinh(dtaup)+zeta*dcosh(dtaup)-zeta
       
       call planck(wave,temp(nx,ny,nz),bbflux)

       I=(I*dcosh(dtaup)+zeta*dsinh(dtaup))*exp(-dtau)
     &+(1.0-dcosh(dtaup)*exp(-dtau))*bbflux

   
     
       em=dsinh(dtaup)*exp(-dtau)*bbflux
      
c this part because of the precision issue  
      cos_beta=cos(beta)
      sin_beta=sin(beta)
c      if (abs(cos_beta).lt.1.0d-14) then
c      cos_beta=0.d0
c      end if 

c      if (abs(sin_beta).lt.1.0d-14) then
c      sin_beta=0.d0
c      end if
      
c      print*,beta,sin_beta,cos_beta
      
       Q=(Q-I_p_i*cos_beta)*exp(-dtau)+em*cos_beta
       
       U=(U-I_p_i*sin_beta)*exp(-dtau)+em*sin_beta
                       
             
 113  end do
c      print*,'done'
c      print*,I,U,Q  
      I_i(nx,ny)=I
      U_i(nx,ny)=U
      Q_i(nx,ny)=Q
  
    
 112   end do
 111   end do
       return
       end


      subroutine planck(wave,temp,bbflux)
      
      double precision bbflux,wave,w,c,miu,val
      double precision c1,c2,temp
       w = wave/1.d8                             
c     ; Angstroms to cm
       c = 2.99792458d10
       miu = c/w
       c1 =  4.6286132d-47                
c             =2*!DPI*h/c/c
       c2 =  4.7959d-11                  
c              ; =h/k
       val =  c2*miu/temp

      bbflux =  c1 * miu**3.0d0 /(exp(val)-1.0d0)

      return
      end


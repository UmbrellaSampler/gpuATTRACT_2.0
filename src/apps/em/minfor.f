      subroutine minfor(
     1 smug, ivmax, x, numModesRec, numModesLig
     2 )
c
c  variable metric minimizer (Harwell subroutine lib.  as in Jumna with modifications)
c     minimizes a single structure


      use, intrinsic :: iso_c_binding, only: c_ptr
      implicit none

c     Parameters
      type(c_ptr) smug
      integer ivmax
      real*8 x
      integer numModesRec, numModesLig
      dimension x(26)

c     Local variables
      real*8  gesa
      integer i,k,ir,isfv,itr,nfun,np,jn
      real*8 c,acc,dff,dgb,fa,fmin,gl1,gl2,gmin,dga,xnull,w
      real*8 fb, delta
      
      real*8 h,g,ga,gb,xaa,xbb,d,step,stepbd,steplb,stmin
      dimension h(26*26)
      dimension g(26),ga(26),gb(26),w(26), delta(26)
      dimension xaa(26), xbb(26), d(26)

      integer, parameter:: ERROR_UNIT = 0

      xnull=0.0d0
      dga=xnull
      numModesRec = 1
       numModesLig = 1
      jn = 6 + numModesRec + numModesLig

      do i=1,jn
       ga(i)=xnull
       d(i)=xnull
       xaa(i)=x(i)
c      write(ERROR_UNIT,*)'xaa',i,xaa(i)
      enddo



      nfun=0
      itr=0
      np=jn+1
      acc=0.000000001D0
c     set the hessian to a diagonal matrix 
      c = 1.0D0
      k=(jn*np)/2
      do i=1,k
       h(i)=xnull
      enddo
      k=1
      do i=1,jn
       h(k)=0.01d0*c
       k=k+np-i
      enddo

c     set some variables for the first iteration
      dff=xnull
c changed g to delta
c      write(*,*), "first iter "
cc      write(*,*), xaa(1)
c      write(*,*), xaa(2)
c      write(*,*), xaa(3)
c      write(*,*), xaa(7)
      call energy_for_fortran_to_call(smug, xaa, gesa, delta)
c  g to delta
      do i=2,jn-2
c       delta(i)=0
      enddo
c      delta(4)=0
c      delta(5)=0
c      delta(6)=0
      delta(8)=0
c     write(*,*),"deltas"
c      write(*,*), delta(1)
c      write(*,*), delta(2)
c      write(*,*), delta(3)
c      write(*,*), delta(7)


c      write(*,*), delta(7)
c      write(*,*), "iter: ", itr
c     write(*,*),"Gradients"
c      write(*,*),energies(1)
c      write(*,*),"rota        ", xaa(1), xaa(2), xaa(3)
c      write(*,*),"translation ", xaa(4), xaa(5), xaa(6)
c            write(*,*), "modeR", xbb(7)
c      write(*,*), "modeL", xbb(8)
c      write(*,*),"#################################forst call "
      
110   fa=gesa
      isfv=1
c store forces, Euler angle, position and ligand and receptor coordinates
      do i=1,jn
c       write(*,*), "delta ", i , " ", g(i)
c changed added line above ga(i) = g(i)
       g(i)=-delta(i)
       ga(i)=g(i)
      enddo
      
      
c     begin the iteration by giving the required printing
135   itr=itr+1
c     calculate the search direction of the iteration
      do i=1,jn
       d(i)=-ga(i)

      enddo
      call mc11e (h,jn,d,w,jn)
      do i=1,jn
c      write(*,*)'d(i)',i,d(i)
      enddo

c     calculate a lower bound on the step-length
c     and the initial directional derivative
      c=xnull
      dga=xnull
      do i=1,jn
      c=max(c,abs(d(i)))
      dga=dga+ga(i)*d(i)
      enddo
c     test if the search direction is downhill
      if (dga.ge.xnull) go to 240
c     set the initial step-length of the line search
      stmin=xnull
      stepbd=xnull
      steplb=acc/c
      fmin=fa
      gmin=dga
      step=1.0d0
      if (dff.le.xnull) step=min(step,1.0d0/c)
      if (dff.gt.xnull) step=min(step,(dff+dff)/(-dga))
170   c=stmin+step
c     test whether func has been called ivmax times
      if (nfun.ge.ivmax) go to 250
      nfun=nfun+1
c     calculate another function value and gradient
c changed added the do loop below
      do i=1,jn
      g(i)=-delta(i)
c      write(ERROR_UNIT,*)'g(i)',g(i),d(i)
      enddo
c     make an Euler rotation + tranlation of ligand center
c      write(*,*), "call next energy ", itr
      do i=1,jn-2
c changed + to minus
       xbb(i)=xaa(i)+c*d(i)
c       write(*,*), xbb(i)
      enddo

      do i=jn-2,jn
c changed + to minus
       xbb(i)=xaa(i)-c*d(i)
c       write(*,*), xbb(i)
      enddo

c      write(*,*), xbb(1)
c      write(*,*), xbb(2)
c      write(*,*), xbb(3)
c      write(*,*), xbb(7)
c changed to delta from gb

      call energy_for_fortran_to_call(smug, xbb, fb, delta)
c      write(*,*), "iter: ", itr

c     store this function value if it is the smallest so far
c      delta(4)=0
c      delta(5)=0
c      delta(6)=0
      delta(8)=0
c      write(*,*),"deltas"
c      write(*,*), delta(1)
c      write(*,*), delta(2)
c      write(*,*), delta(3)
c      write(*,*), delta(7)

      do i=1, jn
        gb(i)=-delta(i)
      enddo
      do i=1,jn-2
c       gb(i)=0
      enddo
c      gb(1)=-delta(1)
      gb(8)=0
      do i=1,jn
c       write(*,*), "E i", i, " ", gb(i)
      enddo
c      write(*,*), "E i", "7", " ", delta(7)
      isfv=min(2,isfv)
      if (fb.gt.gesa) go to 220
      if (fb.lt.gesa) go to 200
      gl1=xnull
      gl2=xnull
      do i=1,jn
       gl1=gl1+(g(i))**2
       gl2=gl2+(gb(i))**2
      enddo
      if (gl2.ge.gl1) go to 220
200   isfv=3
      gesa=fb
c CHANGED

c      write(*,*), "iter: ", itr
c      write(*,*),"Gradients"
c      write(*,*),energies(1)
c      write(*,*),"rota        ", xbb(1), xbb(2),xbb(3)
c      write(*,*),"translation ", xbb(4), xbb(5), xbb(6)
c      write(*,*), "modeR", xbb(7)
c      write(*,*), "modeL", xbb(8)
c      write(*,*)," "
      do i=1,jn
       x(i) = xbb(i)
       g(i)=gb(i)
c       write(*,*), "delta ", i , " ", gb(i)
      enddo
      

220   dgb=xnull
      do i=1,jn
       dgb=dgb+gb(i)*d(i)
      enddo
      
c     branch if we have found a new lower bound on the step-length
      if (fb-fa.le.0.1d0*c*dga) go to 280
c     finish the iteration if the current step is steplb
      if (step.gt.steplb) go to 270
240   if (isfv.ge.2) go to 110

c     at this stage the whole calculation is complete
250   if(nfun.lt.ivmax) then
      nfun=nfun+1
      endif
      return
c       LKHOIEHOHFIUGIGOIZGI
c#######################################################################################
c     calculate a new step-length by cubic interpolation
270   stepbd=step
      c=gmin+dgb-3.0d0*(fb-fmin)/step
      c=gmin/(c+gmin-sqrt(c*c-gmin*dgb))
      step=step*max(0.1d0,c)
      go to 170
c     set the new bounds on the step-length
280   stepbd=stepbd-step
      stmin=c
      fmin=fb
      gmin=dgb
c     calculate a new step-length by extrapolation
      step=9.0d0*stmin
      if (stepbd.gt.xnull) step=0.5d0*stepbd
      c=dga+3.0d0*dgb-4.0d0*(fb-fa)/stmin
      if (c.gt.xnull) step=min(step,stmin*max(1.0d0,-dgb/c))
      if (dgb.lt.0.7d0*dga) go to 170
c     test for convergence of the iterations
      isfv=4-isfv
      if (stmin+step.le.steplb) go to 240
c     revise the second derivative matrix
      ir=-jn
      do i=1,jn
       xaa(i)=xbb(i)
       xbb(i)=ga(i)
       d(i)=gb(i)-ga(i)
       ga(i)=gb(i)
      enddo
      call mc11a(h,jn,xbb,1.0d0/dga,w,ir,1,xnull)
      ir=-ir
      call mc11a (h,jn,d,1.0d0/(stmin*(dgb-dga)),d,ir,0,xnull)
c     branch if the rank of the new matrix is deficient
      if (ir.lt.jn) go to 250
c     begin another iteration
      dff=fa-fb
      fa=fb
      go to 135
      end


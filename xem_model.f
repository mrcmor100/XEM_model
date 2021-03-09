      program run_xem_model
      implicit none

      integer eof
      character*80 rawname, filename, ofilename
      real*8 E0, EP, THETA, A, Z, dis_XS, qe_XS
      integer*4	last_char

      read(*,1968) rawname
 1968 format(a)
CAM Open input file.
      filename = 'input/'//rawname(1:last_char(rawname))//'.inp'
      open(unit=11,status='old',file=filename)

CAM Open output file.
      filename = './output/'//rawname(1:last_char(rawname))//'.out'
      open (unit=14,status='unknown',file=filename)

 40   READ(11,*,IOSTAT=EOF) E0, EP, THETA, A, Z
      if(eof.ge.0) then
        call xem_model(E0, EP, THETA, A, Z, dis_XS, qe_XS)
        goto 40
      endif !if file has stuff
      close(11);
      end

!-----------------------------------------------------------------------------
CAM Main subroutine used in externals.
!-----------------------------------------------------------------------------
      subroutine xem_model(E0, EP, THETA, A, Z, sigdis_new,sig_qe_new)
      implicit none

      COMMON /Y_SCALING/ ag, bg, bigB, f0, alpha1, pfermi, eps
      real*8 ag, bg, bigB, f0, alpha1, pfermi, eps

      include 'constants_dble.inc'
      real*8 E0, EP, THETA, QSQ, NU, THR
      real*8 Y, A, Z, X, N
      real*8 INNP, INNT
      real*8 sig_qe_new, sigdis_new
      real*8 YSCALE
!      include 'model_constants.inc'
      logical first/.true./
      if(first) then
         first=.false.
         call load_parameters(A, Z)
      endif
      INNT = 30
      INNP = 30

      N = A - Z
      THR = THETA*pi/180.
      QSQ = 4*E0*EP*(SIN(THR/2))**2
      NU  = E0 - EP
      X   = QSQ/(2*nuc_mass*NU)
      Y   = YSCALE(E0,EP,THR,A,EPS)
      !Check if Y is right kinda..
      if(Y.lt.1.E19) then 
         call sig_dis_calc(int(A), int(Z), E0, EP, THR, Y, sigdis_new)
         call sig_qe_calc(Y, int(A), int(Z),E0,EP, THR, X, sig_qe_new)
      else
         sigdis_new = 0.0
         sig_qe_new = 0.0
      endif
c      write(236,*) A, Z, X, Y

      return
      end

!-----------------------------------------------------------------------------
CAM Calculate DIS using donal's smearing routine.
!-----------------------------------------------------------------------------
      subroutine sig_dis_calc(a, z, e, ep, thr, y, sigdis)
      implicit none

      COMMON /Y_SCALING/ ag, bg, bigB, f0, alpha1, pfermi, eps
      real*8 ag, bg, bigB, f0, alpha1, pfermi, eps
      COMMON /CORRECTIONS/ johns_fac, fact, corfact, dhxcorfac, emc_corr
      real*8 johns_fac, fact, corfact, dhxcorfac, emc_corr

      include 'constants_dble.inc'
      integer a, z, n, innp, innt
      real*8 e, ep, thr, wsq, x,y, qsq, nu
      real*8 sigdis, pmax
      real*8 f1, f2, w1, w2, tan_2, sig_mott
      real*8 emc_func_xem
!      include 'model_constants.inc'

      n=a-z
      qsq=4.0*e*ep*sin(thr/2)**2
      nu=e-ep
      x=qsq/2.0/nuc_mass/nu
<<<<<<< f7e0ffcda06cc96c2eb90ee2c9c812056d04ae12
      WSQ = -qsq + nuc_mass**2 + 2.0*nuc_mass*nu
=======

      x1=0.8                    !x<x1 --> use emc corrected ld2
      x2=0.9                    !x1<=x<x2 --> smooth transition from emc corrected ld2 to donals smearing

      n=a-z
c      write (*,*) 'New DiS SUB'
      call sig_bar_df(e, ep, thr*180.0/pi, y,dble(0.0), sigma_p,
     +  sigma_n)
>>>>>>> Add most updated xem_model fortran code.  Abstract reading of E0, EP, THETA, A, Z into xem_model.

      PMAX=1.0

      if(a.ge.2) then
         if(WSQ.lt.2.25) then
            innt=30
            innp=30
         else
            innt=30
            innp=10
         endif

<<<<<<< f7e0ffcda06cc96c2eb90ee2c9c812056d04ae12
         CALL smear4all(E,EP,THR,dble(A),dble(Z),dble(N),eps,PMAX,
     +        dble(INNP),dble(INNT),f0,bigB,ag,bg,
     +        alpha1 ,SIGDIS)
=======
c         write(*,*) 'New Sub, passing variables', E,EP,THR,A,Z,N,EPS(a),
c     +        PMAX,INNP,INNT,f0(a),bigB(a),ag(a),bg(a),alpha1(a)
         CALL BDISNEW4HE3(E,EP,THR,dble(A),dble(Z),dble(N),EPS(a),PMAX,
     +        dble(INNP),dble(INNT),f0(a),bigB(a),ag(a),bg(a),
     +        alpha1(a) ,SIGDIS)

         if (x.lt.x1) then
            emc_corr = emc_func_xem(x,a)
         elseif ((x.ge.x1).and.(x.lt.x2)) then
            frac = (x-x1)/(x2-x1)
            emc_corr = 1.0*frac + emc_func_xem(x,a)*(1.-frac)
         elseif(x.ge.x2) then
            emc_corr = 1.0
         endif
>>>>>>> Add most updated xem_model fortran code.  Abstract reading of E0, EP, THETA, A, Z into xem_model.

CAM EMC_FUNC specific correction.
CAM D2 set to F1F221.
         emc_corr = emc_func_xem(x,a)
         sigdis = sigdis*emc_corr

CAM Make a correction to the high_x tail.
         dhxcorfac=1.
         if ((x.gt.0.9)) then
            call  dis_highx_cor(a,x,dhxcorfac)
            sigdis = sigdis*dhxcorfac
         endif

      else if(A.eq.1) then
       call F1F2IN21(dble(Z),dble(A), QSQ, WSQ, F1, F2)
       W1 = F1/.93827231D0
       W2 = F2/nu
       tan_2 = tan(thr/2)**2
CAM    Mott cross section
       sig_mott = 0.3893857*(alpha*cos(thr/2))**2/(2*e*sin(thr/2)**2)**2
       sigdis =  sig_mott*(W2+2.0*W1*tan_2)*1.e-2*1.e9
      endif
      end

!-----------------------------------------------------------------------------
      subroutine sig_qe_calc(y, a, z, e, ep, thr, x, sig_qe)
      implicit none
      include 'constants_dble.inc'
<<<<<<< f7e0ffcda06cc96c2eb90ee2c9c812056d04ae12

      COMMON /Y_SCALING/ ag, bg, bigB, f0, alpha1, pfermi, eps
      real*8 ag, bg, bigB, f0, alpha1, pfermi, eps
      COMMON /JOHNS_FAC/ j_fac
      real*8 j_fac
      COMMON /CORRECTIONS/ johns_fac, fact, corfact, dhxcorfac, emc_corr
      real*8 johns_fac, fact, corfact, dhxcorfac, emc_corr

      integer a, z
      real*8 thr, y, e, ep, x, nu, q4_sq, q3v
      real*8 sig_qe, sigma_n, sigma_p
      real*8 dwdy
      real*8 tail_cor
!      include 'model_constants.inc'

      q4_sq = 4.*e*ep*sin(thr/2.)**2 !4-mom. transf. squared.
      nu = e-ep                 !Energy loss
      q3v = sqrt(q4_sq + nu**2) !3-momentum transfer
      dwdy = q3v/sqrt(nuc_mass**2+q3v**2+y**2+2.0*q3v*y)
      x=q4_sq/2.0/nuc_mass/nu

CAM F(y) from Nadia's Thesis based on:
CAM https://arxiv.org/pdf/nucl-th/9702009.pdf
      y=y*1000
      if(a.eq.2.) then
        sig_qe = (f0-bigB)*alpha1**2*exp(-(ag*y)**2)
     +    /(alpha1**2+y**2)+bigB*exp(-bg*abs(y))
      else
        sig_qe = (f0-bigB)*alpha1**2*exp(-(ag*y)**2)
     +    /(alpha1**2+y**2)+bigB*exp(-(bg*y)**2)
      endif
      y=y/1000.
CAM JOHNS_FAC basically linear factor less than X=1
      if(j_fac.ne.0.0) then
         johns_fac=max(1.,1.+y*1.4*j_fac)
      else
         johns_fac=1.
      endif
      sig_qe=johns_fac*sig_qe

CAM Make get the off-shell contributions based on DeForest:
CAM T. De Forest, Jr. Nuc. Phys. A392 232 (1983)
      call sig_bar_df(e, ep, thr*180/pi, y,dble(0.0), sigma_p, sigma_n)
      fact = (Z*sigma_p+(A-Z)*sigma_n)/dwdy
      sig_qe=sig_qe*fact*1000.

CAM tail_cor applied here
      corfact=tail_cor(x)
      sig_qe=sig_qe*corfact

      return
=======
      integer a,z
      real*8 input, sig_qe, johns_fac
      real*8 q4_sq, q3v, dwdy,thr, fact,nu
      real*8 x, x_cor, corfact, my_frac, x_high, x_low
      real*8 sigma_n, sigma_p, theta
      real*8 y, e, ep, tail_cor
      include 'model_constants.inc'

        y=y*1000
!        write(*,*) 'my params', f0(a), bigB(a), alpha1(a), ag(a), bg(a)
        if(a.eq.2.) then
          sig_qe = (f0(a)-bigB(a))*alpha1(a)**2*exp(-(ag(a)*y)**2)
     +      /(alpha1(a)**2+y**2)+bigB(a)*exp(-bg(a)*abs(y))
        else
          sig_qe = (f0(a)-bigB(a))*alpha1(a)**2*exp(-(ag(a)*y)**2)
     +      /(alpha1(a)**2+y**2)+bigB(a)*exp(-(bg(a)*y)**2)
!       write(*,*) f0(a), ag(a), bg(a),bigB(a), alpha1(a)
       endif
        y=y/1000
        johns_fac=1.
        if (a.eq.3) then
          johns_fac=max(1.,1.+y*1.4*2.5)
        elseif (a.eq.4) then
          johns_fac=max(1.,1.+y*1.4*1.75)
        elseif (a.eq.12) then
          johns_fac=max(1.,1.+y*1.4*2.5)
        elseif ((a.eq.56).or.(a.eq.64)) then
          johns_fac=max(1.,1.+y*1.4*3)
        elseif (a.eq.197) then
          johns_fac=max(1.,1.+y*1.4*4)
        endif

       sig_qe=johns_fac*sig_qe

  	q4_sq = 4.*e*ep*sin(thr/2.)**2 		!4-mom. transf. squared.
	nu = e-ep				!Energy loss
	q3v = sqrt(q4_sq + nu**2)		!3-momentum transfer
	dwdy = q3v/sqrt(nuc_mass**2+q3v**2+y**2+2.0*q3v*y)

        call sig_bar_df(e, ep, thr*180/pi, y,dble(0.0), sigma_p, sigma_n)

        fact = (Z*sigma_p+(A-Z)*sigma_n)/dwdy
        sig_qe=sig_qe*fact*1000.

	if(x.gt.2.0) then
	  x_cor=2.0
        else
          x_cor=x
	endif

	if(a.eq.2) then
	  x_low=1
	  x_high=1.05
	endif

	if(a.gt.2) then
	  if((a.eq.64).or.(a.eq.4).or.(a.eq.197).or.(a.eq.56)) then
	    x_low=1.2
	    x_high=1.4
	  else
	    x_low=1.4
	    x_high=1.6
	  endif
	endif

        corfact=1.
	if((x.ge.(x_low))) then !.and.(A.eq.2)) then
          corfact=tail_cor(x_cor,a)
          if(x.lt.(x_high)) then
	    my_frac=(x-x_low)/(x_high-x_low)
            corfact=my_frac*corfact+(1.-my_frac)
	  endif
          sig_qe=sig_qe*corfact
        endif
>>>>>>> Add most updated xem_model fortran code.  Abstract reading of E0, EP, THETA, A, Z into xem_model.
      end

!---------------------------------------------------------------------
CAM Function to determine Y.  Refer to John A's thesis
!---------------------------------------------------------------------
      real*8 FUNCTION YSCALE(E,EP,THR,A,EPS)
      IMPLICIT NONE
      INCLUDE 'constants_dble.inc'
      real*8 E, EP,THR,A, NU,W,WP,AG,BG,PART1,PART2,PART3,EPS,QSQ,CG
      real*8 backup,RAD

      yscale = 0.0
      NU = E-EP
      QSQ = 4*E*EP*(SIN(THR/2))**2

<<<<<<< f7e0ffcda06cc96c2eb90ee2c9c812056d04ae12
=======
      qv2=qsq+nu**2
>>>>>>> Add most updated xem_model fortran code.  Abstract reading of E0, EP, THETA, A, Z into xem_model.
      W = NU+A*NUC_MASS-EPS
      WP = W**2+((A-1)*NUC_MASS)**2-NUC_MASS**2
                                !write (15,*) 'w, wp, qsq', w,wp,qsq
      AG = 4*W**2-4*(NU**2+QSQ)
      BG = SQRT(NU**2+QSQ)*(4*WP-4*(NU**2+QSQ))
      PART1 = 4*W**2*(((A-1)*NUC_MASS)**2)
      PART2 = 2*WP*(QSQ+NU**2)
      PART3 = (QSQ+NU**2)**2
      CG = PART1+PART2-PART3-WP**2
      rad = BG**2-4*AG*CG
                                !write (15,*) 'A, B, C, RAD', AG,BG,CG,RAD
      if (rad.ge.0) then
         YSCALE = (-BG+SQRT(BG**2-4*AG*CG))/(2*AG)
         backup = (-BG-SQRT(BG**2-4*AG*CG))/(2*AG)
      else
         YSCALE = 1.E20
      endif

<<<<<<< f7e0ffcda06cc96c2eb90ee2c9c812056d04ae12
=======
      RETURN
      END


      real*8 FUNCTION CS_MOTT(THETA,ESTART)
      IMPLICIT NONE
      INCLUDE 'constants_dble.inc'
      real*8 THETA,ESTART,hbarc2

      hbarc2=.38937966
                                !ALPHA = 1/137.02
      CS_MOTT = 0.1973289**2*1e10*(ALPHA**2)*(COS(THETA/2)**2)/(4*(ESTART**2)*(SIN(THETA/2))**4)
 !     CS_MOTT = hbarc2*(ALPHA**2)*(cos(theta/2))**2/(4*(ESTART**2)*(SIN(THETA/2))**4)
                                !(.197**2*1E10)*
                                !(3.893869477*1E9)*
                                !write(6,*)'in cs_mott, got ',cs_mott
>>>>>>> Add most updated xem_model fortran code.  Abstract reading of E0, EP, THETA, A, Z into xem_model.
      RETURN
      END

!-----------------------------------------------------------------------------------------------------

      subroutine dis_highx_cor(anuc,x,cor)
      implicit none
      
      COMMON /DHX_COR/ dhx_xlow,dhx_xhigh,dhx_cor_min,dhx_cor_xalt,
     >     dhx_cor_alt_val, dhx_cor_0, dhx_cor_1
      real*8 dhx_xlow,dhx_xhigh,dhx_cor_min,dhx_cor_xalt,
     >     dhx_cor_alt_val, dhx_cor_0, dhx_cor_1
      
      integer anuc
      real*8 x,cor, frac

      frac=1
      if (x.gt.dhx_xlow) then
         cor=dhx_cor_1*x+ dhx_cor_0
         if((x.ge.dhx_xlow).and.(x.le.dhx_xhigh)) then
            frac = (x-dhx_xlow)/(dhx_xhigh-dhx_xlow)
         endif
      endif
      if (x.gt.dhx_cor_xalt.and.dhx_cor_xalt.ne.0.0) then
         cor=dhx_cor_alt_val
      endif
      cor=frac*cor+1.-frac

      if(cor.lt.0.4) cor=0.4
         
      return
      end

c--------------------------------------------------------------------------------
CAM tail_cor used to fix tail of QE cross-section.
c--------------------------------------------------------------------------------
      real*8 function tail_cor(x)
      implicit none

      COMMON /QE_TAIL_COR/ tc_xlow, tc_xhigh, tc_aa, tc_bb, tc_cc, 
     >       tc_dd, tc_ee, tc_ff, tc_const
      real*8 tc_xlow, tc_xhigh, tc_aa, tc_bb, tc_cc, 
     >       tc_dd, tc_ee, tc_ff, tc_const
      real*8  x, x_cor, corfact, my_frac!, aa, bb, cc, dd, ee, ff

      if(x.gt.tc_const) then
         x_cor = tc_const
      else
         x_cor=x
      endif
      corfact=1
      if((x.ge.(tc_xlow))) then
         corfact=(tc_aa*exp(tc_bb*x_cor) + 
     1   tc_cc*x_cor**6+tc_dd*x_cor**4+tc_ee*x_cor**2+tc_ff)
         if(x.lt.(tc_xhigh)) then
	    my_frac=(x-tc_xlow)/(tc_xhigh-tc_xlow)
            corfact=my_frac*corfact+(1.-my_frac)
         endif
      endif

      tail_cor=corfact
      return
      end

!this is the postthesis iteration of daves
      real*8 function emc_func_xem(x,A) ! now compute the emc effect from our own fits.
      implicit none

      COMMON /EMC_COR/ emc_xlow, emc_xhigh, emc_0, emc_1, emc_2, emc_3,
     >       emc_4, emc_5
      real*8 emc_xlow, emc_xhigh, emc_0, emc_1, emc_2, emc_3,
     >       emc_4, emc_5

      integer a
      real*8 x, xtmp, frac
      real*8 emc, emc_corr
      xtmp = x
      emc = emc_0 + emc_1*xtmp + emc_2*xtmp**2 +
     1     emc_3*xtmp**3 + emc_4*xtmp**4 + emc_5*xtmp**5
c      write(6,*) emc

c      emc = 1.42308 - 2.65087*xtmp + 11.41047*xtmp**2 -
c     1     22.54747*xtmp**3 + 18.47078*xtmp**4 - 5.07537*xtmp**5
c      write(6,*) emc
        
!FROM SIG_DIS_CALC
!CAM EMC_FUNC specific correction.
!CAM D2 set to F1F221.
      if (x.lt.emc_xlow) then
         emc_corr = emc
      elseif ((x.ge.emc_xlow).and.(x.lt.emc_xhigh)) then
         frac = (x-emc_xlow)/(emc_xhigh-emc_xlow)
         emc_corr = 1.0*frac + emc*(1.-frac)
CAM Uncomment when setting emc_func_xem
CAM           emc_corr = emc
      elseif(x.ge.emc_xhigh) then
         emc_corr = 1.0
      endif
      
      emc_func_xem= emc_corr
      return
      end
      

      integer*4 function last_char(string)
C+______________________________________________________________________________
!
! LAST_CHAR - Return the position of the last character in the string which
! is neither a space or a tab. Returns zero for null or empty strings.
C-______________________________________________________________________________

      implicit none
      integer*4 i
      character*(*) string
      character*1 sp/' '/
      character*1 tab/'	'/
      
      save
      
C ============================= Executable Code ================================

      last_char = 0
!     write(*,*) 'LAST CHAR ROUTINE INPUT: ',string
!     write(*,*) 'LAST CHAR ROUTINE SP and TAB: ',sp,tab
      do i = 1,len(string)
         if (string(i:i).ne.sp.and.string(i:i).ne.tab) last_char = i
!     write(*,*) 'LAST CHAR ROUTINE LOOP: ',string(i:i),i
      enddo
!     write(*,*) 'LAST CHAR ROUTINE RESULT: ',last_char
      
      return
      end

      recursive subroutine load_parameters(a, z)
      implicit none

      COMMON /Y_SCALING/ ag, bg, bigB, f0, alpha1, pfermi, eps
      real*8 ag, bg, bigB, f0, alpha1, pfermi, eps
      COMMON /JOHNS_FAC/ j_fac
      real*8 j_fac
      COMMON /QE_TAIL_COR/ tc_xlow, tc_xhigh, tc_aa, tc_bb, tc_cc, 
     >       tc_dd, tc_ee, tc_ff, tc_const
      real*8 tc_xlow, tc_xhigh, tc_aa, tc_bb, tc_cc, 
     >       tc_dd, tc_ee, tc_ff, tc_const
      COMMON /EMC_COR/ emc_xlow, emc_xhigh, emc_0, emc_1, emc_2, emc_3,
     >       emc_4, emc_5
      real*8 emc_xlow, emc_xhigh, emc_0, emc_1, emc_2, emc_3,
     >       emc_4, emc_5
      COMMON /DHX_COR/ dhx_xlow,dhx_xhigh,dhx_cor_min,dhx_cor_xalt,
     >       dhx_cor_alt_val, dhx_cor_0, dhx_cor_1
      real*8 dhx_xlow,dhx_xhigh,dhx_cor_min,dhx_cor_xalt,
     >       dhx_cor_alt_val, dhx_cor_0, dhx_cor_1

      logical found, first/.true./
      integer*4 i,j
      real*8 a, z
      real*8 parmlist(32)

      include 'xem_parameters.inc'

C Get target specific stuff from lookup table.

      found = .false.
      do i = 1,10               !loop over known targets.
         if (lookup(i,1).eq.z.and.float(int(a)).eq.lookup(i,2)) then
	    found = .true.
	    do j = 1,32
               parmlist(j) = lookup(i,j+3)
	    enddo
         endif
      enddo

<<<<<<< f7e0ffcda06cc96c2eb90ee2c9c812056d04ae12
      EPS = parmList(2)
      f0  =  parmList(3)
      bigB =  parmList(4)
      ag =  parmList(5)
      bg =  parmList(6)
      alpha1 =  parmList(7)

      j_fac = parmList(8)

      tc_const = parmList(9)
      tc_xlow = parmList(10)
      tc_xhigh = parmList(11)
      tc_aa = parmList(12)
      tc_bb = parmList(13)
      tc_cc = parmList(14)
      tc_dd = parmList(15)
      tc_ee = parmList(16)
      tc_ff = parmList(17)

      emc_xlow  = parmList(18)
      emc_xhigh = parmList(19)
      emc_0     = parmList(20)
      emc_1     = parmList(21)
      emc_2     = parmList(22)
      emc_3     = parmList(23)
      emc_4     = parmList(24)
      emc_5     = parmList(25)

      dhx_xlow        = parmList(26)
      dhx_xhigh       = parmList(27)
      dhx_cor_min     = parmList(28)
      dhx_cor_0       = parmList(29)
      dhx_cor_1       = parmList(30)
      dhx_cor_xalt    = parmList(31)
      dhx_cor_alt_val = parmList(32)

c      write(6,*) EPS, f0, bigB, ag, bg, alpha1

c      write(6,*) emc_xlow, emc_xhigh, emc_0, emc_1, emc_2, emc_3, emc_4,
c     + emc_5

c      write(6,*) dhx_xlow, dhx_xhigh, dhx_cor_min, dhx_cor_0, dhx_cor_1,
c     + dhx_cor_xalt, dhx_cor_alt_val

      if (.not.found.and.first) then
         write(6,*) 'cant find target in lookup table!'
         write(6,*) 'Using Carbon parameters with given A and Z...'
         first=.false.
         call load_parameters(dble(12.),dble(6.))
         return			!Quit if couldn't find info.
      else if(.not.found.and. .not. first) then
         write(6,*) 'Something is wrong with default load!'
=======
      real*8 function tail_cor(x,a)
      implicit none
      integer a
      real*8  x, aa, bb, cc, dd, ee, ff

      aa=1.
      bb=0.
      cc=0.
      dd=0.
      ee=0.
      ff=0.
      if(a.eq.2) then
        aa =  1.72816025139459
        bb =  2.53114253906281
        cc = -2.72505067059703
        dd = -1.58637090989274
        ee = -16.3972900673533

      elseif(a.eq.3) then

        bb =  0.8
        cc =  0.06864880408328
        dd = -0.320972192919132
        ee =  0.0
        aa =  0.552199789237622

      elseif(a.eq.4) then
        bb = 0.466102123780537
        cc = 0.0156369553828335
        dd = -0.122243059123825
        aa = 0.682462282515971
      elseif(a.eq.9) then
        bb = 0.463011918692135
        cc = 0.0125252624899601
        dd = -0.101843839213646
        aa = 0.674455752091906
      elseif (a.eq.12) then
        bb = 0.222436834975864
        cc = 0.00769270345172033
        dd = -0.060282702596254
        aa = 0.840262866196151
      elseif((a.eq.63).or.(a.eq.64).or.(a.eq.56)) then
        bb = 0.041323394008416
        cc = 0.00447016532137148
        dd = -0.0303635977582275
        aa = 1.00675406673173
      else if (a.eq.197) then
        bb = 0.0667337559531751
        cc = 0.00448892579200859
        dd = -0.0334460588480325
        aa = 0.981274819686673
      else
        write(6,*) ' tail cor, do not want ', a
>>>>>>> Add most updated xem_model fortran code.  Abstract reading of E0, EP, THETA, A, Z into xem_model.
      endif
      
      return
      end

<<<<<<< f7e0ffcda06cc96c2eb90ee2c9c812056d04ae12
      subroutine write_corrections()
      implicit none
      COMMON /CORRECTIONS/ johns_fac, fact, corfact, dhxcorfac, emc_corr
      real*8 johns_fac, fact, corfact, dhxcorfac, emc_corr
      logical first/.true./

      if(first) then
         first=.false.
c         write(6,*) "johns_fac fact*1000. corfact dhxcorfac emc_corr"
      endif
c      write(6,*) johns_fac, fact*1000., corfact, dhxcorfac, emc_corr
      return
      end
=======
!this is the postthesis iteration of daves
	real*8 function emc_func_xem(x,A) ! now compute the emc effect from our own fits.
	implicit none
        integer a
        real*8 x,xtmp
	real*8 emc

	if(x.le.0.9) then
	   xtmp = x
	else
	   xtmp = 0.9
	endif

CDeuterium******************************************************************
	if(A.eq.2) then
C iteration 1 Casey
           emc = 1

CHe3************************************************************************
	else if(A.eq.3) then
C iteration 1 Casey
c           emc = 1
           emc = 0.81349 + 4.17691*xtmp -21.00865*xtmp**2 + 
	1        50.24201*xtmp**3 -57.67029*xtmp**4 + 25.22827*xtmp**5

CHe4************************************************************************
	else if(A.eq.4) then
C iteration 1 Casey
c           emc = 1
           emc = 0.76458 + 4.15387*xtmp -20.56000*xtmp**2 + 
	1        49.19713*xtmp**3 -56.55423*xtmp**4 + 24.57040*xtmp**5

C Be************************************************************************
	else if(A.eq.9) then
C iteration 1 Casey
c           emc = 1
           emc = 0.75577 + 6.31620*xtmp -31.22099*xtmp**2 + 
	1        73.30890*xtmp**3 -83.58729*xtmp**4 + 36.42793*xtmp**5

C Carbon********************************************************************
	else if(A.eq.12) then
C iteration 1 Casey
c           emc = 1
           emc = 0.70352 + 6.44515*xtmp -31.99850*xtmp**2 + 
	1        75.42336*xtmp**3 -86.47596*xtmp**4 + 37.81601*xtmp**5

C Al************************************************************************
	else if(a.eq.27) then
C iteration 1 Casey
c           emc = 1
           emc = 0.63596 + 6.59120*xtmp -32.87832*xtmp**2 + 
	1        76.97413*xtmp**3 -88.46550*xtmp**4 + 39.04195*xtmp**5

C Copper********************************************************************
	else if(A.eq.64) then
C iteration 1 Casey
c           emc = 1
           emc = 0.62897 + 9.34351*xtmp -46.33351*xtmp**2 + 
	1        107.57064*xtmp**3 -122.72402*xtmp**4 + 54.02095*xtmp**5

C Copper********************************************************************
        else if(A.eq.56) then
C iteration 1 Casey
      emc = 1

C Gold**********************************************************************
	else if(A.eq.197) then
C iteration 1 Casey
c           emc = 1
           emc = 0.59609 + 11.06185*xtmp -54.73150*xtmp**2 + 
	1        125.84238*xtmp**3 -142.97502*xtmp**4 + 63.04760*xtmp**5

	else
	   write(*,*) '** in emc_func_xem, unknown target',a
	   stop
	endif

	emc_func_xem= emc
	return
	end
>>>>>>> Add most updated xem_model fortran code.  Abstract reading of E0, EP, THETA, A, Z into xem_model.

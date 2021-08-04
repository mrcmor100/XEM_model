CAM      program run_xem_model
CAM      implicit none
CAM      integer eof
CAM      character*80 rawname, filename, ofilename
CAM      real*8 E0, EP, THETA, A, Z, dis_XS, qe_XS
CAM      integer*4	last_char
CAM
CAM      read(*,1968) rawname
CAM 1968 format(a)
CAMCAM Open input file.
CAM      filename = 'input/'//rawname(1:last_char(rawname))//'.inp'
CAM      open(unit=11,status='old',file=filename)
CAM
CAMCAM Open output file.
CAM      filename = './output/'//rawname(1:last_char(rawname))//'.out'
CAM      open (unit=14,status='unknown',file=filename)
CAM
CAM 40   READ(11,*,IOSTAT=EOF) E0, EP, THETA, A, Z
CAM      if(eof.ge.0) then
CAM        call xem_model(E0, EP, THETA, A, Z, dis_XS, qe_XS)
CAM        goto 40
CAM      endif !if file has stuff
CAM      close(11);
CAM      end
CAM

!-----------------------------------------------------------------------------
CAM Main subroutine used in externals.
!-----------------------------------------------------------------------------
      subroutine SIGMODEL_CALC(E0in,EPin,THETAin,iZ,
     >        iA,avgM,DUM1,DUM2,SIGMApass,XFLAG,FACT)
      implicit none
      include 'constants_dble.inc'
      REAL E0in, EPin, THETAin
      REAL DUM1, DUM2, fact, avgM, SIGMApass
      real*8 sigma_send
      real*8 E0, EP, THETA, QSQ, NU, THR
      real*8 Y, A, Z, X, N
      real*8 INNP, INNT
      real*8 sig_qe_new, sigdis_new
      real*8 YSCALE
      integer iZ, iA, xflag

      
      include 'model_constants.inc'

      E0 = real(E0in,8)
      EP = real(EPin,8)
      THETA = real(THETAin,8)
      A = dble(iA)
      Z = dble(iZ)

      INNT = 30
      INNP = 30
      N = A - Z
      THR = THETA*pi/180.
      QSQ = 4*E0*EP*(SIN(THR/2))**2
      NU  = E0 - EP
      X   = QSQ/(2*nuc_mass*NU)
      Y   = YSCALE(E0,EP,THR,A,EPS(iA))

      sig_qe_new = 0.0
      sigdis_new = 0.0
c      write(*,*) 'Y: ', Y
      IF(Y.lt.1.E19) then
         if((XFLAG.eq.1).or.(XFLAG.eq.2)) then
            call sig_qe_calc(Y,iA, iZ,E0,EP, THR, X, sig_qe_new)
            sig_qe_new = sig_qe_new*1000000.0
         endif
         if((XFLAG.eq.1).or.(XFLAG.eq.3)) then
            call sig_dis_calc(iA, iZ, E0, EP, THR, Y, sigdis_new)
            sigdis_new = sigdis_new
         endif
      ELSE
         sigdis_new = 0.0
         sig_qe_new = 0.0
      ENDIF
c 2003 format (2(E13.5,1x))
c      write(14,2003) sigdis_new,sig_qe_new
      sigma_send = sig_qe_new + sigdis_new
      SIGMApass = real(sigma_send,4)
      return
      end

C____________________________________________________________________________

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

 2002 format (10(E13.5,1x))
      if(y.ne.0) then
         write(14,2002) y,a,z, THETA,EP,x ,sigdis_new,sig_qe_new
      else       !comment out under normal circumstances
         write(14,2002) y,a,z, THETA,EP,x ,sigdis_new,sig_qe_new
      endif

      call write_corrections()

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
      WSQ = -qsq + nuc_mass**2 + 2.0*nuc_mass*nu

      PMAX=1.0

      if(a.ge.2) then
         if(WSQ.lt.2.25) then
            innt=30
            innp=30
         else
            innt=30
            innp=10
         endif

         CALL smear4all(E,EP,THR,dble(A),dble(Z),dble(N),eps,PMAX,
     +        dble(INNP),dble(INNT),f0,bigB,ag,bg,
     +        alpha1 ,SIGDIS)

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
      endif
      
      return
      end

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

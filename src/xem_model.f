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
      include 'constants_dble.inc'
      real*8 E0, EP, THETA, QSQ, NU, THR
      real*8 Y, A, Z, X, N
      real*8 INNP, INNT
      real*8 sig_qe_new, sigdis_new
      real*8 YSCALE

      include 'model_constants.inc'

      INNT = 30
      INNP = 30

      N = A - Z
      THR = THETA*pi/180.
      QSQ = 4*E0*EP*(SIN(THR/2))**2
      NU  = E0 - EP
      X   = QSQ/(2*nuc_mass*NU)
      Y   = YSCALE(E0,EP,THR,A,EPS(int(A)))
      !Check if Y is right kinda..
      if(Y.lt.1.E19) then 
         call sig_dis_calc(int(A), int(Z), E0, EP, THR, Y, sigdis_new)
         call sig_qe_calc(Y, int(A), int(Z),E0,EP, THR, X, sig_qe_new)
      else
         sigdis_new = 0.0
         sig_qe_new = 0.0
      endif
      write(236,*) A, Z, X, Y

 2002 format (10(E13.5,1x))
      if(y.ne.0) then
         write(14,2002) y,a,z, THETA,EP,x ,sigdis_new,sig_qe_new
      else       !comment out under normal circumstances
         write(14,2002) y,a,z, THETA,EP,x ,sigdis_new,sig_qe_new
      endif

      return
      end

!-----------------------------------------------------------------------------
CAM Calculate DIS using donal's smearing routine.
!-----------------------------------------------------------------------------
      subroutine sig_dis_calc(a, z, e, ep, thr, y, sigdis)
      implicit none
      include 'constants_dble.inc'
      integer a, z, n, innp, innt
      real*8 e, ep, thr, wsq, x,y, qsq, nu
      real*8 sigdis, corfac, pmax
      real*8 f1, f2, w1, w2, tan_2, sig_mott
      real*8 frac, x1, x2, emc_corr
      real*8 emc_func_xem
      include 'model_constants.inc'

      n=a-z
      qsq=4.0*e*ep*sin(thr/2)**2
      nu=e-ep
      x=qsq/2.0/nuc_mass/nu
      WSQ = -qsq + nuc_mass**2 + 2.0*nuc_mass*nu

      PMAX=1.0
      x1=0.8 !x<x1 --> use emc corrected ld2
      x2=0.9 !x1<=x<x2 --> smooth transition from emc corrected ld2 to donals smearing

      if(a.ge.2) then
         if(WSQ.lt.2.25) then
            innt=30
            innp=30
         else
            innt=30
            innp=10
         endif

         CALL smear4all(E,EP,THR,dble(A),dble(Z),dble(N),EPS(a),PMAX,
     +        dble(INNP),dble(INNT),f0(a),bigB(a),ag(a),bg(a),
     +        alpha1(a) ,SIGDIS)

CAM EMC_FUNC specific correction.
CAM D2 set to F1F221.
         if (x.lt.x1) then
            emc_corr = emc_func_xem(x,a)
         elseif ((x.ge.x1).and.(x.lt.x2)) then
            frac = (x-x1)/(x2-x1)
            emc_corr = 1.0*frac + emc_func_xem(x,a)*(1.-frac)
CAM Uncomment when setting emc_func_xem
CAM           emc_corr = emc_func_xem(x,a)
         elseif(x.ge.x2) then
            emc_corr = 1.0
         endif
         sigdis = sigdis*emc_corr

CAM Make a correction to the high_x tail.
         corfac=1.
         if ((x.gt.0.9)) then
            call  dis_highx_cor(a,x,corfac)
            sigdis = sigdis*corfac
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
      integer a, z
      real*8 thr, y, e, ep, x, nu, q4_sq, q3v
      real*8 sig_qe, sigma_n, sigma_p
      real*8 dwdy, fact
      real*8 johns_fac
      real*8 x_cor, my_frac, x_high, x_low, corfact, tail_cor
      include 'model_constants.inc'

      q4_sq = 4.*e*ep*sin(thr/2.)**2 !4-mom. transf. squared.
      nu = e-ep                 !Energy loss
      q3v = sqrt(q4_sq + nu**2) !3-momentum transfer
      dwdy = q3v/sqrt(nuc_mass**2+q3v**2+y**2+2.0*q3v*y)
      x=q4_sq/2.0/nuc_mass/nu

CAM F(y) from Nadia's Thesis based on:
CAM https://arxiv.org/pdf/nucl-th/9702009.pdf
      y=y*1000
      if(a.eq.2.) then
        sig_qe = (f0(a)-bigB(a))*alpha1(a)**2*exp(-(ag(a)*y)**2)
     +    /(alpha1(a)**2+y**2)+bigB(a)*exp(-bg(a)*abs(y))
      else
        sig_qe = (f0(a)-bigB(a))*alpha1(a)**2*exp(-(ag(a)*y)**2)
     +    /(alpha1(a)**2+y**2)+bigB(a)*exp(-(bg(a)*y)**2)
      endif
      y=y/1000.

CAM JOHNS_FAC basically linear factor less than X=1
      johns_fac=1.
      if (a.eq.3) then
         johns_fac=max(1.,1.+y*1.4*2.5)
      elseif (a.eq.4) then
         johns_fac=max(1.,1.+y*1.4*1.75)
      elseif (a.eq.10) then
         johns_fac=max(1.,1.+y*1.4*2.5)
      elseif (a.eq.11) then
         johns_fac=max(1.,1.+y*1.4*2.5)
      elseif (a.eq.12) then
         johns_fac=max(1.,1.+y*1.4*2.5)
!         write(*,*) 'johns fac:',y, x, johns_fac
      elseif ((a.eq.56).or.(a.eq.64)) then
         johns_fac=max(1.,1.+y*1.4*3)
      elseif (a.eq.197) then
         johns_fac=max(1.,1.+y*1.4*4)
      endif
      sig_qe=johns_fac*sig_qe

CAM Make get the off-shell contributions based on DeForest:
CAM T. De Forest, Jr. Nuc. Phys. A392 232 (1983)
      call sig_bar_df(e, ep, thr*180/pi, y,dble(0.0), sigma_p, sigma_n)

      fact = (Z*sigma_p+(A-Z)*sigma_n)/dwdy
      sig_qe=sig_qe*fact*1000.

CAM tail_cor applied here
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
      if((x.ge.(x_low))) then   !.and.(A.eq.2)) then
         corfact=tail_cor(x_cor,a)
         if(x.lt.(x_high)) then
	    my_frac=(x-x_low)/(x_high-x_low)
            corfact=my_frac*corfact+(1.-my_frac)
         endif
         sig_qe=sig_qe*corfact
      endif
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

      subroutine sig_bar_df(e1,e2,theta,pl,pt,sig_p,sig_n)
C+______________________________________________________________________________
!
! DESCRIPTION:
!
!   Calculate the DeForest based off-shell cross sections averaged over the
!   angle PHI, which is the angle between the electron scattering plane and
!   the plane containing the initial and final momentum of the struck nucleon.
!
!   We use CIT form factors, based on the four-momentum transfer.
!
!   All energy/momentum variables are in GeV; cross sections are in mB/ster.
!
! ARGUMENTS: (All R*4 data type)
!
!   E1:    Energy of incident electron.
!   E2:    Energy of scattered electron.
!   THETA: Electron scattering angle in deg.
!   PL:    Projection of init. state nucleon's momentum along 3-mom. transfer.
!   PT:    Component of init. state nucleons momentum perp. to 3-mom. trnsfr.
!   SIG_P: Cross section for proton.  (OUTPUT)
!   SIG_N: Cross section for neutron. (OUTPUT)
C-______________________________________________________________________________

      implicit none

      include		'math_physics.inc'

C     Argument data types.

      real*8		e1,e2,theta,pl,pt,sig_p,sig_n
      real*8		gep,gmp,gen,gmn, th, tan_2,nu
      real*8            q4_2, tau, qv_2, qv, pt_2,p_2
      real*8            sig_mott, en_f, e_bar,q4_bar_2
      real*8            tau_bar,tmp, f1p, f1n,f2n,f2p
      real*8            qrat, t1, t2,q4_sq, jerk

C     ============================ Executable Code =====================
C     ============

C     Compute some kinematics.

      th = theta*d_r
 !     th=theta
      tan_2 = tan(th/2)**2
      nu = e1 - e2
      q4_2 = 4*e1*e2*sin(th/2)**2
      tau = q4_2/(4*m_p**2)
      qv_2 = q4_2 + nu**2
      qv = sqrt(qv_2)
      pt_2 = pt*pt
      p_2 = pl**2 + pt_2

!     Mott cross section in mB/sR.

      sig_mott = hc_2*(alpha*cos(th/2))**2/(2*e1*sin(th/2)**2)**2

!     write(*,*) 'sig mott is ', sig_mott

!     Final energy of struck nucleon, assuming it is on shell.

      en_f = sqrt(m_p**2+qv_2+p_2+2*qv*pl)

C     'BAR'-ed quantities of DeForest.

      e_bar = sqrt(p_2 + m_p**2)
      q4_bar_2 = qv_2 - (en_f - e_bar)**2
      tau_bar = q4_bar_2/(4*m_p**2)

!	write(*,*) 'next set ', en_f, e_bar, q4_bar_2, tau_bar, q4_2
C Get form factors. Convert to F1, F2.

!	call nuc_form(real(q4_2),gep,gmp,gen,gmn)

!	gep=0.0
!	gen=0.0
!	gmp=0.0
!	gmn=0.0
	q4_sq=q4_2
        jerk=12.0
!        call nform (dble(14.0),q4_sq,gep,gen,gmp,gmn)
        call nform (dble(12.0),q4_sq,gep,gen,gmp,gmn)
        tmp = gmp

!	write(*,*) 'ff in john ', gep, gen, gmp, gmn
!        call nform (1.0,q4_sq,gep,gen,gmp,gmn)
!        write(28,*) 'sanity check ', gmp,q4_sq,e1,tmp
!        gmp = tmp

c	write(28,*) 'ff in john 2 ', e2,gep, gen, gmp, gmn,tmp
	f1p = (gep + tau*gmp)/(1 + tau)
	f1n = (gen + tau*gmn)/(1 + tau)
	f2p = (gmp - gep)/(1 + tau)
	f2n = (gmn - gen)/(1 + tau)

!	write(*,*) 'others are ', f1p, f1n, f2p, f2n

C DeForest cross sections.

	qrat = q4_2/qv_2
	t1 = q4_bar_2*tan_2/2 + .25*qrat*(q4_bar_2 - q4_2)
	t2 = .25*(qrat**2)*(e_bar + en_f)**2 + (qrat/2 + tan_2)*pt_2
	sig_p = (t1*(f1p + f2p)**2 + t2*(f1p**2 + tau_bar*f2p**2))
	sig_n = (t1*(f1n + f2n)**2 + t2*(f1n**2 + tau_bar*f2n**2))
	sig_p = sig_mott*sig_p/(en_f*e_bar)
	sig_n = sig_mott*sig_n/(en_f*e_bar)

!        write(28,*) 'returning from sigbar', e2, th,sig_p,sig_n, pl,pt
        return
        end

!-----------------------------------------------------------------------------------------------------

	subroutine dis_highx_cor(anuc,x,cor)
        implicit none
        integer anuc
	real*8 x,cor, frac,xlow1,xhigh1

	xlow1=0.9
	xhigh1=0.95

	frac=1.
	cor=1.
	if(anuc.eq.2) then
	  cor=1.
          cor=-3.30482*x+ 4.10442
          if((x.ge.xlow1).and.(x.le.xhigh1)) then
            frac = (x-xlow1)/(xhigh1-xlow1)
          endif
          cor=frac*cor+1.-frac
c          write(*,*) 'No call to dis_highx_cor btwn xlow, xhigh in D2,
c     &     cor:', cor
c          write(*,*) 'No call to dis_highx_cor above x=0.9 in D2,
c     &     cor:', cor
	elseif (anuc.eq.3) then
	  xlow1=1.
	  xhigh1=1.15
	  if (x.gt.xlow1) then
	    if (x.lt.1.15) then
	      cor=-4.80303*x+ 5.74758
	    else
	      cor=0.5
	    endif
	    if((x.ge.xlow1).and.(x.le.xhigh1)) then
	      frac = (x-xlow1)/(xhigh1-xlow1)
	    endif
	  endif
          cor=frac*cor+1.-frac
	 elseif(anuc.eq.4) then
           cor=-1.78944*x+ 2.63272
	   if((x.ge.xlow1).and.(x.le.xhigh1)) then
	     frac = (x-xlow1)/(xhigh1-xlow1)
	   endif
	   cor=frac*cor+1.-frac
	 elseif(anuc.eq.9) then
	   cor=-1.7549631060248*x+ 2.6646629067298
	   if((x.ge.xlow1).and.(x.le.xhigh1)) then
	     frac = (x-xlow1)/(xhigh1-xlow1)
	   endif
	   cor=frac*cor+1.-frac

	 elseif(anuc.eq.12) then
	   cor=-1.29213*x+ 2.2087
	   if((x.ge.xlow1).and.(x.le.xhigh1)) then
	     frac = (x-xlow1)/(xhigh1-xlow1)
	   endif
	   cor=frac*cor+1.-frac
	 elseif((anuc.eq.63).or.(anuc.eq.64)) then
	   cor=-1.65829142599487*x+ 2.48174872208596
	   if((x.ge.xlow1).and.(x.le.xhigh1)) then
	     frac = (x-xlow1)/(xhigh1-xlow1)
	   endif
	   cor=frac*cor+1.-frac
	 elseif(anuc.eq.197) then
	   cor=-1.42430013496752*x+ 2.25789690593227

	   if((x.ge.xlow1).and.(x.le.xhigh1)) then
	     frac = (x-xlow1)/(xhigh1-xlow1)
	   endif
	   cor=frac*cor+1.-frac

	 else
	   cor=1.
	 endif

	 if(cor.lt.0.4) cor=0.4

	return
	end

c--------------------------------------------------------------------------------
CAM tail_cor used to fix tail of QE cross-section.
c--------------------------------------------------------------------------------
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
      elseif (a.eq.10) then
        bb = 0.222436834975864
        cc = 0.00769270345172033
        dd = -0.060282702596254
        aa = 0.840262866196151
      elseif (a.eq.11) then
        bb = 0.222436834975864
        cc = 0.00769270345172033
        dd = -0.060282702596254
        aa = 0.840262866196151
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
      endif

      tail_cor=(aa*exp(bb*x)+cc*x**6+dd*x**4+ee*x**2+ff)
      return
      end

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
c Parameters for F1F221
           emc = 1
           emc = 0.98835 + 0.15439*xtmp - 0.34950*xtmp**2 +
	1        0.34427*xtmp**3 - 0.02874*xtmp**4 - 0.04936*xtmp**5

CHe3************************************************************************
	else if(A.eq.3) then
C iteration 1 Casey
           emc = 1
           emc = 1.31789 - 2.47197*xtmp + 12.03740*xtmp**2 -
	1        26.33859*xtmp**3 + 25.66942*xtmp**4 - 9.14186*xtmp**5

CHe4************************************************************************
	else if(A.eq.4) then
C iteration 1 Casey
           emc = 1
           emc = 1.36112 - 3.41330*xtmp + 15.67047*xtmp**2 -
	1        32.79506*xtmp**3 + 31.50846*xtmp**4 - 11.52059*xtmp**5

C Be************************************************************************
	else if(A.eq.9) then
C iteration 1 Casey
           emc = 1
           emc = 1.52438 - 3.26366*xtmp + 13.78482*xtmp**2 -
	1        26.95815*xtmp**3 + 22.72433*xtmp**4 - 6.66419*xtmp**5

C Carbon********************************************************************
	else if(A.eq.10) then
C iteration 1 Casey
           emc = 1
           emc = 1.42308 - 2.65087*xtmp + 11.41047*xtmp**2 -
	1        22.54747*xtmp**3 + 18.47078*xtmp**4 - 5.07537*xtmp**5

C Carbon********************************************************************
	else if(A.eq.11) then
C iteration 1 Casey
           emc = 1
           emc = 1.42308 - 2.65087*xtmp + 11.41047*xtmp**2 -
	1        22.54747*xtmp**3 + 18.47078*xtmp**4 - 5.07537*xtmp**5

C Carbon********************************************************************
	else if(A.eq.12) then
C iteration 1 Casey
           emc = 1
           emc = 1.42308 - 2.65087*xtmp + 11.41047*xtmp**2 -
	1        22.54747*xtmp**3 + 18.47078*xtmp**4 - 5.07537*xtmp**5

C Al************************************************************************
	else if(a.eq.27) then
C iteration 1 Casey
           emc = 1
           emc = 1.29636 - 1.70678*xtmp + 6.49121*xtmp**2 -
	1        11.52374*xtmp**3 + 6.07310*xtmp**4 + 0.48219*xtmp**5

C Copper********************************************************************
	else if(A.eq.64) then
C iteration 1 Casey
           emc = 1
           emc = 1.41180 - 0.41740*xtmp - 0.41459*xtmp**2 +
	1        5.07170*xtmp**3 - 13.84335*xtmp**4 + 9.81772*xtmp**5

C Copper********************************************************************
        else if(A.eq.56) then
C iteration 1 Casey
           write(*,*) 'No fit for EMC on Iron.'
           emc = 1

C Gold**********************************************************************
	else if(A.eq.197) then
C iteration 1 Casey
           emc = 1
           emc = 1.32738 + 2.06855*xtmp - 13.01543*xtmp**2 +
	1        33.54028*xtmp**3 - 45.36383*xtmp**4 + 23.48632*xtmp**5

	else
	   write(*,*) '** in emc_func_xem, unknown target',a
	   stop
	endif

	emc_func_xem= emc
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
!	write(*,*) 'LAST CHAR ROUTINE INPUT: ',string
!	write(*,*) 'LAST CHAR ROUTINE SP and TAB: ',sp,tab
	do i = 1,len(string)
	   if (string(i:i).ne.sp.and.string(i:i).ne.tab) last_char = i
!	write(*,*) 'LAST CHAR ROUTINE LOOP: ',string(i:i),i
	enddo
!	write(*,*) 'LAST CHAR ROUTINE RESULT: ',last_char

	return
	end

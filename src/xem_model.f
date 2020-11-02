      program xem_model_with_special_sauce
      implicit none
      include 'constants_dble.inc'
      real*8 e, ep, theta, q3v, qsq,nu,thr,cs,cs_err,k
      real*8 sys_err, beta2, tempo, sixGeVData, sixGevDataErr
      real*8 cs_elastic,sig_el,f_y,y,eps(197),yscale,a,z
      real*8 pmax,sigdis!,ag,bg,cg,dg
      real*8 innp, innt, n,fy_err,f2_err,xsi,yalt
      real*8 pfermi(197),R,beta,xsec_mott,cs_mott,f2,x
      real*8 sigma_p, sigma_n, y_pt
      real*8 psi, psi_p, psi_pp,sig_long
      real*8 ag(197), bg(197), bigB(197), f0(197), alpha1(197)
      !New for xem hack
      real*8 ld2_a,ld2_z,ld2_nn,ld2_aux(7),ld2_sig_dis
      real*8 emc_corr,ld2_inta, sigdis_raw
      real*8 x1,x2,sig_dis_emc,sig_donal,sig_dis_donal
      real*8 frac,corfac, x_low, x_high
      real*8 emc_func_xem, sig_qe, f2_fac, label,wsq,cs_elastic_r

      real*8 emc_func_slac, corfact
      real*8 my_frac, sig_before, glob_cor, m_atom(197)
      real*8 m1, yorig, fact, dwdy,q4_sq, tail_cor,x_cor,johns_fac
      real*8 alpha_tn,sign_d,sigp_d, sigep_d, sigen_d,sigm_d,sigm_me

      real*8 err_e, err_ep, err_th, err_cc, err_dummy

      real*8 sig_qe_new, sigdis_new

      integer eof,i

      open(unit=14,file='output/fort.14', status='REPLACE')
!      open(unit=22,file='output/fort.22', status='REPLACE')
!      open(unit=24,file='output/fort.24', status='REPLACE')
!      open(unit=28,file='output/fort.28', status='REPLACE')

      do i=1,197
         pfermi(i)=.2
         EPS(i)=0.01

      enddo

      pfermi(197)=.264
      pfermi(64)=.260
      pfermi(56)=.260
      pfermi(27)=.240
      pfermi(12)=.221
      pfermi(9)=.230
      pfermi(4)=.160
      pfermi(3)=.160
      pfermi(2)=.055
!      eps(197)=.049
      eps(197)=.0058
      eps(64)=0.010
      eps(56)=0.010
      eps(27)=.0083
      eps(12)=.016
      eps(9)=0.0168
      eps(4)=.02
      eps(3)=.0055
      eps(2)=.00225

      m_atom(1)=1.00794*.9314943
      m_atom(2)=2.0141*.9314943
      m_atom(3)=3.0160*.9314943
      m_atom(4)=4.002602*.9314943
      m_atom(9)=9.012182*.9314943
      m_atom(12)=12.0110*.9314943
      m_atom(27)=26.98*.9314943
      m_atom(56)=55.8470*.9314943
      m_atom(64)=63.546*.9314943
      m_atom(197)=196.9665*.9314943

      f0(2)= 0.0087421570999812/.96689
      f0(3)= 0.00530944208609402
      f0(4)= 0.00401973252542711
      f0(9)= 0.00348142965010112
      f0(12)= 0.00318223097961272
      f0(27)=0.0032783
      f0(56)=0.00289
      f0(64)= 0.002874027535171
      f0(197)= 0.00264244988981799

      bigB(2)= 0.000823944541945584/.96689
      bigB(3)= 0.00218428711503937
      bigB(4)= 0.00134486977748062
      bigB(9)= 0.00116084580978475
      bigB(12)= 0.00135911193714437
      bigB(27)=0.0013474
      bigB(56)=0.0014016
      bigB(64)= 0.000886615094402987
      bigB(197)= 0.000763159585000613

      ag(2)= 0.00772721685529826
      ag(3)= 0.00288636824154895
      ag(4)= 0.0026986223386714
      ag(9)= 0.00311951067940091
      ag(12)= 0.00302654042944761
      ag(27)=0.0029698
      ag(56)=0.0031802
      ag(64)= 0.00309586042619543
      ag(197)= 0.00306539017228031

      bg(2)= 0.00939373447970209
      bg(3)= 0.0103491992909402
      bg(4)= 0.00749410002149716
      bg(9)= 0.00783976750998206
      bg(12)= 0.00705046729168206
      bg(27)=0.006576
      bg(56)=0.0072635
      bg(64)= 0.00709447514594065
      bg(197)= 0.00676775547526636

      alpha1(2)= 45.3840728438637
      alpha1(3)= 64.247220058789
      alpha1(4)= 100.255815258545
      alpha1(9)= 110.966740959374
      alpha1(12)= 137.284604974981
      alpha1(27)=131.845
      alpha1(56)=165.7
      alpha1(64)= 132.45766614518
      alpha1(197)= 132.451655636968

      ld2_a=2.00
      ld2_z=1.00
      ld2_nn=1.00

      ld2_aux(3)= 0.008742157099981
      ld2_aux(4)= 0.00082394
      ld2_aux(5)= 0.0077272168552982
      ld2_aux(6)= 0.00939373447970209
      ld2_aux(7)= 45.3840728438637
      ld2_inta=ld2_a

      x1=0.8                    !x<x1 --> use emc corrected ld2
      x2=0.9                    !x1<=x<x2 --> smooth transition from emc corrected ld2 to donals smearing
!      e=5.766!15

      open(unit=11,file='./input/compare_d_h.dat', status
     +  ='old')
      write(*,*) 'opened file'

 40   READ(11,*,IOSTAT=EOF) e,ep,theta,a,z
        if(eof.ge.0) then
c           if (a.eq.64.0) write (*,*) 'doing copper' , e,ep,theta,a,z,cs,cs_err
c           if (z.eq.20.0) then
c              z=29.0
c              a=64.0
c           endif
           n=a-z
c        write(*,*) e,ep,theta,a,z
        PMAX = 2.5*PFERMI(int(a))
        INNT = 30
        INNP = 30

        thr=theta*pi/180.!+0.0003
        qsq=4*E*EP*(SIN(THR/2))**2
        nu=e-ep
        Q3V = sqrt(QSQ+NU**2)
        y=YSCALE(E,EP,THR,A,EPS(int(A)))
        m1=m_atom(int(a))


!        write(*,*) 'my two ys are', yorig, y
        K= Q3V/(sqrt(nuc_mass**2+(Y+Q3V)**2))
        R=0.32/qsq
!        R=0.5*R
        beta=2*(tan(thr/2))**2*(1+nu**2/qsq)/(1+R)
        xsec_mott=cs_mott(thr,e)
        f2=nu*cs/(xsec_mott*(1+beta))
        f2_fac=nu/(xsec_mott*(1+beta))
        beta2=beta
        f2_err=f2*(cs_err/cs)
        x=qsq/(2*nuc_mass*nu)
        xsi=2*x/(1+sqrt(1+4*nuc_mass**2*x**2/qsq))

!
        call sig_dis_calc(int(a), int(z), e, ep, thr, y,sigdis_new)
   !     write (*,*) 'NEW DIS HOTNESS', sigdis, sigdis_new

 2024  format (4(E13.5,1x))
c        write(24, 2024) x, sigdis_raw, emc_corr, corfac
  !      write(*,*) 'EMC RAW ', x, sigdis_raw, emc_corr, corfac
!     ENd of High X tweak.

       call sig_qe_calc(y, int(a), int(z),e,ep,thr, x, sig_qe_new)
 !      write(*,*) 'SUPERMODEL', sig_qe, sig_qe_new

 2002  format (10(E13.5,1x))
 2003  format (3(I4, 1x),13(E13.5,1x))
 2004  format (8(E13.5,1x))
 2005  format (3(I4, 1x),13(E13.5,1x))
      if(y.ne.0) then
        write(14,2002) y,a,z, theta,ep,x ,sigdis_new,sig_qe_new
      else   !comment out under normal circumstances
        write(14,2002) y,a,z, theta,ep,x ,sigdis_new,sig_qe_new
      endif

      goto 40
      endif !if file has stuff
      close(11);



      end

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

      subroutine sig_dis_calc(a, z, e, ep, thr, y,sigdis)
      implicit none
      include 'constants_dble.inc'
      integer a, z, n, innp, innt
      real*8 sig_dis, e, ep, thr, cs_elastic,wsq
      real*8 sigma_n, sigma_p, sigdis_raw
      real*8 x1, x2, emc_corr, x,y, qsq, nu
      real*8 frac, corfac,  emc_func_xem,sig_before, sigdis
      real*8 pmax, f1, f2, w1, w2, tan_2, cs_mott
      include 'model_constants.inc'

      qsq=4.0*e*ep*sin(thr/2)**2
      nu=e-ep
      x=qsq/2.0/nuc_mass/nu

      x1=0.8                    !x<x1 --> use emc corrected ld2
      x2=0.9                    !x1<=x<x2 --> smooth transition from emc corrected ld2 to donals smearing

      n=a-z
c      write (*,*) 'New DiS SUB'
      call sig_bar_df(e, ep, thr*180.0/pi, y,dble(0.0), sigma_p,
     +  sigma_n)

c      write(22,*) a, ep,x, sigma_n, sigma_p
c      write(*,*) 'SIGBAR  ',  a, ep,x, sigma_n, sigma_p

      cs_elastic=sigma_p*z+(a-z)*sigma_n
      cs_elastic=cs_elastic     !*1e9

      PMAX=1.0
!       Here's all the xem-specific model stuff.
      WSQ = -qsq + nuc_mass**2 + 2.0*nuc_mass*nu

      if(a.ge.2) then
         if(WSQ.lt.2.25) then
            innt=30
            innp=30
         else
            innt=30
            innp=10
         endif

c         write(*,*) 'New Sub, passing variables', E,EP,THR,A,Z,N,EPS(a),
c     +        PMAX,INNP,INNT,f0(a),bigB(a),ag(a),bg(a),alpha1(a)
         CALL BDISNEW4HE3(E,EP,THR,dble(A),dble(Z),dble(N),EPS(a),PMAX,
     +        dble(INNP),dble(INNT),f0(a),bigB(a),ag(a),bg(a),
     +        alpha1(a) ,SIGDIS)

c         write(*,*) 'NEW SUB bdis returned' , e, ep, sigdis, eps(a)

         sigdis_raw=sigdis

         if (x.lt.x1) then
            emc_corr = emc_func_xem(x,a)
         elseif ((x.ge.x1).and.(x.lt.x2)) then
            frac = (x-x1)/(x2-x1)
            emc_corr = 1.0*frac + emc_func_xem(x,a)*(1.-frac)
         elseif(x.ge.x2) then
            emc_corr = 1.0
         endif

         sigdis = sigdis*emc_corr

         if (sigdis.lt.0) then
            write(*,*) 'wtf 0', sigdis, x, corfac
         endif

         corfac=1.
         if ((x.gt.0.9)) then
            call  dis_highx_cor(a,x,corfac)
            sig_before=sigdis
            sigdis = sigdis*corfac
            if (sigdis.lt.0) then
               write(*,*) 'wtf', sigdis, x, corfac, sig_before
            endif
         endif

      else if(A.eq.1) then
         call F1F2IN20(dble(Z),dble(A), QSQ, WSQ, F1, F2)
         W1 = F1/.93827231D0
         W2 = F2/nu
         tan_2 = tan(thr/2)**2
C       Mott cross section
c         sigmott=(19732.0/(2.0*137.0388*e1cc*sn**2))**2*cs**2/1.d6
         sigdis =  cs_mott(thr,e)*(W2+2.0*W1*tan_2)
      endif

      end

!-----------------------------------------------------------------------------
      subroutine sig_qe_calc(y, a, z,e,ep,thr, x, sig_qe)
      implicit none
      include 'constants_dble.inc'
      integer a,z
      real*8 input, sig_qe, johns_fac
      real*8 q4_sq, q3v, dwdy,thr, fact,nu
      real*8 x, x_cor, corfact, my_frac, x_high, x_low
      real*8 sigma_n, sigma_p, theta
      real*8 y, e, ep, tail_cor
      include 'model_constants.inc'



c      write (*,*) 'this is fiiiiiine'


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
c        write(*,*) "FYONLY ", x, y, sig_qe, johns_fac
       sig_qe=johns_fac*sig_qe

  	q4_sq = 4.*e*ep*sin(thr/2.)**2 		!4-mom. transf. squared.
	nu = e-ep				!Energy loss
	q3v = sqrt(q4_sq + nu**2)		!3-momentum transfer
	dwdy = q3v/sqrt(nuc_mass**2+q3v**2+y**2+2.0*q3v*y)

        call sig_bar_df(e, ep, thr*180/pi, y,dble(0.0), sigma_p, sigma_n)
c        write (*,*) 'x, qsq', x, q4_sq, sigma_n, sigma_p


        fact = (Z*sigma_p+(A-Z)*sigma_n)/dwdy
        sig_qe=sig_qe*fact*1000.
c       write (*,*) 'Before COR ', sig_qe

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
c        write (*,*) 'Sending back, QE ', sig_qe
      end
!-----------------------------------------------------------------------------
	subroutine y_calc(e1,e2,theta,m1,es,y)



	implicit none

        logical y_calc_ok
	real*8 e1,e2,theta,m1,m_rec,es,y
	real*8 m2,w,th,q4_2,q2,q
	real*8 pp2
	real*8 temp,coeff_a,coeff_b,coeff_c,root2

	include 'constants_dble.inc'

C ============================ Executable Code =================================

	y_calc_ok = .false.			!assume failure

C Compute kinematics
        m_rec=nuc_mass
	m2 = m1 + es - m_rec			!Mass of (A-1) system
	w = e1 - e2	    			!Energy loss
!	th = theta*d_r				!Theta in radians
        th=theta
!        write(*,*) 'okay we have', e1, 21, theta,m1, es
	q4_2 = 4.*e1*e2*(sin(th/2.))**2.	!4-momentum transfer squared
	q2 = q4_2 + w*w				!3-momentum    "        "
	q  = sqrt(q2)

C Compute approximate value for k_perp.

! K_perp now suppressed since two parameters M_rec and E_sep make it
! entirely irrelevent.

C$$$	a = m1/m_p				!Approximate value for A.
C$$$	pf = .22*(1.-exp(-a/8.)) + .04
C$$$	pp2 = pf**2.

	pp2 = 0.				!Suppress k_perp.

        y=0.
C Compute terms in quad. formula.

	temp = q2 + m_rec**2 - m2*m2 - (m1 + w)**2
	coeff_a = 4.*(q2 - (m1 + w)**2.)
	coeff_b = 4.*q*temp
	coeff_c = temp**2. - 4.*(m2*m2 + pp2)*(m1 + w)**2.

C If no real solution, return ERROR.

	root2 = coeff_b**2. - 4.*coeff_a*coeff_c
	if (root2.lt.0..or.coeff_a.eq.0.) return

C Otherwise, return Y, and SUCCESS.

	y = (-coeff_b - sqrt(root2))/(2.*coeff_a)

	y_calc_ok = .true.

	end
!---------------------------------------------------------------------

      real*8 FUNCTION YSCALE(E,EP,THR,A,EPS)
      IMPLICIT NONE
      INCLUDE 'constants_dble.inc'
      real*8 E, EP,THR,A, NU,W,WP,AG,BG,PART1,PART2,PART3,EPS,QSQ,CG
      real*8 backup,RAD
      real*8 qv2


      yscale = 0.0
                                !write (6,*) 'in y scale'
                                !write (15,*) ' nucmass is ', NUC_MASS
      NU = E-EP
      QSQ = 4*E*EP*(SIN(THR/2))**2
                                !write (6,*)' qsq is ', QSQ
!      write (*,*) 'eps, e,ep,nu,a, th', EPS,E,EP,nu,A,thr

!      W = NU+A*NUC_MASS-EPS
!      WP = W**2+((A-1)*NUC_MASS)**2-NUC_MASS**2
                                !write (15,*) 'w, wp, qsq', w,wp,qsq
!      EPS=0.0
!    AG = 4*W**2-4*(NU**2+QSQ)
!     BG = SQRT(NU**2+QSQ)*(4*WP-4*(NU**2+QSQ))
!      PART1 = 4*W**2*(((A-1)*NUC_MASS)**2)
!      PART2 = 2*WP*(QSQ+NU**2)
!      PART3 = (QSQ+NU**2)**2
!      CG = PART1+PART2-PART3-WP**2
!      CG = -WP**2+4*(((A-1)*NUC_MASS)**2)*(A*NUC_MASS+NU)**2
!      rad = BG**2-4*AG*CG
!      write (*,*) 'A, B, C, RAD', AG,BG,CG,RAD
!      if (rad.lt.0) return
!      YSCALE = (-BG+SQRT(BG**2-4*AG*CG))/(2*AG)
!      backup =  (-BG-SQRT(BG**2-4*AG*CG))/(2*AG)


      qv2=qsq+nu**2
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
         backup =  (-BG-SQRT(BG**2-4*AG*CG))/(2*AG)
      endif


!      write(6,*)'got y of', YSCALE
                                !write (6,*) 'got backup y of ', backup
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
      RETURN
      END

       SUBROUTINE CALC_EP_FROM_Y(E,Y,THR,A,EPS,EP)
      IMPLICIT NONE
      INCLUDE 'constants_dble.inc'
      real*8 E,Y,THR,A,EPS,EP,AVAR,BVAR,CVAR,MT,MD,BETA,X,EP2
      MT = A*NUC_MASS
      MD = (A-1)*NUC_MASS
      X = E+MT-EPS-SQRT(MD**2+Y**2)
      BETA = 0.5*(X**2-NUC_MASS**2-Y**2-E**2)
      AVAR = (-E*COS(THR)+X)**2-Y**2
      BVAR = 2*(Y**2*E*COS(THR)-BETA*(-E*COS(THR)+X))
      CVAR = BETA**2-E**2*Y**2
      if (Y.LE.0)THEN
         EP = (-BVAR+SQRT(BVAR**2-4*AVAR*CVAR))/(2*AVAR)
      ELSE
         EP = (-BVAR-SQRT(BVAR**2-4*AVAR*CVAR))/(2*AVAR)
      ENDIF
                                !write (6,*) 'the other root is ', ep2
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

      subroutine y_w (e1, e2, thr,y_ww)
!     from paper by Gurvitz

      implicit none
      include 'constants_dble.inc'
      real*8 e1, e2, thr, y_ww, q3v, qsq, nu

      qsq=4*e1*e2*(sin(thr/2))**2
      nu=e1-e2
      q3v=sqrt(qsq+nu**2)

      y_ww=0.5*(-q3v+nu*sqrt(1+4*nuc_mass**2/qsq))
      return
      end

c-------------------------------------------------------------------------------------------------

	real*8 function emc_func_xem_aji(x,A) ! now compute the emc effect from our own fits.
        real*8 x,a
	real*8 emc

	if(A.eq.3) then

           if (x.lt.1.0) then

	      emc = 3.60574
	1	   -17.9930*x
	2	   +38.5132*x**2
	3	   -7.30356*x**3
	4	   -71.4791*x**4
	5	   +83.6995*x**5
	6	   -27.5481*x**6

	   endif

	else if(A.eq.4) then
	   if (x.lt.1.0) then

	      emc = 4.0905
	1	   -20.5864*x
	2	   +41.1811*x**2
	3	   -4.62533*x**3
	4	   -74.9766*x**4
	5	   +79.3210*x**5
	6	   -22.8134*x**6


	   endif


	else if(A.eq.9) then
	   if (x.lt.1.0) then

	      emc = 4.39669
	1	   -21.8919*x
	2	   +41.8461*x**2
	3	   -2.69928*x**3
	4	   -75.3117*x**4
	5	   +75.1341*x**5
	6	   -19.8517*x**6


	   endif


	else if(A.eq.12) then
	   if (x.lt.1.0) then

	      emc = 4.18521
	1	   -20.8696*x
	2	   +40.8226*x**2
	3	   -3.07714*x**3
	4	   -74.6169*x**4
	5	   +74.8572*x**5
	6	   -19.5719*x**6

	   endif



	else if((A.eq.64).or.(A.eq.56)) then
	   if (x.lt.1.0) then

	      emc = 4.01755
	1	   -19.8687*x
	2	   +39.0684*x**2
	3	   -4.74309*x**3
	4	   -71.9124*x**4
	5	   +78.8172*x**5
	6	   -23.8254*x**6


	   endif


	else if(A.eq.197) then

	   if (x.lt.1.0) then

	      emc = 4.20323
	1	   -21.0955*x
	2	   +40.3644*x**2
	3	   -2.65750*x**3
	4	   -73.9456*x**4
	5	   +75.4186*x**5
	6	   -20.7732*x**6



	   endif
	else
        write(*,*) '** in sigmodel_calc, emc fit not valid', A
	   stop
	endif
!        write(*,*) 'returning from emc with ', x, emc
	emc_func_aji= emc
	return
	end
c-----------------------------------------------------------------------------------------------
c for now  only do he4,be,c,cu

	subroutine highx_cor(anuc,x,cor)

	real*8 x,cor,anuc

         if(anuc.eq.3) then
	   if ((x.gt.0.9).and.(x.lt.1.4)) then
	    cor= -0.908273 + (4.13702*x) -(2.11462*x**2)
	   elseif (x.ge.1.4) then
	      cor=0.74
	   endif

        elseif(anuc.eq.4) then
	   if ((x.gt.0.9).and.(x.lt.1.17)) then
	      cor= 3.24553 - (3.47244*x) +  (1.11309*x**2)
	   elseif (x.ge.1.17) then
	      cor=0.7
	   endif
	elseif(anuc.eq.9) then

	   if ((x.gt.0.9).and.(x.lt.1.26)) then
	      cor= 0.779378 + (1.84808*x) - (1.7588*x**2)
	   elseif (x.ge.1.26) then
	      cor=0.3
	   endif

	elseif(anuc.eq.12) then
	   if ((x.gt.0.9).and.(x.lt.1.26)) then
	      cor=  1.09301 + (0.798708*x) - (0.939027*x**2)
	   elseif (x.ge.1.26) then
	      cor=0.55
	   endif

        elseif((anuc.eq.64).or.(anuc.eq.56)) then
	   if ((x.gt.0.9).and.(x.lt.1.3)) then
	      cor=  3.3275 - (3.94771*x) + (1.496*x**2)
	   elseif (x.ge.1.3) then
	      cor=0.68
	   endif

          elseif(anuc.eq.197) then
	   if ((x.gt.0.9).and.(x.lt.1.3)) then
           cor= -0.389135+ (3.08065*x)- (1.66297*x**2)
	   elseif (x.ge.1.3) then
	      cor=0.8
	   endif

	else
	   cor=1.
	endif

	return
	end

c-----------------------------------------------------------------------------------------------

      subroutine dis_highx_cor_argonne(anuc,x,cor)
      implicit none
      real*8 x,cor,anuc, frac,xlow1,xhigh1,xlow2,xhigh2

      xlow1=0.9
      xhigh1=0.95

      xlow2=1.3
      xhigh2=1.4

      frac=1.
      cor=1.
      if(anuc.eq.3) then
        cor=-2.12112*x+3.03449
        if((x.ge.xlow1).and.(x.le.xhigh1)) then
          frac = (x-xlow1)/(xhigh1-xlow1)
        endif
        cor=frac*cor+1.-frac


      elseif(anuc.eq.4) then
        cor=-1.76466*x+2.68897
        if((x.ge.xlow1).and.(x.le.xhigh1)) then
          frac = (x-xlow1)/(xhigh1-xlow1)
        endif
        cor=frac*cor+1.-frac

      elseif(anuc.eq.9) then
        cor=-1.8383*x+2.77253
        if((x.ge.xlow1).and.(x.le.xhigh1)) then
          frac = (x-xlow1)/(xhigh1-xlow1)
        endif
        cor=frac*cor+1.-frac

      elseif(anuc.eq.12) then
        cor=-1.32193*x+2.28754
        if((x.ge.xlow1).and.(x.le.xhigh1)) then
          frac = (x-xlow1)/(xhigh1-xlow1)
        endif
!     if((x.ge.xlow2).and.(x.le.xhigh2)) then
!     frac = (x-xlow2)/(xhigh2-xlow2)
!     frac=1.-frac
!     endif
!     if(x.gt.xhigh2) frac=0.
        cor=frac*cor+1.-frac
!     write(21,*) 'cor is ', cor, x
      elseif((anuc.eq.64).or.(anuc.eq.56)) then
!     cor=-2.21331*x+3.02106
        cor=1.
        cor=-1.46912*x+2.31581
!        if((x.ge.xlow1).and.(x.le.xhigh1)) then
!          frac = (x-xlow1)/(xhigh1-xlow1)
!        endif
!        cor=frac*cor+1.-frac
!     cor=1.0
      elseif(anuc.eq.197) then
        cor= -1.72192*x+2.65671
        if((x.ge.xlow1).and.(x.le.xhigh1)) then
          frac = (x-xlow1)/(xhigh1-xlow1)
        endif
        cor=frac*cor+1.-frac

      else
        cor=1.
      endif

      if(cor.lt.0.4) cor=0.4
 !     write (*,*) 'returning ', cor

      return
      end
c------------------------------------------------------------------------ !
C     --------Idea's same as aji's, but I did the fit to my x-binned
C     data
!  fort.targ_delta9_allresolved[_andinx]
c-----------------------------------------------------------------------------------------------
	subroutine global_cor(anuc,x,cor)
        implicit none
	real*8 x,cor,anuc,x1,x2,x3,x4,frac, poly_fit
        real*8 cor_x_lt_x1, cor_x_gt_x4, cor_polynom

	if(anuc.eq.2) then

           x1=0.375
	   x2=0.4
	   x3=1.9
	   x4=2.0

	elseif(anuc.eq.3) then

           x1=0.775
	   x2=0.8
	   x3=1.475
	   x4=1.5

	elseif(anuc.eq.4) then

           x1=0.775
	   x2=0.8
	   x3=1.475
	   x4=1.5

	elseif(anuc.eq.9) then

          x1=0.775
          x2=0.8
          x3=2.475
          x4=2.5

	elseif(anuc.eq.12) then

          x1=0.775
          x2=0.8
          x3=1.475
          x4=1.5

	elseif((anuc.eq.64).or.(anuc.eq.56)) then

          x1=0.775
          x2=0.8
          x3=1.975
          x4=2.0

	elseif(anuc.eq.197) then
          x1=0.775
          x2=0.8
          x3=1.475
          x4=1.5

 	else

	   write(*,*) '** in global_cor, unknown target', anuc
	   stop

	endif


	if (x.lt.x1) then
	   cor=1.0
	elseif ((x.ge.x1).and.(x.lt.x2)) then
	   frac = (x-x1)/(x2-x1)
	   cor= poly_fit(anuc,x)*frac + (1.-frac)
	elseif((x.ge.x2).and.(x.lt.x3)) then
	   cor= poly_fit(anuc,x)
	elseif ((x.ge.x3).and.(x.lt.x4)) then
	   frac = (x-x3)/(x4-x3)
	   cor= poly_fit(anuc,x4) *frac + poly_fit(anuc,x) *(1.-frac)
	elseif(x.ge.x4) then
	   cor=poly_fit(anuc,x4)
	endif


	return
	end


      real*8 function poly_fit(anuc,x)
      implicit none
      real*8 x, anuc, p0, p1, p2, p3, p4


	if(anuc.eq.2) then
          p4 = 1.47784948218043
          p3 = -5.23655714293828
          p2 = 6.21510111110596
          p1 = -3.02812744454624
          p0 = 1.52711763635932

	elseif(anuc.eq.3) then
          p4 = 1.02772560984094
          p3 = -6.08172162949192
          p2 = 12.7964815600646
          p1 = -11.5578850603455
          p0 = 4.76700715451466

	elseif(anuc.eq.4) then
          p4 = 1.59668103932533
          p3 = -8.95906115175537
          p2 = 18.1757099596398
          p1 = -15.9580650263325
          p0 = 6.10452712930702

	elseif(anuc.eq.9) then
          p4 = 0.336881628279981
          p3 = -2.47699274917202
          p2 = 6.42234681736013
          p1 = -7.03843528805277
          p0 = 3.74651945659282

	elseif(anuc.eq.12) then
          p4 = 0.774457758375471
          p3 = -4.77263675262629
          p2 = 10.6023709231639
          p1 = -10.1237696277791
          p0 = 4.49996964307447

	elseif((anuc.eq.64).or.(anuc.eq.56)) then
          p4 = 0.947130573149562
          p3 = -5.56050176808789
          p2 = 11.8522362056345
          p1 = -11.0423937551899
          p0 = 4.75470423802242

	elseif(anuc.eq.197) then
          p4 = 10.1348736361699
          p3 = -49.4362625708223
          p2 = 88.9917791787406
          p1 = -70.2751939566697
          p0 = 21.4898267656425

 	else

	   write(*,*) '** in poly_fit, unknown target'
	   stop

	endif

        poly_fit=p4*x**4+p3*x**3+p2*x**2+p1*x+p0




      end

!-----------------------------------------------------------------------------------------------------

	subroutine dis_highx_cor(anuc,x,cor)
        implicit none
        integer anuc
	real*8 x,cor, frac,xlow1,xhigh1,xlow2,xhigh2

	xlow1=0.9
	xhigh1=0.95

	xlow2=1.9
	xhigh2=2.0

	frac=1.
	cor=1.
	if(anuc.eq.2) then
	  cor=1.
          cor=-3.30482*x+ 4.10442
          if((x.ge.xlow1).and.(x.le.xhigh1)) then
            frac = (x-xlow1)/(xhigh1-xlow1)
          endif
          cor=frac*cor+1.-frac
!     xlow1=0.94
!	  xhigh1=.97
!     if (x.gt.xlow1) then
!     if (x.lt.1.15) then
!     cor= -3.1567509*x+ 4.04048547544833
!     else
!     cor=0.5
!     endif
!     if((x.ge.xlow1).and.(x.le.xhigh1)) then
!     frac = (x-xlow1)/(xhigh1-xlow1)
!     endif
!     endif
!     cor=frac*cor+1.-frac
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
!	   if((x.ge.xlow2).and.(x.le.xhigh2)) then
!	     frac = (x-xlow2)/(xhigh2-xlow2)
!	     frac=1.-frac
!	   endif
!	   if(x.gt.xhigh2) frac=0.
	   cor=frac*cor+1.-frac
!	   write(21,*) 'cor is ', cor, x
	 elseif((anuc.eq.63).or.(anuc.eq.64)) then
	   cor=-1.65829142599487*x+ 2.48174872208596
	   if((x.ge.xlow1).and.(x.le.xhigh1)) then
	     frac = (x-xlow1)/(xhigh1-xlow1)
	   endif
	   cor=frac*cor+1.-frac
!	   cor=1.0
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
	real*8 function emc_func_xem_prethesis(x,A) ! now compute the emc effect from our own fits.
	implicit none
        real*8 x,a,xtmp
	real*8 emc


c	if (x.le.1.0) then
	if(x.le.0.9) then
	   xtmp = x
	else
	   xtmp = 0.9
	endif

	   if(A.eq.2) then
	      emc =1.0
C it 2
	      if(xtmp.lt.0.2) then
		 emc=1.06
	      else
		 emc = 0.79515 +1.9777*xtmp - 3.9724*xtmp**2 -0.66967*xtmp**3
	1	   +8.3082*xtmp**4 - 5.5263*xtmp**5
	      endif
CHe3***********************************************************************************
	   else if(A.eq.3) then
C it 2
	      emc = 1.0118 +1.1029*xtmp -2.5081*xtmp**2 - 0.22033*xtmp**3
	1	   + 4.8120*xtmp**4 - 3.2865*xtmp**5
CHe4***********************************************************************************
	   else if(A.eq.4) then
C it2
	      emc = 0.84622 + 2.2462*xtmp - 4.7909*xtmp**2
	1	   + 0.065713*xtmp**3 + 7.6154*xtmp**4 - 5.2029*xtmp**5

	   else if(A.eq.9) then
C it 2
	      emc = 0.80887 + 3.9354*xtmp - 8.6056*xtmp**2 -0.16342*xtmp**3
	1	   + 14.074*xtmp**4 -9.3065*xtmp**5
	   else if(A.eq.12) then
C it 2
	      emc = 0.8279 + 3.5070*xtmp -7.5807*xtmp**2
	1	   -0.60935*xtmp**3 +13.081*xtmp**4 -8.5083*xtmp**5
	   else if(a.eq.27) then
	      emc = 0.98645 + 3.0385*xtmp - 22.072*xtmp**2 + 74.981*xtmp**3
	1	   - 132.97*xtmp**4 + 113.06*xtmp**5 -35.612*xtmp**6
	   else if((A.eq.64).or.(a.eq.63).or.(a.eq.56)) then
C it 2
	      emc = 1.1075 + 2.7709*xtmp - 6.5395*xtmp**2 -0.46848 *xtmp**3
	1	   +10.534*xtmp**4 - 6.6257*xtmp**5

	   else if(A.eq.197) then
C it 2
	      emc = 1.1404 + 4.0660*xtmp -10.318*xtmp**2 -1.9036*xtmp**3
	1	   + 21.969*xtmp**4 - 14.461*xtmp**5


	   else
	      write(*,*) '** in emc_func_xem, unknown target', a
	      stop
	   endif


	emc_func_xem_prethesis= emc
	return
	end


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
         aa = 1.72816025139459
        bb = 2.53114253906281
        cc = -2.72505067059703
        dd = -1.58637090989274
        ee = -16.3972900673533

      elseif(a.eq.3) then

        bb = 0.8
        cc = 0.06864880408328
        dd = -0.320972192919132
        ee = 0
        aa = 0.552199789237622

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
      endif

      tail_cor=(aa*exp(bb*x)+cc*x**6+dd*x**4+ee*x**2+ff)
      return
      end

      real*8 FUNCTION SIG_EL(E,EP,THR,A,Z,rosen_p,rosen_n,mott)
      IMPLICIT NONE
      INCLUDE 'constants.inc'
      real*8 E, EP, THR,MOTT,ELASTIC,GEP,GMP,TAU,W1,W2,QSQ,A,Z,N
      real*8 ROSEN_N,ROSEN_P,CS_MOTT,factor,factor2,factor3
      QSQ = 4*E*EP*(SIN(THR/2))**2
      N = A-Z
                                !write(6,*) 'in sigel, Z N A', Z,N,A
      MOTT=CS_MOTT(THR,E)
                                !factor = 0.524275
      factor = sqrt(qsq + nuc_mass**2)/nuc_mass ! epf/m
      factor3=(1+qsq*(tan(thr/2))**2/(2*nuc_mass**2))
      ELASTIC= MOTT!/factor
!      write (6,*) 'My elastic is ', ELASTIC
      factor2 = (1+2*E*(sin(THR/2))**2/nuc_mass)
!      write(*,*) 'factor 2 is ', factor2
      CALL ROSEN(ELASTIC,QSQ,THR,ROSEN_N,ROSEN_P)
      SIG_EL = N*ROSEN_N+Z*ROSEN_P
                                !write (6,*) 'factor and factor 2 and factor3', factor,factor2,factor3
!      write (6,*) 'sig el is ', sig_el
      RETURN
      END
      SUBROUTINE ROSEN(CS_NS,QSQ,THR,CSN,CSP)
      IMPLICIT NONE
      real*8 QSQ,GEP,GMP,GEN,GMN,KP,TAU,H(0:5),gep1,gmp1,gen1,gmn1
      real*8 P,MP,PRODUCT,CS_NS,CSN,CSP,THR
      REAL*4 Q2,one,two
      real*8 GEP4,GMP4,GEN4,GMN4
      INTEGER I,J,ig
      H(0) = 1.0007
      H(1) = 1.01807
      H(2) = 1.05584
      H(3) = 0.836380
      H(4) = 0.6864584
      H(5) = 0.672830
      MP = 0.938256
      KP = 1.7927
      P=0
                                !write (6,*) ' in rosen, cs_ns is ', cs_ns
      DO I=0,5
         PRODUCT = 1
         DO J=0,5
            IF (I.NE.J)THEN
               PRODUCT = PRODUCT*(SQRT(QSQ)-J)/(I-J)
            ENDIF
         END DO
         P = P+H(I)*PRODUCT
      END DO
                                !write (6,*) ' got P of ', P
                                !P = 1.0
      GEP1 = P/(1+QSQ/0.71)**2
      GMP1 = GEP1*(1.0+KP)
      TAU = QSQ/(4*MP**2)
      GEN1 =  TAU*(1.91)*GEP1/(1+5.6*TAU)
      GMN1 = GEP1*(-1.91)
                                !write (6,*)'IN ROSEN, GEP GMP  GEN GMN', GEP,GMP,GEN,GMN
!      write(15,*) thr*180/3.14159, qsq, gep1, gmp1, gen1,gmn1

      Q2 = (QSQ)/(.1973289**2)
                                !CALL BBG_FF(real*8(QSQ),GEP4,GMP4,GEN4,GMN4)
                                !GMP4 = GMP4*2.79
                                !GMN4 = GMN4*(-1.91)
                                !write (6,*)'BBG_FF, GEP GMP  GEN GMN', GEP4,GMP4,GEN4,GMN4
      one = 0.0
      two = 0.0
      ig = 13
!      write(*,*) 'sending nform', q2, ig,qsq
      CALL nform(dble(13.0),QSQ,GEP4,GEN4,GMP4,GMN4)
!      write(16,*) thr*180/3.14159, qsq, gep4, gmp4, gen4,gmn4
      GEP = GEP4
      GMP = GMP4
      GEN = GEN4
      GMN = GMN4
 !     GEP = GEP1
 !     GMP = GMP1
 !     GEN = GEN1
 !     GMN = GMN1
!      write (6,*)'Nform, GEP GMP  GEN GMN', GEP,GMP,GEN,GMN
!      write (6,*)'Nform4, GEP GMP  GEN GMN', GEP4,GMP4,GEN4,GMN4
      CSN= CS_NS*(((GEN**2+TAU*GMN**2)/(1+TAU))+2*tau*GMN**2*(tan(THR/2))**2)
      CSP=CS_NS*(((GEP**2+TAU*GMP**2)/(1+TAU))+2*TAU*GMP**2*(tan(THR/2))**2)
!      write (6,*)' SIGP, SIGN', CSP,CSN
      END

      real*8 function get_alpha_tn(e,ep,thr)

      include 'constants.inc'
      real*8 e, ep, thr


      q2=4*e*ep*(sin(thr/2))**2
      nu=e-ep
      x=q2/2/nuc_mass/nu
      q3v=sqrt(q2+q2**2/4/nuc_mass**2/x**2)
      q_minus=nu-q3v
      inv_mass=sqrt(4*nuc_mass**2+4*nu*nuc_mass-q2)
      get_alpha_tn=2-(q_minus+2*nuc_mass)/(2*nuc_mass)*(1+sqrt(inv_mass
     +  **2-4*nuc_mass**2)/inv_mass)



      end

  !this is the postthesis iteration of daves
	real*8 function emc_func_xem(x,A) ! now compute the emc effect from our own fits.
	implicit none
        integer a
        real*8 x,xtmp
	real*8 emc


c	if (x.le.1.0) then
	if(x.le.0.9) then
	   xtmp = x
	else
	   xtmp = 0.9
	endif

	emc =1.0
CDeuterium******************************************************************************
	if(A.eq.2) then
C it 2
c	   if(xtmp.lt.0.2) then
c	      emc=1.06
c	   else
c	      emc = 0.79515 +1.9777*xtmp - 3.9724*xtmp**2 -0.66967*xtmp**3
c	1	   +8.3082*xtmp**4 - 5.5263*xtmp**5
c	   endif
c	   emc = emc*0.96689
c it 3
	   emc = 0.70443 +2.3742*xtmp - 4.6566*xtmp**2 -0.78540*xtmp**3
	1	+9.3838*xtmp**4 - 6.1256*xtmp**5
CHe3***********************************************************************************
	else if(A.eq.3) then
C it 2
c	      emc = 1.0118 +1.1029*xtmp -2.5081*xtmp**2 - 0.22033*xtmp**3
c	1	   + 4.8120*xtmp**4 - 3.2865*xtmp**5
C it 3
	   emc = 0.92170 +1.7544*xtmp -3.7324*xtmp**2 - 0.24293*xtmp**3
	1	+ 6.7613*xtmp**4 - 4.6089*xtmp**5
CHe4***********************************************************************************
	else if(A.eq.4) then
C it2
c            emc = 0.84622 + 2.2462*xtmp - 4.7909*xtmp**2
c	1	   + 0.065713*xtmp**3 + 7.6154*xtmp**4 - 5.2029*xtmp**5
C it3
	   emc = 0.70050 + 3.1241*xtmp - 6.1738*xtmp**2
	1	- 0.049988*xtmp**3 + 9.3053*xtmp**4 - 6.1348*xtmp**5
C Be**********************************************************************************
	else if(A.eq.9) then
C it 2
c	      emc = 0.80887 + 3.9354*xtmp - 8.6056*xtmp**2 -0.16342*xtmp**3
c	1	   + 14.074*xtmp**4 -9.3065*xtmp**5
C it 3
	   emc = 0.46324 + 6.1220*xtmp - 12.184*xtmp**2 -1.0956*xtmp**3
	1	+ 20.316*xtmp**4 -12.899*xtmp**5
C Carbon**********************************************************************************
	else if(A.eq.12) then
C it 2
c         emc = 0.8279 + 3.5070*xtmp -7.5807*xtmp**2
c	1	   -0.60935*xtmp**3 +13.081*xtmp**4 -8.5083*xtmp**5
C it 3
	   emc = 0.63653 + 4.6458*xtmp -9.2994*xtmp**2
	1	-1.2226*xtmp**3 +16.157*xtmp**4 -10.236*xtmp**5
C Al**********************************************************************************
	else if(a.eq.27) then
	   emc = 0.98645 + 3.0385*xtmp - 22.072*xtmp**2 + 74.981*xtmp**3
	1	- 132.97*xtmp**4 + 113.06*xtmp**5 -35.612*xtmp**6
C Copper**********************************************************************************
	else if((A.eq.64).or.(A.eq.56)) then
C it 2
c	      emc = 1.1075 + 2.7709*xtmp - 6.5395*xtmp**2 -0.46848 *xtmp**3
c	1	   +10.534*xtmp**4 - 6.6257*xtmp**5
c it 3
	   emc = 0.58372 + 6.0358*xtmp - 11.988*xtmp**2 -1.0211*xtmp**3
	1	+18.567*xtmp**4 - 11.482*xtmp**5
C Gold**********************************************************************************
	else if(A.eq.197) then
C it 2
c	      emc = 1.1404 + 4.0660*xtmp -10.318*xtmp**2 -1.9036*xtmp**3
c	1	   + 21.969*xtmp**4 - 14.461*xtmp**5
C it 3
	   emc = 0.44132 + 8.1232*xtmp -16.141*xtmp**2 -5.6562*xtmp**3
	1	+ 35.606*xtmp**4 - 22.008*xtmp**5

	else
	   write(*,*) '** in emc_func_xem, unknown target',a
	   stop
	endif

	emc_func_xem= emc
	return
	end


      subroutine donal_rosen(e0,ep,th,ig,sigep,sigen,sigmot)
      implicit none
      real*8 e0, ep, th, ig, sigep, sigen, thr, rads, rmp,cs2
      real*8 s2, t2, frec, elos, qq, hbarc,fscnst, tau, qqf, gep2
      real*8 gmp2, gen2, gmn2, sigmot, c1, c2, c3, d1, d2, d3
      real*8 gen, gmn, gep, gmp

!      ig=dble(13.0)
!      write (*,*) 'got ', e0, ep, th, ig
      rmp =.93827
      rads = 3.14159/180.
      thr = th*rads/2.

      cs2 = (cos(thr))**2
      s2 = (sin(thr))**2
      t2 = (tan(thr))**2
      frec = 1. + 2.*e0*s2/rmp
!     ef = e0/frec
      elos = e0-ep
      qq = 4.*e0*ep*s2
      hbarc = 0.1973289
      fscnst = 1./137.02

      tau=qq/(4.*(rmp)**2)

      qqf =  qq/hbarc**2

      call nform(dble(13.0),qq,gep,gen,gmp,gmn)
      gep2 = gep**2

      gmp2 = gmp**2

      gen2 = gen**2

      gmn2 = gmn**2

      sigmot =    (fscnst/(2.*e0*s2))**2*cs2
      sigmot = sigmot * 1.0e4 * hbarc**2 !micro barns
      c1=gep2

      c2=gmp2*tau

      c3=c2*2.*(1.+tau)*t2

      sigep = sigmot*(c1+c2+c3)/(1.+tau)
c
      d1 = gen2

      d2 = gmn2*tau

      d3 = d2*2.*(1. + tau)*t2

      sigen = sigmot*(d1+d2+d3)/(1.+tau)

      return
      end

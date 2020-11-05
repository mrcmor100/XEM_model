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

      open(unit=11,file='./input/helium3.inp', status
     +  ='old')
      write(*,*) 'opened file'

 40   READ(11,*,IOSTAT=EOF) e,ep,theta,a,z
        if(eof.ge.0) then
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

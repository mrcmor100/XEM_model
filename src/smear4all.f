	subroutine smear4all(eic1,ep1,theta1,aa1,zz1,ann1,esep1,pmax1,
     &	innt1,innp1,f01,bigB1,ag1,bg1,alpha1,sigdeep)

	implicit none
	real*8 wp,ww1,ww2, w11, w22, nu,p
	real*8 eic1,ep1,theta1,aa1,zz1,ann1,esep1,pmax1
        real*8  innt1,innp1,ag1,bg1,cg1,dg1,sigdeep
        real*8 alpha1, f01,bigB1, alpha, bigB, f0

        real*8 eic, ep, theta, aa, zz, ann, esep, pmax
        real*8 innt, innp, bg, ag, cg, dg, agd, bgd, cgd
        real*8 dgd, tt2, sigm, w2a, w1a, fackf, weight1
        real*8 weight2, arg1, arg2, arg21, arg22
        real*8 pi, s2, anu, q2, amp, amn, x, amp2, ama_1
        real*8 ama, q3, w1ap, w2ap, w1an, w2an, du, dp, akf
        real*8 wp_max, w1ap1, w2ap1, w1an1, w2an1, u, ei
        real*8 w1, w2, inp, wpmin2, anup, radical

	real*8 rk,rho
	integer iu, ip
	real*8 rq2		!to pass to w1w2

	eic = eic1
	ep = ep1
	theta = theta1
	aa = aa1
	zz = zz1
	ann = ann1
	esep = esep1
	pmax = pmax1
	innt = innt1
	innp = innp1
	ag = ag1
	bg = bg1
	f0 = f01
	BigB = bigB1
	alpha=alpha1

	f0=f0*1000
	bigB=bigB*1000
	alpha=alpha/1000
	ag=ag*1000
	bg=bg*1000

	pi = 3.14159265
	s2 = sin(theta/2.)**2
	anu = eic - ep
	q2 = 4.*eic*ep*s2
	rq2=q2		! to pass to w1w2

	amp = 0.938273
	x = q2/2./amp/(eic-ep)
	amp2 = amp*amp
	amn = 0.939565
	if(aa.eq.2) then
	  ama_1=0.93827
	elseif(aa .gt. 2 .and. aa .lt. 4)then !he3
	  ama_1 = 1.87609
	elseif (aa.eq.4) then
	  ama_1 = 2.8094
	elseif (aa.eq.9) then
	  ama_1 = 7.5402
	elseif (aa.eq.10) then
	  ama_1 = 10.25553-2.0*0.938272
	elseif (aa.eq.11) then
	  ama_1 =  10.25553-0.938272
	elseif (aa.eq.12) then
	  ama_1 = 10.25553
	elseif (aa.eq.27) then
	  ama_1 = 24.205
	elseif (aa.eq.56) then
	  ama_1 = 51.1743
	elseif ((aa.eq.63).or.(aa.eq.64)) then
	  ama_1 = 57.6886
	elseif (aa.eq.197) then
	  ama_1 = 182.5394
	endif
	ama = ama_1 + amp - esep

	q3 = sqrt(q2 +anu*anu)
	w1ap = 0.0
	w2ap = 0.0
	w1an = 0.0
	w2an = 0.0
	du = 2./innt
	dp = pmax/innp

	wp_max = 0.0
!	do smearing correction
	do ip = 1,int(innp)
	  p = (ip - .5)*dp
	  w1ap1 = 0.0
	  w2ap1 = 0.0
	  w1an1 = 0.0
	  w2an1 = 0.0

	  do 10 iu = 1,int(innt)
	    u = (iu - .5)*du
	    u = u - 1.
	    ei = ama - sqrt(p*p+ama_1*ama_1)
	    anup = (ei*anu - p*u*q3)/amp
	    wpmin2 = (amp + .135)**2
	    if((ei*ei - p*p + 2.*ei*anu - 2.*p*u*q3 - q2).le
	1     .wpmin2) goto 10
	    radical =(ei*ei - p*p + 2.*ei*anu - 2.*p*u*q3 - q2)
	    wp = sqrt(radical)
	    w11=0.
	    w22=0.
	    nu = (wp**2 - amp**2 + q2) / 2. / amp
	    call F1F2IN21(dble(1.),dble(1.),dble(rq2),dble(wp**2),w11,w22)
	    if((w11.eq.0).or.(w22.eq.0)) goto 10
	    w11=w11/amp
	    w22=w22/nu
	    w1=w11
	    w2=w22
	    if (wp .gt. wp_max) wp_max = wp
	    arg1 = w1 + (1 -u*u)*p*p*w2/2./amp2
	    arg21 = (1. + p*u*q2/amp/anup/q3)**2*anup*anup/anu/anu
	    arg22 = (1.- u*u)*p*p*q2/2./amp2/q3/q3
	    arg2 = w2*(arg21 + arg22)
	    w1ap1 = w1ap1 + arg1*du
	    w2ap1 = w2ap1 + arg2*du
	    w11=0.
	    w22=0.
	    call F1F2IN21(dble(0.),dble(1.),dble(rq2),dble(wp**2),w11,w22)
	    if((w11.eq.0).or.(w22.eq.0)) then
c	      write(*,*) 'got w11 or w22 as zero: ', aa1, zz1, rq2, ip, iu 
	      goto 10
	    else
!	      write(*,*) 'got w11 and w22 n zero: ', aa1, zz1, rq2 
	    endif
	    w11=w11/amp
	    w22=w22/nu
	    w1=w11
	    w2=w22
	    arg1 = w1 + (1 -u*u)*p*p*w2/2./amp2
	    arg21 = (1. + p*u*q2/amp/anup/q3)**2*anup*anup/anu/anu
	    arg22 = (1.- u*u)*p*p*q2/2./amp2/q3/q3
	    arg2 = w2*(arg21 + arg22)
	    w1an1 = w1an1 + arg1*du
	    w2an1 = w2an1 + arg2*du
 10	  continue

	  weight1 = ((exp(-(ag*p)**2)*(f0-bigB)*alpha**2)*(((ag**2)
	1   *(alpha**2+p**2)+1)/(alpha**2+p**2)**2))
	  if (aa.eq.2) then
	    if(p.lt.0) then
	       write(*,*) 'No idea why p<0.'
	      weight1=weight1-(bg*bigB*exp(-bg*abs(p))/(2*p))
	    else
	      weight1=weight1+(bg*bigB*exp(-bg*abs(p))/(2*p))
	    endif
	  else
	      weight1=weight1+(bg**2*bigB*exp(-(bg*p)**2))
	  endif
	  weight1=weight1*p*p*dp/pi!*2.48176
	  w1ap = w1ap + w1ap1*weight1
	  w2ap = w2ap + w2ap1*weight1
	  w1an = w1an + w1an1*weight1
	  w2an = w2an + w2an1*weight1
	enddo

	fackf = 2.*pi
	w1a = fackf*(zz*w1ap + ann*w1an)
	w2a = fackf*(zz*w2ap + ann*w2an)
	sigm = cos(theta/2.)/(2.*137.*eic*s2)
!        sigm=0.389379*sigm**2  !(GeV**2 mb )
!sig_mott = hc_2*(alpha*cos(th/2))**2/(2*e1*sin(th/2)**2)**2
	sigm = (0.1973*sigm)**2
	tt2 = tan(theta/2.)**2
	sigdeep = 1.e-2*1.e9*sigm*(w2a + 2.*w1a*tt2) !fm2 to b, then b to nb!
	return
	end

      program run_xem_model
      implicit none

      integer eof
      character*80 rawname, filename, ofilename
      !Passed from input file to xem_model
      real*8 E0, EP, THETA, A, Z
      !Passed from xem_model to run_xem_model
      real*8 Y, X, dis_XS, qe_XS
      logical reload_params/.true./
      integer*4	last_char

      read(*,1968) rawname
 1968 format(a)
 2002 format (8(E13.5,1x))
CAM Open input file.
      filename = 'input/'//rawname(1:last_char(rawname))//'.inp'
      open(unit=11,status='old',file=filename)

CAM Open output file.
      filename = './output/'//rawname(1:last_char(rawname))//'.out'
      open (unit=14,status='unknown',file=filename)

 40   READ(11,*,IOSTAT=EOF) E0, EP, THETA, A, Z
      if(eof.ge.0) then
        call xem_model(E0, EP, THETA, A, Z, Y, X, dis_XS, qe_XS, 
     >        reload_params)

        !Redundant check on y. Don't add unphysical points to table
        if(Y.lt.1.E19.and.y.ne.0) then
           write(14,2002) y, a, z, THETA, EP, x, dis_XS, qe_XS
        endif
        
        goto 40
      endif !if file has stuff
      close(11);
      end


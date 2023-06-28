! gfortran -c  -O2 -fallow-argument-mismatch -std=legacy test_l3dlocloc.f
! gfortran -o int2 test_l3dlocloc.o laprouts3d.o yrecursion.o l3dtrans.o prini.o rotviarecur3.o
! gfortran -o int2 test_l3dlocloc.o 
      implicit None
      real *8 ztrg(3),source(3)
      real *8 c1(3),c2(3)
      real *8, allocatable :: dc(:,:)
      real *8, allocatable :: pots(:),flds(:,:),hesss(:,:)
      real *8, allocatable :: opots(:),oflds(:,:),ohesss(:,:)
      complex *16, allocatable :: locexp1(:,:),locexp2(:,:)
      real *8, allocatable :: charge(:),dipvec(:,:)
      real *8 wlege(100000),err_out(3)
      real *8 scarray_loc(100000)
      real *8 bsize, rscale1, rscale2, shift
      integer nd,idim,i,ns,nt,nlege,nterms1,nterms2,lw7,lused7,nn,ier
      real *8 eps, thresh

c     set scale for center 1 and center 2
      bsize = 1.0d0
      rscale1 = bsize/4
      rscale2 = rscale1

      allocate(charge(1),dipvec(3,1))

c     local exp center 1
      c1(1)=0.248d0
      c1(2)=0.251d0
      c1(3)=0.249d0
c     source position
      source(1)=c1(1)-0.47
      source(2)=c1(2)+0.51
      source(3)=c1(3)+0.49
			charge(1) = 1.1d0
			dipvec(1,1) = 0.31
			dipvec(2,1) = 0.33
			dipvec(3,1) = 0.43
			ns = 1

c     local exp center 2
      shift = 0.001
      c2(1) = c1(1) + shift
      c2(2) = c1(2) + shift*1.01
      c2(3) = c1(3) - shift*0.99

c     target position
      ztrg(1) = c1(1) + 0.01
      ztrg(2) = c1(2) + 0.02
      ztrg(3) = c1(3) - 0.01
      nt = 1

      allocate(opots(nt),pots(nt))
      allocate(oflds(3,nt),flds(3,nt))
      allocate(ohesss(6,nt),hesss(6,nt))


c
c     direct calculation:
c
      do i=1,nt
				opots(i) = 0
				oflds(1,i) = 0
				oflds(2,i) = 0
				oflds(3,i) = 0
				ohesss(1,i) = 0
				ohesss(2,i) = 0
				ohesss(3,i) = 0
				ohesss(4,i) = 0
				ohesss(5,i) = 0
				ohesss(6,i) = 0
      enddo
			
      thresh = 1.0d-15
			call lpotfld3dall_targ(0,source,charge,ns,ztrg,
     1    	opots,oflds)


c     local exp order for center 1
      nterms1 = 39
c     local exp order for center 2
      nterms2 = 7
			
			ier = 0
			allocate(locexp1(0:nterms1,-nterms1:nterms1))
			locexp1 = 0.0d0
			call l3dformta(ier,rscale1,source,charge,ns,c1,
     1		           nterms1,locexp1)
			
		  ! print *, "local expansion: ", locexp1(1,:)
		  
      call l3dtaeval(rscale1,c1,locexp1,nterms1,
     1		ztrg,pots,0,flds,ier)	
			
			print *,"local exp1 eval error:"
      call errprinth2(nt,pots,opots,flds,oflds,hesss,ohesss,
     1   err_out)

cc    shift local exp from center 1 to center 2
      allocate(locexp2(0:nterms2,-nterms2:nterms2))
			locexp2 = 0.0d0
			call l3dloclocquadu(rscale1,c1,locexp1,nterms1,
     1           rscale2,c2,locexp2,nterms2,ier)
		  
      do i=1,nt
				pots(i) = 0
				flds(1,i) = 0
				flds(2,i) = 0
				flds(3,i) = 0
				hesss(1,i) = 0
				hesss(2,i) = 0
				hesss(3,i) = 0
				hesss(4,i) = 0
				hesss(5,i) = 0
				hesss(6,i) = 0
      enddo

			call l3dtaeval(rscale2,c2,locexp2,nterms2,ztrg,
     1		pots,0,flds,ier)	
		  
		  print *,"local exp2 eval error:"
			call errprinth2(nt,pots,opots,flds,oflds,hesss,ohesss,
     1   err_out)

      end

c
      subroutine errprinth2(nt,pot,opot,fld,ofld,hess,ohess,errs)
      implicit real *8 (a-h,o-z)
      real *8 pot(nt),opot(nt),fld(3,nt)
      real *8 ofld(3,nt),hess(6,nt),ohess(6,nt)
      real *8 errs(3)
 1000  format(4D15.5) 
      
      errp = 0
      dddp = 0
      errf = 0
      dddf = 0
      errh = 0
      dddh = 0

      do i=1,nt
				errp = errp + abs(pot(i)-opot(i))**2
				dddp = dddp + abs(opot(i))**2
      enddo

c

      errs(1) = sqrt(errp/dddp)
			errs(2) = 0
			errs(3) = 0
      write(*,'(a,e11.4,a,e11.4,a,e11.4)') 
     1     'pot error=',errs(1),'   grad error=',errs(2),
     2     '   hess error=',errs(3)   
      return
      end
c			
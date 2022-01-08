program femimp

!-------------------------------------------------------------------------
!Two  -dimensional analysis of plane or axisymmetric strain analysis
!Three-dimensional analysis of an elastic solid
!-------------------------------------------------------------------------
use main
use geom
implicit none

integer,parameter::iwp=selected_real_kind(15)

!------------------------------------------------------------------------------
!1.0      Global Variables
!------------------------------------------------------------------------------

!------------------------------Scalars Integers---------------------------

 integer	::	fixed_freedoms		!Number of fixed displacements
 integer	::	i					!Simple counter
 integer	::	iel					!Simple counter
 integer	::	k					!Simple counter
 integer	::	loaded_nodes		!Number of loaded nodes
 integer	::	ndim=2				!Number of dimensions
 integer	::	ndof				!Number of degrees of freedom per element
 integer	::	nels				!Number of elements
 integer	::	neq					!Number of degrees of freedom in the mesh
 integer	::	nip					!Number of integration point per element
 integer	::	nim					!Number of integration point per element
 integer	::	nn					!Number of nodes in the mesh 
 integer	::	nod					!Number of nodes per element
 integer	::	nodof=2    			!Number of degrees of freedom per node
 integer	::	nprops=30			!Number of material properties
 integer	::	np_types			!Number of different property types
 integer	::	nr					!Number of restrained nodes
 integer	::	nst=3  				!Number of stress (strain) terms
 integer	::	nxe					!Number of elements in x(r)-direction
 integer	::	nye,nze				!Number of elements in y(z)-direction
 integer	::	nlen				!Number of letters
 integer	::	j

!------------------------------Others-------------------------------------

 integer	::	iy					!Number of elements in y(z)-direction
 integer	::	iters				!Number of letters
 integer	::	kgaus				!Number of letters
 integer	::	nstep
 integer	::	npri
 integer	::	nres
 integer	::	miter
 integer	::	geox = 0

!------------------------------Scalars Reals------------------------------

 real (iwp)	:: det					!Determinant of the jacobian matrix
 real (iwp)	:: one=1.0_iwp			
 real (iwp)	:: penalty=1.0e20_iwp	
 real (iwp)	:: zero=0.0_iwp		
 real (iwp)	:: facto=0.25_iwp		
 real (iwp)	:: toler=0.0001_iwp	
 real (iwp)	:: ratio				
 real (iwp)	:: resid				
 real (iwp)	:: retot				
 real (iwp)	:: oneh=100_iwp		
 real (iwp)	:: three=3.0_iwp		
 real (iwp)	:: pi=3.14159265358979_iwp	
 real (iwp)	:: twopi=6.283185307179586_iwp 
 real (iwp)	:: half=0.5_iwp			
 real (iwp)	:: two=2.0_iwp		
 real (iwp)	:: pt5=0.50_iwp		

 real (iwp)	:: sigm				
 real (iwp)	:: dsbar				
 real (iwp)	:: theta				
 real (iwp)	:: sy					
 real (iwp)	:: f					
 real (iwp)	:: xgash				
 real (iwp)	:: xgish				
 real (iwp)	:: s1,s2,s3			

 real (iwp)	:: area
 real (iwp)	:: a0,a1,a2,a3,a4,a5,a6,a7
 real (iwp)	:: dtim
 real (iwp)	:: delta
 real (iwp)	:: gaama
 real (iwp)	:: aalfa
 real (iwp)	:: beeta
 real (iwp)	:: time
 real (iwp)	:: hards
 real (iwp)	:: dvolu
 real (iwp)	:: extra
 real (iwp)	:: anu,PP,QQ,RR,sxx,syy,sxy
 real (iwp)	:: em
 real (iwp)	:: pv
 real (iwp)	:: rhow

!------------------------------Logical------------------------------------

 LOGICAL	  :: consistent=.TRUE.
 LOGICAL	  :: dynamic   =.FALSE.

!------------------------------Scalars Characters-------------------------

 character (len=15)		::	element			!Element type
 character (len=15)		::	dir				!Element and node numbering direction
 character (len=15)		::	type_2d			!Plane' or 'Axisymmetry'
 character (len=15)		::	argv			!File
 character (len=15)		::  generate		!Generate mesh
 character (len=45)		::  title			!Problem title

!------------------------------Dynamic Integer Arrays---------------------

 integer,allocatable	::	etype   (:)		!Element property type vector
 integer,allocatable	::	g       (:)		!Element steering vector
 integer,allocatable	::	g_g   (:,:)		!Global element steering matrix
 integer,allocatable	::  g_num (:,:)		!Global element node numbers matrix
 integer,allocatable	::	kdiag   (:)		!Diagonal term location vector
 integer,allocatable	::	nf    (:,:)		!Nodal freedom matrix
 integer,allocatable	::	no      (:)		!Fixed freedom numbers vector
 integer,allocatable	::	node    (:)		!Fixed nodes vector
 integer,allocatable	::	num     (:)		!Element node numbers vector
 integer,allocatable	::	sense   (:)		!Sense of freedoms to be fixed vector

!------------------------------Dynamic Real Arrays-------------------------

 real (iwp),allocatable	::	bee		(:,:)	!Strain-displacement matrix
 real (iwp),allocatable	::	coord	(:,:)	!Element nodal coordinates
 real (iwp),allocatable	::	dee		(:,:)	!Stress  strain matrix
 real (iwp),allocatable	::	der		(:,:)	!Shape function derivates local  coord.
 real (iwp),allocatable	::	deriv	(:,:)	!Shape function derivates global coord.
 real (iwp),allocatable	::	eld       (:)	!Element nodal displacements
 real (iwp),allocatable	::	fun       (:)	!Shape functions
 real (iwp),allocatable	::	gc        (:)	!Integrating point coordinates
 real (iwp),allocatable	::	g_coord (:,:)	!Global nodal coordinates
 real (iwp),allocatable	::	jac     (:,:)	!Jacobian matrix
 real (iwp),allocatable	::	km      (:,:)	!Element stiffness matrix
 real (iwp),allocatable	::	kv        (:)	!Global  stiffness matrix
 real (iwp),allocatable	::	loads     (:)	!Nodal loads and displacements
 real (iwp),allocatable	::	points  (:,:)	!Integrating point local coordinates
 real (iwp),allocatable	::	prop    (:,:)	!Element properties (E and v for each el.)
 real (iwp),allocatable	::	sigma     (:)	!Stress terms
 real (iwp),allocatable	::	value     (:)	!Fixed values of displacements
 real (iwp),allocatable	::	weights   (:)	!Weighting coefficients
 real (iwp),allocatable	::	x_coords  (:)	!x(r)-coordinates of mesh layout
 real (iwp),allocatable	::	y_coords  (:)	!y(r)-coordinates of mesh layout
 real (iwp),allocatable	::	z_coords  (:)	!y(r)-coordinates of mesh layout

 real (iwp),allocatable	::	force     (:)	!Nodal loads and displacements
 real (iwp),allocatable	::	eload     (:)	!Nodal loads and displacements
 real (iwp),allocatable	::	bload     (:)	!Nodal loads and displacements
 real (iwp),allocatable	::	tload     (:)	!Nodal loads and displacements
 real (iwp),allocatable	::	rload     (:)	!Nodal loads and displacements
 real (iwp),allocatable	::	asdis     (:)	!Nodal loads and displacements
 real (iwp),allocatable	::	tdisp     (:)	!Nodal loads and displacements
 real (iwp),allocatable	::	strsg   (:,:)	!Nodal loads and displacements
 real (iwp),allocatable	::	strag   (:,:)	!Nodal loads and displacements
 real (iwp),allocatable	::	stran     (:)	!Stress terms
 real (iwp),allocatable	::	desig     (:)	!Stress terms
 real (iwp),allocatable	::	epstn     (:)	!Stress terms
 real (iwp),allocatable	::	stres     (:)	!Stress terms
 real (iwp),allocatable	::	pl      (:,:)	!Stress terms
 real (iwp),allocatable	::	effst     (:)	!Stress terms

 real (iwp),allocatable	::	ks        (:)	!Global  stiffness matrix
 real (iwp),allocatable	::	mm      (:,:)
 real (iwp),allocatable	::	mv        (:)
 real (iwp),allocatable	::	mc        (:)

 real (iwp),allocatable	::	ecm     (:,:)
 real (iwp),allocatable	::	cv        (:)
 real (iwp),allocatable	::	u0        (:)
 real (iwp),allocatable	::	u10       (:)
 real (iwp),allocatable	::	u20       (:)
 real (iwp),allocatable	::	d         (:)
 real (iwp),allocatable	::	v         (:)
 real (iwp),allocatable	::	a         (:)
 real (iwp),allocatable	::	vc		  (:)
 real (iwp),allocatable	::	kd		  (:)
 real (iwp),allocatable	::	rt		  (:)
 real (iwp),allocatable	::	kp        (:)
 real (iwp),allocatable	::	ve        (:)

!------------------------------------------------------------------------------
!1.0      Input and Initialisation
!------------------------------------------------------------------------------

 call getname(argv,nlen)
 open(10,file=argv(1:nlen)//'.dat')
 open(11,file=argv(1:nlen)//'.res1')
 open(12,file=argv(1:nlen)//'.res2')

 read(10,*)title
 write(11,'(a,a)') "Problem: ",title
 read(10,*)type_2d,element,nod,nels,nn,nip,nim,nodof,nst,ndim,generate

 if(generate=='yes') then
    IF(ndim==2) THEN 
 	    read(10,*) nxe,nye
        call mesh_size(element,nod,nels,nn,nxe,nye)
    ELSEIF(ndim==3) THEN
	    read(10,*) nxe,nye,nze
	    call mesh_size(element,nod,nels,nn,nxe,nye,nze) 
    END IF
 end if
 ndof=nod*nodof
 if(type_2d=='axisymmetric')nst=4

 allocate (nf(nodof,nn))
 allocate (points(nip,ndim))
 allocate (dee(nst,nst))
 allocate (g_coord(ndim,nn))
 allocate (coord(nod,ndim))
 allocate (jac(ndim,ndim))
 allocate (weights(nip))
 allocate (num(nod))
 allocate (g_num(nod,nels))
 allocate (der(ndim,nod))
 allocate (deriv(ndim,nod))
 allocate (bee(nst,ndof))
 allocate (km(ndof,ndof))
 allocate (eld(ndof))
 allocate (sigma(nst))
 allocate (g(ndof))
 allocate (g_g(ndof,nels))
 allocate (mm(ndof,ndof))
 allocate (gc(ndim))
 allocate (fun(nod))
 allocate (etype(nels))
 allocate (x_coords(nxe+1))
 allocate (y_coords(nye+1))
 allocate (z_coords(nze+1))

 if(generate=='yes') then
    IF(ndim==2) read(10,*)x_coords,y_coords
	IF(ndim==3) read(10,*)x_coords,y_coords,z_coords
 else
    read(10,*)g_coord
    read(10,*)g_num
 end if

 nf=1
 read(10,*)nr,(k,nf(:,k),i=1,nr)
 call formnf(nf)
 neq=MAXVAL(nf)
 allocate(loads(0:neq),kdiag(neq))

 read(10,*)np_types
 allocate (prop(nprops,np_types))
 read(10,*)prop
 etype=1
 if(np_types>1)read(10,*)etype
 read(10,*) toler,nstep,miter,dtim,aalfa,beeta,gaama,delta,npri,nres


 loads=zero
 read(10,*)loaded_nodes,(k,loads(nf(:,k)),i=1,loaded_nodes)
 read(10,*)fixed_freedoms
 loads(0)=zero
 if(fixed_freedoms/=0)then
   allocate(node(fixed_freedoms),sense(fixed_freedoms),                   &
     value(fixed_freedoms),no(fixed_freedoms))
   read(10,*)(node(i),sense(i),value(i),i=1,fixed_freedoms)
   do i=1,fixed_freedoms
     no(i)=nf(sense(i),node(i))
   end do
   kv(kdiag(no))=kv(kdiag(no))+penalty
   loads(no)=kv(kdiag(no))*value
 end if

!----------------Loop Elements to Find Global Arrays Sizes---------------------

 kdiag=0
 elements_1: do iel=1,nels
   num=g_num(:,iel)
   call num_to_g(num,nf,g)
   g_g(:,iel)=g
   call fkdiag(kdiag,g)
 end do elements_1
 do i=2,neq;kdiag(i)=kdiag(i)+kdiag(i-1);end do
 
!----------------------------Print initial data--------------------------------

 write(11,'(a)') "---------------------------------------------------------------------------"
 write(11,'(a)') "                             Global Coordinates                            "
 write(11,'(a)') "---------------------------------------------------------------------------"
 do k=1,nn
   write(11,'(a,i5,a,<ndim>e12.4)') "Node",k,"        ",g_coord(:,k)
 end do 
 write(11,'(a)') "---------------------------------------------------------------------------"
 write(11,'(a)') "                             Global Node Numbers                           "
 write(11,'(a)') "---------------------------------------------------------------------------"
 do k=1,nels
   write(11,'(a,i5,a,<nod>i5)') "Element",k,"        ",g_num(:,k)
 end do 
 write(11,'(2(a,I5))')                                                    &
   " There are",neq," equations and the skyline storage is",kdiag(neq)
 write(11,'(a)') "---------------------------------------------------------------------------"
 write(11,'(a)') "                             Material Parameters                           "
 write(11,'(a)') "---------------------------------------------------------------------------"
 do i=1,np_types
    write(11,'(a,i5)') "Material Label     ",i
	write(11,'(a,1e12.4)') "Young Modulus  =",	prop(23,i)
	write(11,'(a,1e12.4)') "Poisson        =",	prop(24,i)
    write(11,'(a,1e12.4)') "SpecificWeight =",	prop(26,i)
    write(11,'(a,1e12.4)') "UniaxialStress =",	prop(27,i)    
    write(11,'(a,1e12.4)') "Hardening      =",	prop(28,i)
 end do 
 write(11,'(a)') "---------------------------------------------------------------------------"
 write(11,'(a)') "                          Time Stepping Parameters                         "
 write(11,'(a)') "---------------------------------------------------------------------------"

 write(11,'(a,1e12.4)') "Tolerance   =",	toler
 write(11,'(a,i5)')     "Steps       =",	nstep
 write(11,'(a,i5)')     "MaxIters    =",	miter
 write(11,'(a,1e12.4)') "StepTime    =",	dtim
 write(11,'(a,1e12.4)') "Alfa        =",	aalfa
 write(11,'(a,1e12.4)') "Beta        =",	beeta
 write(11,'(a,1e12.4)') "Gamma       =",	gaama
 write(11,'(a,1e12.4)') "Delta       =",	delta
 write(11,'(a,i5)')     "StepPrint   =",	npri
 write(11,'(a,i5)')     "NodeToPrint =",	nres

!-----------------------------Allocate Real Arrays-----------------------------

 allocate(ks(kdiag(neq)))
 allocate(force(0:neq))
 allocate(eload(0:neq))
 allocate(tload(0:neq))
 allocate(tdisp(0:neq))
 allocate(asdis(0:neq))
 allocate(bload(ndof))
 allocate(strsg(nst,nip*nels))
 allocate(strag(nst,nip*nels))
 allocate(stran(nst))
 allocate(desig(nst))
 allocate(epstn(nip*nels))
 allocate(effst(nip*nels))
 
 allocate(stres(nst))
 allocate(pl(nst,nst))

 allocate(kv(kdiag(neq)))
 allocate(mv(kdiag(neq)))
 allocate(cv(kdiag(neq)))
 allocate(ecm(ndof,ndof))
 allocate(kp(kdiag(neq)))
 allocate(d(0:neq))
 allocate(v(0:neq))
 allocate(a(0:neq))
 allocate(u0(0:neq))
 allocate(u10(0:neq))
 allocate(u20(0:neq))
 allocate(vc(0:neq))
 allocate(kd(0:neq))
 allocate(rt(0:neq))
 allocate(ve(0:neq))
 allocate(mc(kdiag(neq)))

 tload=zero
 eload=zero
 asdis=zero
 tdisp=zero
 strsg=zero
 strag=zero
 epstn=zero
 vc=zero
 kd=zero
 cv=zero

 a0=dtim*dtim*(0.5_iwp-delta)
 a1=dtim*(one-gaama)
 a2=dtim*dtim*delta
 a3=dtim*gaama
 a4=one/a2
 a5=beeta*gaama*dtim
 a6=aalfa*gaama*dtim
 a7=one+a6


!------------------------------------------------------------------------------
!2.0      Mass matrix assembly with another number of integration points
!------------------------------------------------------------------------------
 DEALLOCATE (points,weights)
 ALLOCATE   (points(nim,ndim),weights(nim))
 CALL sample(element,points,weights)

 mv = zero
 gc = one
 elements_2: DO iel=1,nels
	num=g_num(:,iel)
	coord = TRANSPOSE(g_coord(:,num))
	g     = g_g(:,iel)
	mm    = zero
	area  = zero
	rhow  = prop(26,etype(iel))
	gauss_pts_2: DO i=1,nim
		CALL shape_fun(fun,points,i)
		CALL shape_der(der,points,i)
		jac   = MATMUL(der,coord)
		det   = determinant(jac)
		dvolu = det*weights(i)
		IF(type_2d=='axisymmetric')THEN
			gc    = MATMUL(fun,coord)
			dvolu = dvolu*gc(1)*twopi
		END IF
		area   = area + det*weights(i)
		IF(consistent)THEN
			CALL ecmat(ecm,fun,ndof,nodof)
			mm = mm + ecm*rhow*dvolu
		END IF
   END DO gauss_pts_2
   IF(.NOT.consistent) CALL elmat(area,rhow,mm)
   CALL fsparv(mv,mm,g,kdiag)
 END DO elements_2 
 mc = mv

!------------------------------------------------------------------------------
!3.0      Initial stiffness assembly for getting initial conditions
!------------------------------------------------------------------------------
 DEALLOCATE (points,weights)
 ALLOCATE   (points(nip,ndim),weights(nip))
 CALL sample(element,points,weights)

 kv    = zero
 kgaus = 0
 gc    = one
 elements_3: DO iel=1,nels
    sy = prop(27,etype(iel)) 
    em = prop(23,etype(iel))
	pv = prop(24,etype(iel))
	CALL deemat(dee,em,pv,type_2d)
	num   = g_num(:,iel)
	coord = TRANSPOSE(g_coord(:,num))
	g=g_g(:,iel)
	km    = zero
	gauss_pts_3: DO i=1,nip
	    kgaus = kgaus + 1
        effst(kgaus) = sy
		CALL shape_fun(fun,points,i)
		CALL shape_der(der,points,i)
		jac   = matmul(der,coord)
		det   = determinant(jac)
		CALL invert(jac)
		deriv = MATMUL(jac,der)
		gc    = MATMUL(fun,coord)
		call beemat(bee,deriv,type_2d,gc,fun)
		dvolu = det*weights(i)
		if(type_2d=='axisymmetric') dvolu=dvolu*gc(1)*twopi
		km = km + MATMUL(MATMUL(TRANSPOSE(bee),dee),bee)*dvolu
   END DO gauss_pts_3
   CALL fsparv(kv,km,g,kdiag)
 END DO elements_3 

!------------------------------------------------------------------------------
!4.0      Initial Condition: d,v,a (displacement, velocity, aceleration)
!------------------------------------------------------------------------------

 d  = zero
 v  = zero
 a  = zero
 cv = aalfa*mv + beeta*kv
 CALL sparin(mc,kdiag)
 CALL linmul_sky(cv,v,vc,kdiag)
 CALL linmul_sky(kv,d,kd,kdiag)
 rt = zero
 rt = loads*load(zero)
 a  = rt - vc - kd
 CALL spabac(mc,a,kdiag) 
 d(0) = zero
 a(0) = zero
 v(0) = zero
 
 strsg = zero
 eload = zero
 strag = zero


!------------------------------------------------------------------------------
!5.0      Time Stepping Loop
!------------------------------------------------------------------------------
 write(11,'(a)') "---------------------------------------------------------------------------"
 write(11,'(a)') "                          Time Stepping Process                            "
 write(11,'(a)') "---------------------------------------------------------------------------"

 time = zero
 u0   = zero    
 timesteps: DO j=1,nstep
    
    WRITE(11,'(a,i5)') "TimeStep: ",j
	time = time + dtim

!------------------------------------------------------------------------------
!6.0      Stiffness matrix assembly
!------------------------------------------------------------------------------
	!IF(j==1) THEN
	kv    = zero
	cv    = zero
	kgaus = 0
	gc    = one
	elements_4: DO iel=1,nels
		num   = g_num(:,iel)
		coord = TRANSPOSE(g_coord(:,num))
		g     = g_g(:,iel)
		eld   = u0(g)
		km    = zero
		gauss_pts_4: DO i=1,nip
		    kgaus = kgaus+1
			em    = prop(23,etype(iel))
			pv    = prop(24,etype(iel))
			CALL deemat(dee,em,pv,type_2d)
			CALL shape_fun(fun,points,i)
			CALL shape_der(der,points,i)
			jac   = MATMUL(der,coord)
			det   = determinant(jac)
			CALL invert(jac)
			deriv = MATMUL(jac,der)
			gc    = MATMUL(fun,coord)
			CALL beemat(bee,deriv,type_2d,gc,fun)
			dvolu = det*weights(i)
			IF(type_2d=='axisymmetric') dvolu=dvolu*gc(1)*twopi
			if(epstn(kgaus)/=zero.and.nst/=3) then
				call vmdpl(dee,strsg(:,kgaus),pl)
				dee=dee-pl
			end if
			km = km + MATMUL(MATMUL(TRANSPOSE(bee),dee),bee)*dvolu
		END DO gauss_pts_4
		CALL fsparv(kv,km,g,kdiag)
	END DO elements_4

	!END IF
!------------------------------------------------------------------------------
!7.0      Predict Displacement, Velocity and Acceleration
!------------------------------------------------------------------------------

	 u0 = d		
	u10 = v
      d = d + dtim * v + a0 * a	    !new predictor
	  v = v + a1 * a		        !new predictor
	u20 = a4 * (u0 - d)
 
    !  d = d + dtim * v + a0 * a   !new predictor
	!  v = v + a1 * a		       !new predictor
	! u0 = d		
	!u10 = v
	!u20 = a4 * (u0 - d)

	rt  = zero
	rt  = loads * load(time)
	cv  = aalfa*mv + beeta * kv
	kp  = a7 * a4 * mv + kv * (a5 * a4 + one)
	CALL sparin(kp,kdiag)

!------------------------------------------------------------------------------
!8.0      Iteration Loop
!------------------------------------------------------------------------------

	iteration1: DO iters=1,miter

!------------------------------------------------------------------------------
!9.0      Update Stresses
!------------------------------------------------------------------------------

		eload = zero
		kgaus = 0
		gc    = one
		elements_5: DO iel=1,nels
		    sy    = prop(27,etype(iel))
            hards = prop(28,etype(iel))
            em    = prop(23,etype(iel))
			pv    = prop(24,etype(iel))
			CALL deemat(dee,em,pv,type_2d)
			num   = g_num(:,iel)
			coord = TRANSPOSE(g_coord(:,num))
			g     = g_g(:,iel)
			eld   = u0(g)
			bload = zero
			gauss_pts_5: DO i=1,nip
				kgaus = kgaus+1
				CALL shape_fun(fun,points,i)
				CALL shape_der(der,points,i)
				jac   = MATMUL(der,coord)
				det   = determinant(jac)
				CALL invert(jac)
				deriv = MATMUL(jac,der)
				gc    = MATMUL(fun,coord)
				CALL beemat(bee,deriv,type_2d,gc,fun)
				dvolu = det*weights(i)
				IF(type_2d=='axisymmetric')dvolu=dvolu*gc(1)*twopi
				stran = MATMUL(bee,eld)
				stran = stran-strag(:,kgaus)
				strag(:,kgaus)=strag(:,kgaus)+stran
				desig = MATMUL(dee,stran)
                
				!plasticity 
				stres=strsg(:,kgaus)+desig
				call invar(stres,sigm,dsbar,theta)
				f=dsbar-effst(kgaus)
    
				IF(f>0) THEN
					sy=f/(three*dee(3,3)+hards)
		            epstn(kgaus)=epstn(kgaus)+sy
					sy=effst(kgaus)+hards*sy
                    
                    IF(nst==3) THEN
					    stres(1)=(sy/dsbar)*(stres(1)-sigm)+sigm
					    stres(2)=(sy/dsbar)*(stres(2)-sigm)+sigm
					    stres(3)=(sy/dsbar)* stres(3)
                    ELSEIF(nst==4) THEN
					    stres(1)=(sy/dsbar)*(stres(1)-sigm)+sigm
					    stres(2)=(sy/dsbar)*(stres(2)-sigm)+sigm
					    stres(4)=(sy/dsbar)*(stres(4)-sigm)+sigm
					    stres(3)=(sy/dsbar)* stres(3) 
                    ELSEIF(nst==6) THEN
					    stres(1)=(sy/dsbar)*(stres(1)-sigm)+sigm
					    stres(2)=(sy/dsbar)*(stres(2)-sigm)+sigm
					    stres(3)=(sy/dsbar)*(stres(3)-sigm)+sigm
					    stres(4)=(sy/dsbar)* stres(4)
					    stres(5)=(sy/dsbar)* stres(5)
					    stres(6)=(sy/dsbar)* stres(6)
                    END IF
					!epstn(kgaus)=one
					effst(kgaus)=sy
				    call invar(stres,sigm,dsbar,theta)	
				END IF
				strsg(:,kgaus)=stres
				!end plasticity 1

				bload = bload + MATMUL(TRANSPOSE(bee),strsg(:,kgaus))*dvolu
			END DO gauss_pts_5
			eload(g)=eload(g)+bload
		END DO elements_5

!------------------------------------------------------------------------------
!10.0      Update displacement and Check Convergence
!------------------------------------------------------------------------------

		eload(0) = zero
		CALL linmul_sky(cv,u10,vc,kdiag)
		CALL linmul_sky(mv,u20,kd,kdiag)
		vc  = rt - vc
		kd  = vc + kd - eload
		CALL spabac(kp,kd,kdiag)
		u0  =  u0 + kd
		u20 = (u0 - d) * a4
		u10 =  v  + a3 * u20

		kd(0) = zero
		u0(0) = zero ; u10(0) = zero ; u20(0) = zero
		resid = DOT_PRODUCT(kd,kd)
		resid = DSQRT(resid)
		retot = DOT_PRODUCT(u0,u0)
		retot = DSQRT(retot)
		ratio = resid / retot
		WRITE(11,'(a,i5,a,a,1e12.6)') "Iter: ",iters,"  ","Conv. ",ratio
		IF(ratio<=toler) EXIT
		IF(iters==miter) STOP

	END DO iteration1

	d = u0
	v = u10
	a = u20

!------------------------------------------------------------------------------
!11.0      Output Results
!------------------------------------------------------------------------------

	WRITE(12,'(I5,1E12.4,<ndim>E12.4)')nres,time,d(nf(:,nres))

END DO timesteps


STOP
CONTAINS


FUNCTION load(t) RESULT(load_result)

!-----------------------Load-time function--------------------------------
 IMPLICIT NONE     
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::t
 REAL(iwp)::load_result
 load_result=1.0_iwp


RETURN
END FUNCTION load

end program


			
		!
		! Reference implementation of the KT17 model of Mercury's magnetospheric
		! magnetic field in FORTRAN 90 programming language.
		!
		program kt17
		implicit none
		integer*4 id
		integer*4 line,stat
		real*8 f,x(4),y(4),z(4),bx(4),by(4),bz(4),vx(4),vy(4),vz(4),r(1),bt(1),kx(4),ky(4),kz(4),bxz(1)
		real*8 m,m_electron,m_proton,delta_t,c,q,Kinect_Energy,Pitch_Angle,Velocity,Current_Energy
		real*8 Current_mu,Initial_mu,Current_Pitch,Gyro_Angle
		real*8 grav
		integer*4 i,j,k
		include 'kt17_common.f90'
		include 'kt17_param.f90'
		rhel=0.387
		act=50.0
		x=-1.20d0
		y=0.0d0
		z=0.0d0
		r=sqrt(x(1)**2+y(1)**2+(z(1)+0.2)**2)
		f=2.06873-0.00279*act
		rss=f*rhel**(1.0/3.0)
		print *,rss
		tamp1=6.4950+0.0229*act
		tamp2=1.6245+0.0088*act
		m_electron=9.1d-31!kg
		m_proton=1.67d-27!kg
		q=1.6d-19!C
		c=3.0d8!m/s
		grav=3.7!m/s^2
		m=m_proton
		Current_Energy=0.0d0
		Kinect_Energy=5d3*(1.6d-19)!(eV to)Joule
		Pitch_Angle=50.0d0!degree
		Gyro_Angle=-45.0d0!degree, angle between vx and vy
		Velocity=c*sqrt(1-(m**2)*(c**4)/((m*(c**2)+Kinect_Energy)**2))
		call kt17_bfield(1,x(1),y(1),z(1),bx(1),by(1),bz(1))!determine the magnetic field line vector
		bt(1)=sqrt(bx(1)**2+by(1)**2+bz(1)**2)
		bxz(1)=sqrt(bx(1)**2+bz(1)**2)
		!assume the particle's velocity's y component is 0.
		vx(1)=Velocity*&
		(cos(Pitch_Angle*4*atan(1.0d0)/180)*bxz(1)*bx(1)+&
		sqrt(bz(1)**4+(bz(1)*bx(1))**2-&
		(cos(Pitch_Angle*4*atan(1.0d0)/180)*bxz(1)*bz(1))**2))/&
		(bx(1)**2+bz(1)**2)*cos(Gyro_Angle*4*atan(1.0d0)/180)
		vy(1)=vx(1)*tan(Gyro_Angle*4*atan(1.0d0)/180)
		vz(1)=sqrt(Velocity**2-vx(1)**2-vy(1)**2)
		
		print *,vx(1),vy(1),vz(1)
		Initial_mu=m*(vx(1)**2+vy(1)**2+vz(1)**2-((bx(1)*vx(1)+by(1)*vy(1)+bz(1)*vz(1))/bt(1))**2)/(bt(1)*1d-9*(2)*1.6d-19)
		open(1,file="KT17_test_particle_result_moderate_pdyn.txt",status='old')
		do i=1,5000000
			Current_Energy=m*(vx(1)**2+vy(1)**2+vz(1)**2)/(2*1.6d-19)
			Current_mu=m*(vx(1)**2+vy(1)**2+vz(1)**2-((bx(1)*vx(1)+by(1)*vy(1)+bz(1)*vz(1))/bt(1))**2)/(bt(1)*1d-9*(2)*1.6d-19)/Initial_mu
			Current_Pitch=acos(((bx(1)*vx(1)+by(1)*vy(1)+bz(1)*vz(1))/(bt(1)*sqrt(vx(1)**2+vy(1)**2+vz(1)**2))))*180/(3.1415926)
			print *,x(1),y(1),z(1),vx(1)/1e3,vy(1)/1e3,vz(1)/1e3,bx(1),by(1),bz(1),Current_mu,Current_Energy
			write(1,*) x(1),y(1),z(1),bt(1),Current_mu,Current_Energy
			do j=1,1 !output per 1 step
				!write(1,*) x(1),y(1),z(1),vx(1)/1e3,vy(1)/1e3,vz(1)/1e3,bx(1),by(1),bz(1),Current_mu,Current_Energy
				call kt17_bfield(1,x(1),y(1),z(1),bx(1),by(1),bz(1))
				
				r=sqrt(x(1)**2+y(1)**2+(z(1)+0.196)**2)
				bt=sqrt(bx(1)**2+by(1)**2+bz(1)**2)
				
				delta_t=0.06558/bt(1) !0.1% of gyro-period
				
				kx(1)=(vy(1)*bz(1)-vz(1)*by(1))*q/(1d9*m) - grav * x(1) / (r(1)**3)
				ky(1)=(vz(1)*bx(1)-vx(1)*bz(1))*q/(1d9*m) - grav * y(1) / (r(1)**3)
				kz(1)=(vx(1)*by(1)-vy(1)*bx(1))*q/(1d9*m) - grav * (z(1)+0.196) / (r(1)**3)
				
				vx(2)=vx(1)+kx(1)*delta_t/2
				vy(2)=vy(1)+ky(1)*delta_t/2
				vz(2)=vz(1)+kz(1)*delta_t/2
				
				x(2)=x(1)+vx(1)*delta_t/(2*2440d3)
				y(2)=y(1)+vy(1)*delta_t/(2*2440d3)
				z(2)=z(1)+vz(1)*delta_t/(2*2440d3)
				
				call kt17_bfield(1,x(2),y(2),z(2),bx(2),by(2),bz(2))
				
				kx(2)=(vy(2)*bz(2)-vz(2)*by(2))*q/(1d9*m) - grav * x(2) / (r(1)**3)
				ky(2)=(vz(2)*bx(2)-vx(2)*bz(2))*q/(1d9*m) - grav * y(2) / (r(1)**3)
				kz(2)=(vx(2)*by(2)-vy(2)*bx(2))*q/(1d9*m) - grav * (z(2)+0.196) / (r(1)**3)
				
				vx(3)=vx(1)+kx(2)*delta_t/2
				vy(3)=vy(1)+ky(2)*delta_t/2
				vz(3)=vz(1)+kz(2)*delta_t/2
				
				x(3)=x(1)+vx(2)*delta_t/(2*2440d3)
				y(3)=y(1)+vy(2)*delta_t/(2*2440d3)
				z(3)=z(1)+vz(2)*delta_t/(2*2440d3)
				
				call kt17_bfield(1,x(3),y(3),z(3),bx(3),by(3),bz(3))
				
				kx(3)=(vy(3)*bz(3)-vz(3)*by(3))*q/(1d9*m) - grav * x(3) / (r(1)**3)
				ky(3)=(vz(3)*bx(3)-vx(3)*bz(3))*q/(1d9*m) - grav * y(3) / (r(1)**3)
				kz(3)=(vx(3)*by(3)-vy(3)*bx(3))*q/(1d9*m) - grav * (z(3)+0.196) / (r(1)**3)
				
				vx(4)=vx(1)+kx(3)*delta_t
				vy(4)=vy(1)+ky(3)*delta_t
				vz(4)=vz(1)+kz(3)*delta_t
				
				x(4)=x(1)+vx(3)*delta_t/(2440d3)
				y(4)=y(1)+vy(3)*delta_t/(2440d3)
				z(4)=z(1)+vz(3)*delta_t/(2440d3)
				
				call kt17_bfield(1,x(4),y(4),z(4),bx(4),by(4),bz(4))
				
				kx(4)=(vy(4)*bz(4)-vz(4)*by(4))*q/(1d9*m) - grav * x(4) / (r(1)**3)
				ky(4)=(vz(4)*bx(4)-vx(4)*bz(4))*q/(1d9*m) - grav * y(4) / (r(1)**3)
				kz(4)=(vx(4)*by(4)-vy(4)*bx(4))*q/(1d9*m) - grav * (z(4)+0.196) / (r(1)**3)
				
				x(1)=x(1)+(vx(1)+2*vx(2)+2*vx(3)+vx(4))*delta_t/(6*2440d3)
				y(1)=y(1)+(vy(1)+2*vy(2)+2*vy(3)+vy(4))*delta_t/(6*2440d3)
				z(1)=z(1)+(vz(1)+2*vz(2)+2*vz(3)+vz(4))*delta_t/(6*2440d3)
				
				vx(1)=vx(1)+delta_t/(6)*(kx(1)+2*kx(2)+2*kx(3)+kx(4))
				vy(1)=vy(1)+delta_t/(6)*(ky(1)+2*ky(2)+2*ky(3)+ky(4))
				vz(1)=vz(1)+delta_t/(6)*(kz(1)+2*kz(2)+2*kz(3)+kz(4))
				
				
				if((bx(1)**2 .lt. 1.0d-3) .and. (by(1)**2 .lt. 1.0d-3) .and. (bz(1)**2 .lt. 1.0d-3)) then
					print *,"magnetopause shadowing"
					exit
				end if
				if((bx(2)**2 .lt. 1.0d-3) .and. (by(2)**2 .lt. 1.0d-3) .and. (bz(2)**2 .lt. 1.0d-3)) then
					print *,"magnetopause shadowing"
					exit
				end if
				if((bx(3)**2 .lt. 1.0d-3) .and. (by(3)**2 .lt. 1.0d-3) .and. (bz(3)**2 .lt. 1.0d-3)) then
					print *,"magnetopause shadowing"
					exit
				end if
				if((bx(4)**2 .lt. 1.0d-3) .and. (by(4)**2 .lt. 1.0d-3) .and. (bz(4)**2 .lt. 1.0d-3)) then
					print *,"magnetopause shadowing"
					exit
				end if
				
				if(x(1)**2+y(1)**2+(z(1)+0.196)**2 .lt. 1.0d0) then
					print *,"loss cone,solid planet"
					exit
				end if
			end do
			if((bx(1)**2 .lt. 1.0d-3) .and. (by(1)**2 .lt. 1.0d-3) .and. (bz(1)**2 .lt. 1.0d-3)) then
				exit
			end if
			if(x(1)**2+y(1)**2+(z(1)+0.196)**2 .lt. 1.0d0) then
				exit
			end if
		end do
		end
		!
		!-----------------------------------------------------------------------
		! Subroutine KT17_BFIELD
		!-----------------------------------------------------------------------
		!
		subroutine kt17_bfield(n,x_a,y_a,z_a,bx_a,by_a,bz_a)
		implicit none
		include 'kt17_common.f90'
		integer*4 n,id_a(n),mode,msm,noshield,id,i
		real*8 x_a(n),y_a(n),z_a(n)
		real*8 x,y,z
		real*8 x_mso,y_mso,z_mso
		real*8 bx_mso,by_mso,bz_mso
		real*8 bx_msm,by_msm,bz_msm
		real*8 bx_dcf,by_dcf,bz_dcf
		real*8 bx_dsk,by_dsk,bz_dsk
		real*8 bx_slb,by_slb,bz_slb
		real*8 bx_a(n),by_a(n),bz_a(n)
		real*8 kappa,kappa3
		real*8 fi,gradfix,gradfiy,gradfiz
		real*8 fx,fy,fz,hx,hy,hz
		!
		! initialize variables
		!
		kappa=r0/rss
		kappa3=kappa**3
		!
		! magnetic field computation
		!
		do i=1,n
		x=x_a(i)
		y=y_a(i)
		z=z_a(i)
		x=x*kappa
		y=y*kappa
		z=z*kappa
		call kt17_mpdist(0,x,y,z,fi,id,gradfix,gradfiy,gradfiz)
		if (fi .lt. mptol) id=1
		if (noshield .eq. 1) id=1
		if (id .eq. 1) then
		bx_dcf=0.0d0
		by_dcf=0.0d0
		bz_dcf=0.0d0
		bx_dsk=0.0d0
		by_dsk=0.0d0
		bz_dsk=0.0d0
		bx_slb=0.0d0
		by_slb=0.0d0
		bz_slb=0.0d0
		fx=0.0d0
		fy=0.0d0
		fz=0.0d0
		hx=0.0d0
		hy=0.0d0
		hz=0.0d0
		call kt17_dipole(x,y,z,fx,fy,fz)
		call kt17_shield(n_dipshld,r_dipshld,x,y,z,hx,hy,hz)
		bx_dcf=kappa3*(fx+hx)
		by_dcf=kappa3*(fy+hy)
		bz_dcf=kappa3*(fz+hz)
		fx=0.0d0
		fy=0.0d0
		fz=0.0d0
		hx=0.0d0
		hy=0.0d0
		hz=0.0d0
		call kt17_taildisk(x,y,z,fx,fy,fz)
		call kt17_shield(n_diskshld,r_diskshld,x,y,z,hx,hy,hz)
		bx_dsk=tamp1*(fx+hx)
		by_dsk=tamp1*(fy+hy)
		bz_dsk=tamp1*(fz+hz)
		fx=0.0d0
		fy=0.0d0
		fz=0.0d0
		hx=0.0d0
		hy=0.0d0
		hz=0.0d0
		call kt17_tailslab(x,y,z,fx,fy,fz)
		call kt17_shield(n_slabshld,r_slabshld,x,y,z,hx,hy,hz)
		bx_slb=tamp2*(fx+hx)
		by_slb=tamp2*(fy+hy)
		bz_slb=tamp2*(fz+hz)
		bx_msm=bx_dcf+bx_dsk+bx_slb
		by_msm=by_dcf+by_dsk+by_slb
		bz_msm=bz_dcf+bz_dsk+bz_slb
		bx_a(i)=bx_msm
		by_a(i)=by_msm
		bz_a(i)=bz_msm
		else
		bx_a(i)=1.0d-8
		by_a(i)=1.0d-8
		bz_a(i)=1.0d-8
		endif
		enddo
		return
		end
		!
		!-----------------------------------------------------------------------
		! Subroutine KT17_DIPOLE
		!-----------------------------------------------------------------------
		!
		subroutine kt17_dipole(xmsm,ymsm,zmsm,bx,by,bz)
		! calculates components of dipole field
		!
		! input parameters: x,y,z - msm coordinates in rm (1 rm = 2440 km)
		!
		! output parameters: bx,by,bz - field components in msm system, in nanotesla.
		implicit none
		include 'kt17_common.f90'
		real*8 xmsm,ymsm,zmsm
		real*8 bx,by,bz
		real*8 psi,sps,cps
		real*8 p,u,v,t,q
		! dipole tilt
		psi=0.0d0
		sps=sin(psi/57.29577951d0)
		cps=sqrt(1.0d0-sps**2)
		! compute field components
		p=xmsm**2
		u=zmsm**2
		v=3.0d0*zmsm*xmsm
		t=ymsm**2
		q=mu/sqrt(p+t+u)**5
		bx=q*((t+u-2.0d0*p)*sps-v*cps)
		by=-3.0d0*ymsm*q*(xmsm*sps+zmsm*cps)
		bz=q*((p+t-2.0d0*u)*cps-v*sps)
		return
		end
		!
		!-----------------------------------------------------------------------
		! Subroutine KT17_MPDIST
		!-----------------------------------------------------------------------
		!
		subroutine kt17_mpdist(mode,x,y,z,fi,id,gradfix,gradfiy,gradfiz)
		implicit none
		include 'kt17_common.f90'
		integer*4 mode,id
		real*8 x,y,z,fi,gradfix,gradfiy,gradfiz
		real*8 rho2,r,rho
		real*8 ct,st,t,sp,cp
		real*8 rm,drm_dt
		real*8 gradfir,gradfit,gradfip
		rho2=y**2+z**2
		r=sqrt(x**2+rho2)
		rho=sqrt(rho2)
		id=1
		if (rho .gt. 1.0d-8) then ! not on the x-axis - no singularities to worry
		! about
		ct=x/r
		st=rho/r
		t=atan2(st,ct)
		sp=z/rho
		cp=y/rho
		else ! on the x-axis
		if (x .gt. 0.0d0) then ! on the dayside
		ct=x/r
		st=1.0d-8/r ! set rho=10**-8, to avoid singularity of
		! grad_fi (if mode=1, see gradfip=... below)
		t=atan2(st,ct)
		sp=0.0d0
		cp=1.0d0
		else ! on the tail axis! to avoid singularity:
		fi=-1000.0d0 ! assign rm=1000 (a conventional substitute
		! value)
		return ! and exit
		endif
		endif
		rm=r0/sqrt(alfa*(1.0d0+ct)) ! standard form of shue et al.,1997,
		! magnetopause model
		if (rm .lt. r) id=-1
		fi=r-rm
		if (mode .eq. 0) return ! skip calculation of the gradient of fi
		drm_dt=0.25d0*rm**3/r0**2*st
		gradfir=1.0d0
		gradfit=-drm_dt/r
		gradfip=0.0d0 ! axial symmetry
		gradfix=gradfir*ct-gradfit*st
		gradfiy=(gradfir*st+gradfit*ct)*cp-gradfip*sp
		gradfiz=(gradfir*st+gradfit*ct)*sp+gradfip*cp
		return
		end
		!
		!-----------------------------------------------------------------------
		! Subroutine KT17_SHIELD
		!-----------------------------------------------------------------------
		!
		subroutine kt17_shield(n,r,x,y,z,bx,by,bz)
		implicit none
		integer*4 jmax,kmax,j,k,n,o
		real*8 r(n),x,y,z
		real*8 c(n),p(n),cypj,sypj,szpk,czpk,sqpp,epp
		real*8 hx,hy,hz
		real*8 bx,by,bz
		o=nint(-0.5+sqrt(n+0.25))
		c(1:o*o)=r(1:o*o)
		p(1:o)=r(o*o+1:o*o+o)
		jmax=o
		kmax=o
		bx=0.0d0
		by=0.0d0
		bz=0.0d0
		do j=1,jmax
		do k=1,kmax
		cypj=cos(y*p(j))
		sypj=sin(y*p(j))
		szpk=sin(z*p(k))
		czpk=cos(z*p(k))
		sqpp=sqrt(p(j)**2+p(k)**2)
		epp=exp(x*sqpp)
		hx=-sqpp*epp*cypj*szpk
		hy=+epp*sypj*szpk*p(j)
		hz=-epp*cypj*czpk*p(k)
		bx=bx+hx*c((j-1)*kmax+k)
		by=by+hy*c((j-1)*kmax+k)
		bz=bz+hz*c((j-1)*kmax+k)
		enddo
		enddo
		return
		end
		!
		!-----------------------------------------------------------------------
		! Subroutine KT17_TAILDISK
		!-----------------------------------------------------------------------
		!
		subroutine kt17_taildisk(xmsm,ymsm,zmsm,bx,by,bz)
		! calculates msm components of the field from a t01-like 'long-module' equatorial
		! current disk with a 'hole' in the center and a smooth inner edge
		! (see tsyganenko, jgra, v107, no a8, doi 10.1029/2001ja000219, 2002, fig.1, right
		! panel).
		!
		!------------input parameters:
		!
		! d0 - basic (minimal) half-thickness
		! deltadx - sunward expansion factor for the current sheet thickness
		! deltady - flankward expansion factor for the current sheet thickness
		! x,y,z - msm coordinates
		!
		!------------output parameters:
		! bx,by,bz - field components in msm system, in nanotesla.
		implicit none
		include 'kt17_common.f90'
		integer*4 nr3,i
		real*8 xmsm,ymsm,zmsm,bx,by,bz
		real*8 f(n_taildisk),b(n_taildisk),c(n_taildisk)
		real*8 xshift,sc,x,y,z,d0_sc,deltadx_sc,deltady_sc
		real*8 rho,drhodx,drhody
		real*8 dex,d,dddy,dddx
		real*8 dzeta,ddzetadx,ddzetady,ddzetadz
		real*8 bi,ci,s1,s2,ds1drho
		real*8 ds2drho,ds1ddz,ds2ddz
		real*8 ds1dx,ds1dy,ds1dz,ds2dx,ds2dy,ds2dz
		real*8 s1ts2,s1ps2,s1ps2sq,fac1,as,dasds1,dasds2
		real*8 dasdx,dasdy,dasdz
		nr3=n_taildisk/3
		f(1:nr3)=r_taildisk(1:nr3)
		b(1:nr3)=r_taildisk(nr3+1:2*nr3)
		c(1:nr3)=r_taildisk(2*nr3+1:3*nr3)
		xshift=0.3d0 ! shift the center of the disk to the dayside by xshift
		sc=7.0d0 ! renormalize length scales
		x=(xmsm-xshift)*sc
		y=ymsm*sc
		z=zmsm*sc
		d0_sc=d0*sc
		deltadx_sc=deltadx*sc
		deltady_sc=deltady*sc
		rho=sqrt(x**2+y**2)
		drhodx=x/rho
		drhody=y/rho
		dex=exp(x/7.0d0)
		d=d0_sc+deltady_sc*(y/20.0d0)**2+deltadx_sc*dex ! the last term makes the
		! sheet thicken sunward, to
		dddy=deltady_sc*y*0.005d0 ! avoid problems in the
		! subsolar region
		dddx=deltadx_sc/7.0d0*dex
		dzeta=sqrt(z**2+d**2) ! this is to spread out the
		! sheet in z direction over
		! finite thickness 2d
		ddzetadx=d*dddx/dzeta
		ddzetady=d*dddy/dzeta
		ddzetadz=z/dzeta
		bx=0.0d0
		by=0.0d0
		bz=0.0d0
		do i=1,5
		bi=b(i)
		ci=c(i)
		s1=sqrt((rho+bi)**2+(dzeta+ci)**2)
		s2=sqrt((rho-bi)**2+(dzeta+ci)**2)
		ds1drho=(rho+bi)/s1
		ds2drho=(rho-bi)/s2
		ds1ddz=(dzeta+ci)/s1
		ds2ddz=(dzeta+ci)/s2
		ds1dx=ds1drho*drhodx+ds1ddz*ddzetadx
		ds1dy=ds1drho*drhody+ds1ddz*ddzetady
		ds1dz=ds1ddz*ddzetadz
		ds2dx=ds2drho*drhodx+ds2ddz*ddzetadx
		ds2dy=ds2drho*drhody+ds2ddz*ddzetady
		ds2dz=ds2ddz*ddzetadz
		s1ts2=s1*s2
		s1ps2=s1+s2
		s1ps2sq=s1ps2**2
		fac1=sqrt(s1ps2sq-(2.0d0*bi)**2)
		as=fac1/(s1ts2*s1ps2sq)
		dasds1=(1.0d0/(fac1*s2)-as/s1ps2*(s2*s2+s1*(3.0d0*s1+4.0d0*s2)))/(s1*s1ps2)
		dasds2=(1.0d0/(fac1*s1)-as/s1ps2*(s1*s1+s2*(3.0d0*s2+4.0d0*s1)))/(s2*s1ps2)
		dasdx=dasds1*ds1dx+dasds2*ds2dx
		dasdy=dasds1*ds1dy+dasds2*ds2dy
		dasdz=dasds1*ds1dz+dasds2*ds2dz
		bx=bx-f(i)*x*dasdz
		by=by-f(i)*y*dasdz
		bz=bz+f(i)*(2.0d0*as+x*dasdx+y*dasdy)
		enddo
		return
		end
		!
		!-----------------------------------------------------------------------
		! Subroutine KT17_TAILSLAB
		!-----------------------------------------------------------------------
		!
		subroutine kt17_tailslab(xmsm,ymsm,zmsm,bx,by,bz)
		! calculates msm components of the field from an equatorial harris-type current
		! sheet, slowly expanding sunward
		!
		!------------input parameters:
		!
		! d0 - basic (minimal) half-thickness
		! deltadx - sunward expansion factor for the current sheet thickness
		! deltady - flankward expansion factor for the current sheet thickness
		! scalex - e-folding distance for the sunward expansion of the current sheet
		! scaley - scale distance for the flankward expansion of the current sheet
		! zshift - z shift of image sheets
		! x,y,z - msm coordinates
		!
		!------------output parameters:
		! bx,by,bz - field components in msm system, in nanotesla.
		implicit none
		include 'kt17_common.f90'
		real*8 xmsm,ymsm,zmsm,bx,by,bz
		real*8 d,dddx,zpzi,zmzi
		d=d0+deltadx*exp(xmsm/scalex)+deltady*(ymsm/scaley)**2
		dddx=deltadx/scalex*exp(xmsm/scalex)
		zpzi=zmsm+zshift
		zmzi=zmsm-zshift
		bx=(tanh(zmsm/d)-0.5d0*(tanh(zmzi/d)+tanh(zpzi/d)))/d
		by=0.0d0
		bz=(zmsm*tanh(zmsm/d)-0.5d0*(zmzi*tanh(zmzi/d)+zpzi*tanh(zpzi/d)))*dddx/d**2
		return
		end



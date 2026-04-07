!hi
program jpl
	implicit double precision (a-h,o-z)
	parameter (nn=50000,n0=50,n1=2,n2=4,n3=18,n4=8,n5=3,n6=95,nh=50,dnk=6.0d0,dnkconv=1.15)
	
	dimension xw(nn),yw(nn),x(nn),y(nn)
	x=0.0d0
	y=0.0d0
	yw=0.0d0
	xw=0.0d0
	
	!nx=n0+n1+n2+n3+n4+n5+n6+1
	dist0=7.5d0*25.4d0
	dist1=0.31d0*25.4d0
	dist2=2.2d0*25.4d0
	dist3=7.285d0*25.4d0
	ang1=45.0d0*4.0*datan(1.0d0)/180.0
	ang2=15.0d0*4.0*datan(1.0d0)/180.0
	h1=2.5d0*25.4d0
	rc=0.5d0*25.4d0
	r1=0.8d0*25.4d0

	dx=(dist0)/dfloat(n0)
	do i =1,n0+1
		xw(i)=-dist0+dfloat(i-1)*dx
		yw(i)=h1
	enddo

	dx=(dist1)/dfloat(n1)
	x0=xw(n0+1)
	print*,x0
	do i =n0+2,n0+n1+1
		xw(i)=x0+dfloat(i-(n0+1))*dx
		yw(i)=h1
	enddo

	da=ang1/dfloat(n2)
	do i=n0+n1+2,n0+n1+n2+1
		dx=r1*dsin(da*dfloat(i-(n0+n1+1)))
		dy=-r1*(1.0d0-dcos(da*dfloat(i-(n0+n1+1))))
		xw(i)=xw(n0+n1+1)+dx
		yw(i)=yw(n0+n1+1)+dy
	enddo
	xfull=(dist2-r1*dsin(ang1)-dist1)
	do i=n0+n1+n2+2,n0+n1+n2+n3+1
		x0=xw(n0+n1+n2+1)
		y0=yw(n0+n1+n2+1)
		eta=dfloat(i-(n0+n1+n2+1))/dfloat(n3)
		feta=(dexp(-dnkconv*eta)-1.0d0)/(dexp(-dnkconv)-1.0d0)
		xw(i)=x0+xfull*feta
		yw(i)=y0-(xw(i)-x0)
	enddo
	da=ang1/dfloat(n4)
	do i=n0+n1+n2+n3+2,n0+n1+n2+n3+n4+1
		x0=xw(n0+n1+n2+n3+1)
		y0=yw(n0+n1+n2+n3+1)
		xw(i)=x0+rc*(dsin(ang1)-dsin(ang1-da*dfloat(i-(n0+n1+n2+n3+1))))
		yw(i)=y0+rc*(dcos(ang1)-dcos(ang1-da*dfloat(i-(n0+n1+n2+n3+1))))
	enddo
	da=ang2/dfloat(n5)
	do i=n0+n1+n2+n3+n4+2,n0+n1+n2+n3+n4+n5+1
		x0=xw(n0+n1+n2+n3+n4+1)
		y0=yw(n0+n1+n2+n3+n4+1)
		xw(i)=x0+rc*dsin(da*dfloat(i-(n0+n1+n2+n3+n4+1)))
		yw(i)=y0+rc*(1.0d0-dcos(da*dfloat(i-(n0+n1+n2+n3+n4+1))))
	enddo
	dx=(dist3-dist2-rc*dsin(ang1)-rc*dsin(ang2))/dfloat(n6)
	do i=n0+n1+n2+n3+n4+n5+2,n0+n1+n2+n3+n4+n5+n6+1
		x0=xw(n0+n1+n2+n3+n4+n5+1)
		y0=yw(n0+n1+n2+n3+n4+n5+1)
		xw(i)=x0+dx*dfloat(i-(n0+n1+n2+n3+n4+n5+1))
		yw(i)=y0+dtan(ang2)*(xw(i)-x0)
	enddo
	open(21,file='jplwall.dat')
	do i =1,nx
		write(21,*) i,xw(i),yw(i),0.0d0
	enddo
	
	nx=n0+n1+n2+n3+n4+n5+n6+1
	ny=nh+1
	nodes=nx*ny
	neles=(nx-1)*(ny-1)
	print*, nodes,neles
	print*, 2*(nx-1+ny-1)
	open(22,file='tape7.neu')
	write(22,*) nodes,neles
	do i = 1,nx
		do j=1,ny
			n=(i-1)*ny+j
			eta=dfloat(j-1)/dfloat(ny-1)
			feta=(dexp(-dnk*eta)-1.0d0)/(dexp(-dnk)-1.0d0)
			x(n)=xw(i)
			y(n)=yw(i)*feta
			write(22,99) n,x(n),y(n),0.0d0
99			format(i8,3f18.10)
		enddo
	enddo
	
	do i=1,nx-1
		do j=1,ny-1
			n=(i-1)*(ny-1)+j
			nl1=n+(i-1)
			nl2=nl1+ny
			nl3=nl2+1
			nl4=nl1+1
			write(22,98) n,nl1,nl2,nl3,nl4
98			format(5i8)
		enddo
	enddo	
	print*,'first grid pt at: ',y(ny)-y(ny-1),' mm'
	close(21)
	close(22)
endprogram jpl


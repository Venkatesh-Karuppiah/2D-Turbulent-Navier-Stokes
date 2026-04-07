!hi
program jpl
	implicit double precision (a-h,o-z)
	parameter (nn=50000,n1=10,n2=40,n3=10,n4=5,n5=50,nh=50,dnk=9.0d0,dnkconv=1.15)
	
	dimension xw(nn),yw(nn),x(nn),y(nn)
	x=0.0d0
	y=0.0d0
	yw=0.0d0
	xw=0.0d0
	
	!nx=n0+n1+n2+n3+n4+n5+n6+1
	!dist0=7.5d0*25.4d0
	dist1=1.0d1
	dist2=4.74d1
	dist3=5.25d1
	ang1=22.33d0*datan(1.0d0)/45.0d0
	ang2=11.24d0*datan(1.0d0)/45.0d0
	h1=3.52d1
	rc=2.74d1

	dx=dist1/dfloat(n1)
	do i =1,n1+1
		xw(i)=dx*dfloat(i-1)
		yw(i)=h1
	enddo
	
	x0=xw(n1+1)
	y0=yw(n1+1)
	dx=dist2/dfloat(n2)
	do i=n1+2,n1+n2+1
	xw(i)=x0+dx*dfloat(i-(n1+1))
	yw(i)=y0-dtan(ang1)*(xw(i)-x0)
	enddo
	
	da=ang1/dfloat(n3)
	x0=xw(n1+n2+1)
	y0=yw(n1+n2+1)
	do i=n1+n2+2,n1+n2+n3+1
		xw(i)=x0+rc*(dsin(ang1)-dsin(ang1-da*dfloat(i-(n1+n2+1))))
		yw(i)=y0+rc*(dcos(ang1)-dcos(ang1-da*dfloat(i-(n1+n2+1))))
	enddo
	xt=xw(n1+n2+n3+1)
	
	da=ang2/dfloat(n4)
	x0=xw(n1+n2+n3+1)
	y0=yw(n1+n2+n3+1)
	do i=n1+n2+n3+2,n1+n2+n3+n4+1
		xw(i)=x0+rc*dsin(da*dfloat(i-(n1+n2+n3+1)))
		yw(i)=y0+rc*(1.0d0-dcos(da*dfloat(i-(n1+n2+n3+1))))
	enddo
	
	x0=xw(n1+n2+n3+n4+1)
	y0=yw(n1+n2+n3+n4+1)
	dx=(dist3)/dfloat(n5)
	do i=n1+n2+n3+n4+2,n1+n2+n3+n4+n5+1
		xw(i)=x0+dx*dfloat(i-(n1+n2+n3+n4+1))
		yw(i)=y0+dtan(ang2)*(xw(i)-x0)
	enddo
	open(21,file='nasawall.dat')
	do i =1,n1+n2+n3+n4+n5+1
	xw(i)=xw(i)-xt
		write(21,*) i,xw(i),yw(i),0.0d0
	enddo
	
	nx=n1+n2+n3+n4+n5+1
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



program comp
	implicit double precision (a-h,o-z)
	parameter (nn=60000,n0=42,n1=24,nh=75,dnk=10.0d0)
	
	dimension xw(nn),yw(nn),x(nn),y(nn),xt(nn),yt(nn)
	
	
	x=0.0d0
	y=0.0d0
	yw=0.0d0
	xw=0.0d0
	xt=0.0d0
	yt=0.0d0
	
	del=0.83!in
	height=5.0*del
	do i =1,n0+1
		xw(i)=7.0d0*del*dfloat(i-(n0+1))/dfloat(n0)
		yw(i)=0.0d0
		xt(i)=xw(i)
		yt(i)=height
	enddo
	do i=n0+2,n0+n1+1
		xw(i) = 4.0d0*del*dfloat(i-(n0+1))/dfloat(n1)
		yw(i) = xw(i)*dtan(24.0d0*dacos(-1.0d0)/180.0d0)
		xt(i)=xw(i)
		yt(i)=height
	enddo
	
	neles=(n0+n1)*nh
	nodes=(n0+n1+1)*(nh+1)
	ind=0
	open(22,file="tape7.neu")
	write(22,*) nodes,neles
	do i=1,n0+n1+1
		do j=1,nh+1
			ind=(i-1)*(nh+1)+j		
			eta=dfloat(nh-(j-1))/dfloat(nh)
			feta=1.0d0-(dexp(-dnk*eta)-1.0d0)/(dexp(-dnk)-1.0d0)
			x(ind)=xw(i)
			y(ind)=yw(i)+(yt(i)-yw(i))*feta
5			format(i8,3f18.8)
			x(ind)=x(ind)*25.4d0
			y(ind)=y(ind)*25.4d0
			write(22,5) ind,x(ind),y(ind),0.0d0
		enddo
	enddo
	print*,'dely min= ',y(2),'mm'
	print*,'y+ =',y(2)/1000.0*22.197851004721684/(6.255405724073682e-05)
	nx=n0+n1+1
	ny=nh+1
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
	
	
	
	
endprogram comp


      subroutine hbin2(x,y,cell,cnt,xcm,ycm, size, shape,
     *                rx,ry, bnd, n, cellid, weight)

C	Copyright 1991
C	Version Date:	September 16, 1994
C	Programmer:	Dan Carr
C	Indexing:	Left to right, bottom to top
C			bnd(1) rows, bnd(2) columns
C	Output:	cell ids for non empty cells, revised bnd(1)

c			optionally also return cellid(1:n)
c     Copyright (2004) Nicholas Lewin-Koh and Martin Maechler

      implicit none

      integer n, nc, cell(*), bnd(2), cellid(*)
c     cellid(*): length 1 or n
      double precision x(n),y(n),cnt(*),xcm(*),ycm(*),rx(2),ry(2)
      double precision  size, shape, weight(n)
      integer i, i1, i2, iinc
      integer j1, j2, jinc
      integer L1, L2, L, lmax, lat
      double precision c1, c2, con1, con2, dist1, dist2, w1, w2
      double precision sx, sy, xmin, ymin, xr, yr
      logical keepID

      keepID = (cellid(1) .eq. 0)
C_______Constants for scaling the data_____________________________
      xmin = rx(1)
      ymin = ry(1)
      xr = rx(2)-xmin
      yr = ry(2)-ymin
      c1 = size/xr
      c2 = size*shape/(yr*sqrt(3.))

      jinc= bnd(2)
      lat=jinc+1
      iinc= 2*jinc
      lmax=bnd(1)*bnd(2)
      con1=.25
      con2=1.0/3.0
C_______Binning loop________________________________________
      do i=1,n
        sx = c1 * (x(i) - xmin)
        sy = c2 * (y(i) - ymin)
        j1 = int(sx+.5)
        i1 = int(sy+.5)
        dist1 = (sx-j1)**2 + 3.*(sy-i1)**2
        L1 = i1*iinc + j1+1
        
        j2 = int(sx)
        i2 = int(sy)
        dist2=(sx-j2 -.5)**2 + 3.*(sy-i2 -.5)**2
        L2 = i2*iinc+ j2+lat
        
        w1 = dist2/(dist1+dist2)
        w2 = 1. - w1

        cnt(L1) = cnt(L1) + w1 * weight(i)
        cnt(L2) = cnt(L2) + w2 * weight(i)
        
        if (keepID) then
            if (w1 .gt. w2) then
                cellid(i)=L1
            else
                cellid(i)=L2
            endif
        endif
C_______Weighted average________________________________________
C        xcm(L)=xcm(L)+ weight(i)*(x(i)-xcm(L))/cnt(L)
C        ycm(L)=ycm(L)+ weight(i)*(y(i)-ycm(L))/cnt(L)
      enddo
C_______Compression of output________________________________________
      nc=0
      do L=1,lmax
        nc=nc+1
        cell(nc)=L
        cnt(nc)=cnt(L)
        xcm(nc)=xcm(L)
        ycm(nc)=ycm(L)
      enddo
      n=nc
      bnd(1)=(cell(nc)-1)/bnd(2)+1
      return
      end
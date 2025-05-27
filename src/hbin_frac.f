      subroutine hbin_frac(x,y,cell,cnt,xcm,ycm, size, shape, 
     *                     rx,ry, bnd, n, cellid, weight)
     
      implicit none

      integer n, nc, cell(*), bnd(2), cellid(*)
      double precision x(n), y(n), cnt(*), xcm(*), ycm(*)
      double precision rx(2), ry(2), size, shape, weight(n)
      integer i, i1, i2, i3, i4, j1, j2, j3, j4, iinc, jinc, lat, lmax
      integer L, L1, L2, L3
      double precision sx, sy, xmin, ymin, xr, yr
      double precision c1, c2
      double precision dist3, dist4
      double precision xA, yA, xB, yB, xC, yC
      double precision denom, w1, w2, w3
      logical keepID

      keepID = (cellid(1) .eq. 0)
C_______Constants for scaling the data_____________________________
      xmin = rx(1)
      ymin = ry(1)
      xr   = rx(2) - xmin
      yr   = ry(2) - ymin
      c1   = size/xr
      c2   = size*shape/(yr*sqrt(3.))

      jinc = bnd(2)
      lat  = jinc + 1
      iinc = 2*jinc
      lmax = bnd(1)*bnd(2)

C_______Binning loop________________________________________
      do i = 1, n
        sx = c1 * (x(i) - xmin)
        sy = c2 * (y(i) - ymin)

        !__Two original candidates
        j1 = int(sx + .5)
        i1 = int(sy + .5)
        L1 = i1*iinc + j1 + 1

        j2 = int(sx)
        i2 = int(sy)
        L2 = i2*iinc + j2 + lat

        !__Determine third candidate
        i3 = i1
        i4 = i2
        if (j2 .lt. j1) then
          j3 = j1 - 1
          j4 = j2 + 1
        else
          j3 = j1 + 1
          j4 = j2 - 1
        endif

        dist3 = (sx-j3)**2+3.*(sy-i3)**2
        dist4 = (sx-j4-.5)**2+3.*(sy-i4-.5)**2
        if (dist3 .gt. dist4) then
          if((j4>=0).and.(j4<jinc).and.(i4>=0).and.(2*i4<bnd(1)))then
            L3 = i4*iinc + j4 + lat
          else
            L3 = -1
          endif
        else
          if((j3>=0).and.(j3<jinc).and.(i3>=0).and.(2*i3<bnd(1)))then
            L3 = i3*iinc + j3 + 1
          else
            L3 = -1
          endif
        endif

        !__compute hex-center coordinates in scaled space
        xA = j1
        yA = i1
        xB = j2 + .5
        yB = i2 + .5
        if (L3 == (i3*iinc + j3 + 1)) then
          xC = j3
          yC = i3
        else
          xC = j4 + .5
          yC = i4 + .5
        endif

        !__barycentric weights
        denom = (yB - yC)*(xA - xC) + (xC - xB)*(yA - yC)
        w1 = ((yB - yC)*(sx - xC) + (xC - xB)*(sy - yC)) / denom
        w2 = ((yC - yA)*(sx - xC) + (xA - xC)*(sy - yC)) / denom
        w3 = 1. - w1 - w2
    
        if (L3 .eq. -1) then
          w3 = w1 + w2
          w1 = w1/w3
          w2 = w2/w3
          w3 = 0
        endif

        !__distribute into bins
        cnt(L1) = cnt(L1) + w1 * weight(i)
        cnt(L2) = cnt(L2) + w2 * weight(i)
        if (L3 .ne. -1) then
          cnt(L3) = cnt(L3) + w3 * weight(i)
        endif
       
        !__assign primary id
        if (keepID) then
          if (w1 >= w2 .and. w1 >= w3) then
            cellid(i) = L1
          elseif (w2 >= w3) then
            cellid(i) = L2
          else
            cellid(i) = L3
          endif
        endif
        
C_______Weighted average________________________________________
        if (cnt(L1) .gt. 0) then
          xcm(L1)=xcm(L1)+ w1*weight(i)*(x(i)-xcm(L1))/cnt(L1)
          ycm(L1)=ycm(L1)+ w1*weight(i)*(y(i)-ycm(L1))/cnt(L1)
        endif
        
        if (cnt(L2) .gt. 0) then
          xcm(L2)=xcm(L2)+ w2*weight(i)*(x(i)-xcm(L2))/cnt(L2)
          ycm(L2)=ycm(L2)+ w2*weight(i)*(y(i)-ycm(L2))/cnt(L2)
        endif
        if(L3 .ne. -1 .and. cnt(L3) .gt. 0) then
          xcm(L3)=xcm(L3)+ w3*weight(i)*(x(i)-xcm(L3))/cnt(L3)
          ycm(L3)=ycm(L3)+ w3*weight(i)*(y(i)-ycm(L3))/cnt(L3)
    	endif    
      enddo

C_____ Output
      nc = 0
      do L = 1, lmax
        nc = nc + 1
        cell(nc) = L
        cnt(nc) = cnt(L)
        xcm(nc) = xcm(L)
        ycm(nc) = ycm(L)
      enddo
      n = nc
      bnd(1) = (cell(nc)-1)/bnd(2) + 1
      return
      end
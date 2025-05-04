      subroutine hbin3(x,y,cell,cnt,xcm,ycm, size, shape, 
     *                rx,ry, bnd, n, cellid, weight)
      implicit none

      integer n, nc, cell(*), bnd(2), cellid(*)
      double precision x(n), y(n), cnt(*), xcm(*), ycm(*)
      double precision rx(2), ry(2), size, shape, weight(n)
      integer i, i1, i2, i3, i4, j1, j2, j3, j4, iinc, jinc, lat, lmax
      integer L1, L2, L3, L
      double precision sx, sy, xmin, ymin, xr, yr
      double precision c1, c2
      double precision dist1, dist2, dist3, dist4
      double precision d1, d2, d3, w1, w2, w3, wsum
      logical keepID

      keepID = (cellid(1) .eq. 0)

C_____ Scaling constants
      xmin = rx(1)
      ymin = ry(1)
      xr   = rx(2) - xmin
      yr   = ry(2) - ymin
      c1   = size      / xr
      c2   = size*shape/(yr*sqrt(3.0d0))

      jinc = bnd(2)
      lat  = jinc + 1
      iinc = 2*jinc
      lmax = bnd(1)*bnd(2)

C_____ Main binning loop
      do i = 1, n
        sx = c1 * (x(i) - xmin)
        sy = c2 * (y(i) - ymin)

        !— two original candidates —
        j1    = int(sx + 0.5d0)
        i1    = int(sy + 0.5d0)
        dist1 = (sx - j1)**2 + 3.0d0*(sy - i1)**2
        L1    = i1*iinc + j1 + 1

        j2    = int(sx)
        i2    = int(sy)
        dist2 = (sx - j2 - 0.5d0)**2 + 3.0d0*(sy - i2 - 0.5d0)**2
        L2    = i2*iinc + j2 + lat

        !— “else-branch” fallback as third candidate —
        if(j2 .lt. j1) then
          j3 = j1 - 1
        else
          j3 = j1 + 1
        endif
        i3 = i1
        
        if(j2 .lt. j1) then
          j4 = j2 + 1
        else
          j4 = j2 - 1
        endif
        i4 = i2
        
        dist3 = (sx - j3)**2 + 3.0d0*(sy - i3)**2
        dist4 = (sx - j4 - 0.5d0)**2 + 3.0d0*(sy - i4 - 0.5d0)**2
        if (dist3 .gt. dist4) then
          dist3 = dist4
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
        
        
        d1 = max(dist1, 1d-12)
        d2 = max(dist2, 1d-12)
        d3 = max(dist3, 1d-12)

        !— inverse-distance weights —
        w1 = 1d0/d1
        w2 = 1d0/d2
        
        if (L3 .ne. -1) then
          w3 = 1d0/d3
        else
          w3 = 0
        endif
        
        wsum = w1 + w2 + w3
        w1 = w1/wsum
        w2 = w2/wsum
        w3 = w3/wsum

        !— distribute into all three bins —
        cnt(L1) = cnt(L1) + w1*weight(i)
        cnt(L2) = cnt(L2) + w2*weight(i)
        
        if (L3 .ne. -1) then
          cnt(L3) = cnt(L3) + w3*weight(i)
        endif

        !— optional: primary cell for ID —
        if (keepID) then
          if (w1 .ge. w2 .and. w1 .ge. w3) then
            cellid(i) = L1
          elseif (w2 .ge. w1 .and. w2 .ge. w3) then
            cellid(i) = L2
          else
            cellid(i) = L3
          endif
        endif
      enddo

C_____ Output
      nc = 0
      do L = 1, lmax
          nc = nc + 1
          cell(nc) = L
          cnt(nc)  = cnt(L)
          xcm(nc)  = xcm(L)
          ycm(nc)  = ycm(L)
      enddo

      n        = nc
      bnd(1)   = (cell(nc)-1) / bnd(2) + 1
      return
      end

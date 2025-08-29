#
#  (C) Brian J Soher, 2020
#

# NB. this file was originally zmgs.pentium.f but was renamed due to python
#     import name constraints, AND since we aren't doing OpenMP speedups, the
#     pentium and risc versions are pretty much the same -- bjs

# FIXME bjs - may want to move this back inside zreorth?

def zmgs(n,k,V,ldv,vnew,index):
    
    #  Check for quick return
    if (k<=0) or (n<=0): return
    
    iblck = 0               # iblck = 1
    p = index[iblck]
    q = index[iblck+1]
    
    while(p <= k and p > 0 and p <= q):
        
        # NB. ndot never declared, but also never used, so ignore
        # ndot = ndot + (q-p+1)
        
        for i in range(p,q+1):        # do i=p,q
            
            # s = 0.0+0.0j
            # for j in range(n):
            #     s += V[j, i].conjugate() * vnew[j]
            #
            # bjs - Here's the fast way.
            
            s = (V[:n,i].conjugate() * vnew[:n]).sum()

            # for j in range(n):
            #     vnew[j] = vnew[j] - (s * V[j, i])
            
            vnew[:n] -= (s*V[:n,i])

        iblck += 2
        p = index[iblck]
        q = index[iblck+1]

    return

#     subroutine zmgs(n,k,V,ldv,vnew,index)
#     implicit none
#     include 'stat.h'
#     integer n,k,ldv,index(*)
#     complex*16 V(ldv,*),vnew(*)
#     integer i,j,p,q,iblck
#     complex*16 s
#     
#     c     Check for quick return
#     if ((k.le.0).or.(n.le.0)) return
#     iblck = 1
#     p = index(iblck)
#     q = index(iblck+1)
#     do while(p.le.k .and.p .gt.0 .and. p.le.q)
#        ndot = ndot + (q-p+1)
#        do i=p,q
#           s = dcmplx(0d0,0d0)
#           do j=1,n
#              s = s + dconjg(V(j,i))*vnew(j)
#           enddo
#           do j=1,n
#              vnew(j) = vnew(j) - s*V(j,i)
#           enddo
#        enddo
#        iblck = iblck + 2
#        p = index(iblck)
#        q = index(iblck+1)
#     enddo
#     end

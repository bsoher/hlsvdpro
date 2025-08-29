#
#  (C) Brian J Soher, 2020
#


#****************************************************************************
#     This simple version of MGS is faster on Pentium machines.
#
#  NB. In checking the original Makefile flags in make.inc include file, I
#      don't think this module is ever compiled in complex*16 version of
#      the lansvdp() Make, but zmgs.pentium.py is compilee.  Maybe this file
#      was miscopied, as it seems to get compiled in the 'double' version.

def dmgs(n,k,V,ldv,vnew,index):
    
    #  Check for quick return
    if k<=0 or n<=0:
        return
    
    iblck = 0               # iblck = 1
    p = index[iblck]
    q = index[iblck+1]
    
    while(p <= k and p > 0 and p <= q):
        
        # NB. ndot never declared, but also never used, so ignore
        # ndot = ndot + (q-p+1)
        
        for i in range(p,q+1):        # do i=p,q
            # s = 0.0
            # for j in range(n):
            #     s = s + V[j,i]*vnew[j]

            s = (V[:n,i]*vnew[:n]).sum()

            # for j in range(n):
            #     vnew[j] = vnew[j] - s*V[j,i]
                
            vnew[:n] -= (s * V[:n, i])

        iblck = iblck + 2
        p = index[iblck]
        q = index[iblck+1]

    return
    
#     implicit none
#     include 'stat.h'
#     integer n,k,ldv,index(*)
#     double precision V(ldv,*),vnew(*)
#     integer i,j,p,q,iblck
#     double precision s
#     
#     c     Check for quick return
#     if ((k.le.0).or.(n.le.0)) return
#     iblck = 1
#     p = index(iblck)
#     q = index(iblck+1)
#     do while(p.le.k .and.p .gt.0 .and. p.le.q)
#        ndot = ndot + (q-p+1)
#        do i=p,q
#           s = 0d0
#           do j=1,n
#              s = s + V(j,i)*vnew(j)
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

      subroutine zreorth(n,k,V,ldv,vnew,normvnew,index,alpha,work,
     c     iflag)
c
c     Orthogonalize the N-vector VNEW against a subset of the columns of
c     the N-by-K matrix V(1:N,1:K) using iterated classical or modified
c     Gram-Schmidt. LDV is the leading dimension of the array containing
c     V.
c     
c     Which columns to orthogonalize against is decided by the integer
c     array INDEX = [s_1,e_1, s_2,e_2,..., s_k,e_l, s_{l+1}], which
c     selects the columns V(:,[s_1:e_1 s_2:e_2 ... s_l:e_l]). s_{l+1}
c     must be larger than k and marks the end of INDEX.
c
c     The reorthogonalization is repeated until
c
c       ||VNEW'|| > ALPHA * ||VNEW|| , 
c
c     where VNEW' is the vector obtained by orthogonalizing VNEW.  If
c     VNEW' fails to satisfy this after 4 tries, VNEW is deemed to lie
c     numerically in the span of V(:,[s_1:e_1 s_2:e_2 ... s_l:e_l]), and
c     is set to the zero vector.
c
c     On return NORMVNEW contains ||VNEW||.
c
c     WORK is a workspace array of length at least 
c
c       max e_i-s_i+1, i=1...l.
c     
c     WORK is only used if IFLAG==1.
c
c     If IFLAG==0 then iterated modified Gram-Schmidt is used.
c     If IFLAG==1 then iterated classical Gram-Schmidt is used.
c

c     References: 
c       Aake Bjorck, "Numerical Methods for Least Squares Problems",
c       SIAM, Philadelphia, 1996, pp. 68-69.
c     
c       J.~W. Daniel, W.~B. Gragg, L. Kaufman and G.~W. Stewart, 
c       ``Reorthogonalization and Stable Algorithms Updating the
c       Gram-Schmidt QR Factorization'', Math. Comp.,  30 (1976), no.
c       136, pp. 772-795.
c
c       B. N. Parlett, ``The Symmetric Eigenvalue Problem'', 
c       Prentice-Hall, Englewood Cliffs, NJ, 1980. pp. 105-109

c     Rasmus Munk Larsen, Stanford, 1999.

c     %-----------%
c     | Arguments |
c     %-----------%
      implicit none
      include 'stat.h'
      integer n,k,ldv,iflag,index(*)
      double complex V(ldv,*),vnew(*),work(*)
      double precision normvnew

c     %------------%
c     | Parameters |
c     %------------%
      integer NTRY
      double precision one, zero
      parameter(one = 1.0d0, zero = 0.0d0, NTRY=4)
      
c     %-----------------%
c     | Local variables |
c     %-----------------%
      integer i,itry
      double precision alpha,normvnew_0
      real t2,t3
      
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
      external zgemv
      
c     %--------------------%
c     | External Functions |
c     %--------------------%
      double precision dznrm2
      external dznrm2      

      if (k.le.0 .or. n.le.0) return

      call second(t2)

      do itry=1,NTRY
         normvnew_0 = normvnew         
         if (iflag.eq.1 ) then
            call zCGS(n,k,V,ldv,vnew,index,work)
         else  
            call zMGS2(n,k,V,ldv,vnew,index)
         endif

c         ndot = ndot + k
         normvnew = dznrm2(n,vnew,1)
         if (normvnew.gt.alpha*normvnew_0) goto 9999
      enddo
 8888 normvnew = zero
c     vnew is numerically in span(V) => return vnew = (0,0,...,0)^T
      do i=1,n
         vnew(i) = zero
      enddo

 9999 call second(t3)
      treorth = treorth + (t3-t2)
      nreorth = nreorth + 1
      return
      end



c
c****************************************************************************
c

      subroutine zCGS(n,k,V,ldv,vnew,index,work)

c     Block  Gram-Schmidt orthogonalization:
c     FOR i= 1:l
c         vnew = vnew - V(:,[s_i:e_i])*(V(:,[s_i:e_i])'*vnew)
c      
c     If l=1 and s_1=1 and e_1=k then this becomes classical Gram-Schmidt.

 
c     %-----------%
c     | Arguments |
c     %-----------%
      implicit none
      include 'stat.h'
      integer n,k,ldv,index(*)
      double complex V(ldv,*),vnew(*),work(*)
c     %------------%
c     | Parameters |
c     %------------%
      double complex one, zero
      parameter(one = (1.0d0,0.0d0), zero = (0.0d0,0.0d0))
      integer i,p,q,l


      i=1
      do while(index(i).le.k .and. index(i).gt.0)
c
c     Select the next block of columns from V
c         
         p = index(i)
         q = index(i+1)
         l = q-p+1

c         ndot = ndot + l


         call zgemv('C',n,l,one,V(1,p),ldv,vnew(1),1,zero,work(1),1)

         call zgemv('N',n,l,-one,V(1,p),ldv,work(1),1,one,vnew(1),1)
         i = i+2
      enddo
      end      



c
c****************************************************************************
c


      subroutine zMGS2(n,k,V,ldv,vnew,index)
c      subroutine MGS2(n,k,V,ldv,vnew)
      implicit none
      include 'stat.h'
      integer n,k,ldv,index(*)
      double complex V(ldv,*),vnew(*)
      integer i,j,p,q,iblck
      double complex s,zdotc
      external zdotc



c     Check for quick return
      if ((k.le.0).or.(n.le.0)) return
      iblck = 1
      do while(index(iblck).le.k.and.index(iblck).gt.0)
         p = index(iblck)
         q = index(iblck+1)
         ndot = ndot + (q-p+1)
         do i=p,q
            s = (0.0d0,0.0d0)
            do j=1,n
               s = s + conjg(V(j,i))*vnew(j)
            enddo
Cc             s = - zdotc(n,V(1,i),1,vnew(1),1)

            do j=1,n
               vnew(j) = vnew(j) - s*V(j,i)
            enddo
Cc             call zaxpy(n,s,V(1,i),1,vnew(1),1)     
         enddo
         iblck = iblck + 2
      enddo
      end

c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c

      subroutine zreorth2(n,k,V,ldv,vnew,nrm,index)
      implicit none
      integer n,k,ldv,index(*)
      complex*16 V(ldv,*),vnew(*),h,s
      double precision  nrm

c     Local variables
      include 'stat.h'
      integer i,j,p,q,iblck
      double precision gamma,nrm0,thr
      parameter(gamma = 0.98d0)
      real t2,t3

c     External subroutines
      double precision dznrm2
      external dznrm2


c     
c     Modified Gram-Schmidt orthogonalization:
c     Orthogalizes vnew against the k vectors in V by the
c     iterative process     
c     
c     FOR i=1...k DO          
c       vnew = vnew - DOT( V(:,i), vnew ) * V(:,i) 
c

c     This simple version is faster on Pentium machines.
c     Compile with "g77 -O6 -funroll-all-loops -fomit-frame-pointer"
     

c     Check for quick return
      if (k.le.0 .or. n.le.0) return
      call second(t2)

      nrm0 = nrm*nrm   
      thr = gamma*nrm0
      iblck = 1
      do while(index(iblck).le.k.and.index(iblck).gt.0)
         p = index(iblck)
         q = index(iblck+1)

         ndot = ndot + (q-p+1)
         do i=p,q
            s = (0.0d0,0.0d0)
            do j=1,n
               s = s + conjg( V(j,i))*vnew(j)
            enddo
            h = s

            do j=1,n
               vnew(j) = vnew(j) - s*V(j,i)
            enddo

   
            if  ((dble(conjg(s)*s)).gt.thr)  then

               ndot = ndot+1
               s = (0.0d0,0.0d0)
      
               do j=1,n
                  s = s + conjg(V(j,i))*vnew(j)

               enddo
               h = s + h

               do j=1,n
                  vnew(j) = vnew(j) - s*V(j,i)
               enddo

            endif
            nrm0 = nrm0 - dble(conjg(h)*h)
            thr = nrm0*gamma
         enddo
         iblck = iblck + 2
      enddo
      nrm = dznrm2(n,vnew,1)
      call second(t3)
      treorth = treorth + (t3-t2)
      nreorth = nreorth + 1



      end
c
c****************************************************************************
c

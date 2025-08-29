

      subroutine zgemmina(m,n,k,A,lda,B,ldb,dwork,ldwork)

c     
c     compute  A <- A*B'
      implicit none
c      character*1 transb
      integer m,n,k,lda,ldb,ldwork 
      double precision B(ldb,*)
      double complex A(lda,*),dwork(ldwork)
      integer i,j,l,blocksize

      if((m.le.0).or.(n.le.0).or.(k.le.0)) return
      if (ldwork.lt.n) stop 'Too little workspace in DGEMM_OVWR_LEFT'
c      if (n.gt.k) stop 'n>k in DGEMM_OVWR_LEFT'
      blocksize = int(ldwork/n)

      do i=1,m-blocksize+1,blocksize
         call zdgemma(blocksize,n,k,A(i,1),lda,
     c        B,ldb,dwork,blocksize)
         do j=0,n-1
            do l=0,blocksize-1
               A(i+l,j+1) = dwork(j*blocksize+1+l)
            enddo
c            call dcopy(blocksize,dwork(j*blocksize+1),1,A(i,j+1),1)
         enddo
      enddo
      
      call zdgemma(m-i+1,n,k,A(i,1),lda,B,ldb,dwork,m-i+1)
c      call zdgemma(m,n,k,A(i,1),lda,B,ldb,dwork,m-i+1)
      do j=0,n-1
         do l=0,m-i

            A(i+l,j+1) = dwork(j*(m-i+1)+1+l)
         enddo
c         call zcopy(m-i+1,dwork(j*(m-i+1)+1),1,A(i,j+1),1)
      enddo
      return
      end


      subroutine zdgemma(m,n,k,a,lda,b,ldb,c,ldc)
c     compute  C <- A*B'     
      implicit none

      integer m,n,k,i,j,l,lda,ldb,ldc
      double precision b(ldb,*),temp
      double complex c(ldc,*),a(lda,*)
      do j=1,n
         do i=1,m
            c(i,j)=0
         end do
         do l=1,k
c la prossima se  e' trasposta 
            temp=b(j,l)
c            temp=b(l,j)
            do i=1,m
               c(i,j)=c(i,j)+temp*a(i,l)
            end do
         end do
      end do

c      do j=1,n
c         do i=1,m
c            c(i,j)=0
c         end do
c         do l=1,k
c            temp=b(l,j)
c            do i=1,m
c               c(i,j)=c(i,j)+temp*a(i,l)
c            end do
c         end do
c      end do


      end



      subroutine zgemmin(m,n,k,A,lda,B,ldb,C,ldwork)

c     
c     compute  A <- A*B'
      implicit none
c      character*1 transb
      integer m,n,k,lda,ldb,ldwork 
      double precision B(ldb,*),temp
      double complex A(lda,*), C(lda,*)
      integer i,j,l

      if((m.le.0).or.(n.le.0).or.(k.le.0)) return
      if (ldwork.lt.n) stop 'Too little workspace in DGEMM_OVWR_LEFT'


      do j=1,n
         do i=1,m
            C(i,j)=0
         end do
         do l=1,k
c la prossima se  e' trasposta 
            temp=B(j,l)
c            temp=B(l,j)
            do i=1,m
               A(i,j)=C(i,j)+temp*A(i,l)
            end do
         end do
      end do



      end



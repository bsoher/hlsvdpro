
    
       subroutine zgemmina_python(m,n,k,aaa_r,aaa_i,lda,oda,B,ldb)
C         same as zgemmina() but with dynamic allocation of dwork
      implicit none
      complex*16, allocatable :: dwork(:)
      complex*16, DIMENSION(:, :), allocatable :: A
      integer m,n,k,lda,oda,ldb,ldwork 
      real*8 B(ldb,*)
      real*8 aaa_r(lda,*), aaa_i(lda,*)
      integer ips, jps


C     PS - Not sure big to make this. I'm not convinced the original
C     Fortran gets it right, either. 
      ldwork = 90145

C       print *, "zgemmina_python: lda = ", lda, ", oda = ", oda

      allocate(dwork(ldwork))

      allocate(A(lda, oda))

      do ips =1, lda
        do jps = 1, oda
          A(ips, jps) = dcmplx(aaa_r(ips,jps), aaa_i(ips,jps))
C           print *, "A(", ips, ",", jps, ") = ", A(ips, jps)
        end do
      end do
      call FLUSH()

      call zgemmina(m, n,k,A,lda,B,ldb,dwork,ldwork)

      do ips =1, lda
        do jps = 1, oda
          aaa_r(ips, jps) = REALPART(A(ips, jps))
          aaa_i(ips, jps) = IMAGPART(A(ips, jps))
C           print *, "A(", ips, ",", jps, ") = ", A(ips, jps)
        end do
      end do

      deallocate(dwork)

      return
      end


      subroutine zgemmina(m,n,k,A,lda,B,ldb,dwork,ldwork)

c     
c     compute  A <- A*B'
      implicit none
c      character*1 transb
      integer m,n,k,lda,ldb,ldwork 
      double precision B(ldb,*)
      double complex A(lda,*),dwork(ldwork)
      integer i,j,l,blocksize
      integer ips, jps

C       print *, "inside zgemmina"
C       print *, "zgemmina: m = ", m
C       print *, "zgemmina: n = k = ", n
C       print *, "zgemmina: ldwork = ", ldwork
      if((m.le.0).or.(n.le.0).or.(k.le.0)) return
      if (ldwork.lt.n) stop 'Too little workspace in DGEMM_OVWR_LEFT'
c      if (n.gt.k) stop 'n>k in DGEMM_OVWR_LEFT'
      blocksize = int(ldwork/n)
      
C       print *, "inside zgemmina, blocksize = ", blocksize


      do i=1,m-blocksize+1,blocksize
C          print *, "inside zgemmina, i = ", i
C          call FLUSH()
C          print *, "before zdgemma1"
         call zdgemma(blocksize,n,k,A(i,1),lda,
     c        B,ldb,dwork,blocksize)

C          print *, "after zdgemma1"
C          do ips=1,ldwork
C             print *, "dwork(", ips, ") = ", dwork(ips)
C          end do

         do j=0,n-1
            do l=0,blocksize-1
C                print *, "j = ", j, "l = ", l
C                  print *, "index = ", j*blocksize+1+l
               A(i+l,j+1) = dwork(j*blocksize+1+l)
C                print *, "A(", i+l,",",j+1,") = ", A(i+l,j+1)
            enddo
c            call dcopy(blocksize,dwork(j*blocksize+1),1,A(i,j+1),1)
         enddo
      enddo

C       print *, "before zdgemma2, i = ", i
C       call FLUSH()
C       do ips = 1,25
C         do jps = 1, 25
C             print *, "A(", ips, ",", jps, ") = ", A(ips, jps)
C         end do
C       end do

C       do ips = 1,20
C         do jps = 1, 20
C             print *, "B(", ips, ",", jps, ") = ", B(ips, jps)
C         end do
C       end do

C      call FLUSH()
      

      
      call zdgemma(m-i+1,n,k,A(i,1),lda,B,ldb,dwork,m-i+1)
C       print *, "done zdgemma2"
C       call FLUSH()
C       do ips=1,ldwork
C          print *, "dwork(", ips, ") = ", dwork(ips)
C       end do
c      call zdgemma(m,n,k,A(i,1),lda,B,ldb,dwork,m-i+1)
      do j=0,n-1
         do l=0,m-i

            A(i+l,j+1) = dwork(j*(m-i+1)+1+l)
         enddo
c         call zcopy(m-i+1,dwork(j*(m-i+1)+1),1,A(i,j+1),1)
      enddo
C       do ips = 1,lda
C         do jps = 1, 50
C             print *, "A(", ips, ",", jps, ") = ", A(ips, jps)
C         end do
C       end do

      return
      end


      subroutine zdgemma(m,n,k,a,lda,b,ldb,c,ldc)
c     compute  C <- A*B'     
      implicit none

      integer m,n,k,i,j,l,lda,ldb,ldc
      double precision b(ldb,*),temp
      double complex c(ldc,*),a(lda,*)

C       print *, "inside zdgemma, m = ", m
C       print *, "inside zdgemma, n = ", n
C       print *, "inside zdgemma, k = ", k
      do j=1,n
         do i=1,m
C            print *, "zdgemma: zeroing(", i, ", ", j, ")"
            c(i,j)=0
         end do
         do l=1,k
c la prossima se  e' trasposta 
            temp=b(j,l)
C              print *, "b(", j, ", ", l, ") = ", b(j,l)
            do i=1,m
C                if (i.lt.20) then 
C                   if (l.lt.20) then
C                     print *, "c(", i, ", ", j, ") = ", c(i,j)
C                   end if
C                end if
               c(i,j)=c(i,j)+temp*a(i,l)
C                if (i.lt.25) then 
C                   if (l.lt.25) then
C                      print *, "a(", i, ", ", l, ") = ", a(i,l)
C                   end if
C                end if
            end do
         end do
      end do

C     call FLUSH()

C      print *, "supposedly exiting"
C      call EXIT()

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


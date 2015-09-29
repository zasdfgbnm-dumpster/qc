C*************************************************************
      subroutine Get1EInt(ndim)
C     
C     Get one electron integrals
C         
C*************************************************************
      implicit double precision (a-h,o-z)
c
      logical FileExist
c
      dimension oneh((NDIM*(NDIM+1))/2),ovrlp((NDIM*(NDIM+1))/2),
     &          buf(600),ibuf(600)
      DIMENSION oneh_full(NDIM,NDIM),ovrlp_full(NDIM,NDIM)
c
      ilnbuf=600
      FileExist=.false.
      inquire(file='IIII',exist=FileExist)
      if (FileExist) then
c         write(*,*) 'IIII file exists'
         open(unit=10,file='IIII',form='UNFORMATTED',
     &        access='SEQUENTIAL')
         rewind 10
         call locate(10,'ONEHAMIL')
         call zero(oneh,ldim)
         nut = 1
         do while (nut.gt.0)
            read(10) buf, ibuf, nut
            do int = 1, nut
               oneh(ibuf(int)) = buf(int)
            end do
         end do
c
         call locate(10,'OVERLAP ')
         call zero(ovrlp,ldim)
         nut = 1
         do while (nut.gt.0)
            read(10) buf, ibuf, nut
            do int = 1, nut
               ovrlp(ibuf(int)) = buf(int)
            end do
         end do
         close(10)
      else
         write(*,*) 'IIII file does not exist'
         stop
      end if
c unpack
      ITHRU = 0
      DO I = 1, NDIM
         DO J = 1, I
            ITHRU = ITHRU + 1
            oneh_full(I,J) = oneh(ITHRU)
            ovrlp_full(I,J) = ovrlp(ITHRU)
            oneh_full(J,I) = oneh(ITHRU)
            ovrlp_full(J,I) = ovrlp(ITHRU)
         END DO
      END DO
c output
      write (*,*) "one electron h"
      DO I = 1, NDIM
         DO J = 1, NDIM
            write(*,*) oneh_full(I,J)
         END DO
      END DO
      write (*,*) "overlap integrals"
      DO I = 1, NDIM
         DO J = 1, NDIM
            write(*,*) ovrlp_full(I,J)
         END DO
      END DO
      return
      end

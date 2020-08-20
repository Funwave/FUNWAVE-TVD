         REAL,DIMENSION(:,:),ALLOCATABLE :: DEPTH,X,Y
         INTEGER :: Mglob, Nglob, I, J
         REAL :: dep_channel,slope,width_channel,DX,DY
         CHARACTER(LEN = 80) :: FileName

         FileName = 'bathy_param.txt'
         OPEN(1,File=TRIM(FileName))

         READ(1,*) Mglob
         READ(1,*) Nglob
         READ(1,*) DX
         READ(1,*) DY
         READ(1,*) dep_channel
         READ(1,*) width_channel
         READ(1,*) slope
         CLOSE(1)

         ALLOCATE (DEPTH(Mglob,Nglob),X(Mglob,Nglob),Y(Mglob,Nglob))

         DO J=1,Nglob
         DO I=1,Mglob
           Depth(I,J)=dep_channel
           X(I,J)=(I-1.0)*DX
           Y(I,J)=(J-Nglob*0.5)*DY
         ENDDO
         ENDDO

         DO J=1,Nglob
         DO I=1,Mglob
            IF(abs(Y(I,J))>width_channel*0.5)THEN
              DEPTH(I,J) = dep_channel -
     &               slope*(abs(Y(I,J))-width_channel*0.5)
            ENDIF
         ENDDO
         ENDDO

         OPEN(2,FILE='depth.txt')
         DO J=1,Nglob
          WRITE(2,100)(Depth(I,J),I=1,Mglob)
         ENDDO
100      FORMAT(9000f12.6)
         CLOSE(2)

         END



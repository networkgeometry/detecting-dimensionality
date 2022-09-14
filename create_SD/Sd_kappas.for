* Program to generate a SD version of a given network with an specific
* Beta and number of dimensions. This program estimates the hidden
* degrees of each class of degree in order to force the degree
* distribution.
* cutoff added
       PROGRAM SD
        IMPLICIT DOUBLE PRECISION (x)
        CHARACTER*300 FILEK,filenet,fileoutk,arg3,arg4,arg5

        INTEGER SEED,D,N,DMAX,Nedgesmax,maxdegree
        INTEGER I,J
        DOUBLE PRECISION B,A,AUX,maxerror,PROB
        DOUBLE PRECISION PI,M,R,TA
        double precision result,c
        double precision ran2

        PARAMETER (N=100000)               !Number of nodes
        PARAMETER (Nedgesmax=20000000)               !Number of edges
        PARAMETER (DMAX=20)                  !Similarity space dimension
        PARAMETER (PI=4.D0*DATAN(1.D0))

        INTEGER deg(1:N) !la lista de grado
        INTEGER xd(1:N) !la distribucion de grado
        dimension ndegree(1:N)
        dimension nedge(1:Nedgesmax,1:2)
        dimension npunteroini(1:N)
        dimension npunterofin(1:N)
        dimension xk(1:N) !la kappa de cada clase de degree
        
        dimension xexpk(1:N) !el expected degree de cada clase de degree
        dimension X(1:N,1:DMAX) !Coordenadas de cada nodo
        LOGICAL :: file_exists
        xe=1.0D0   !Precision observed degrees
        CALL get_command_argument(1, FILEK)
        CALL get_command_argument(2, filenet)
        CALL get_command_argument(3, arg3)
        CALL get_command_argument(4, arg4)
        CALL get_command_argument(5, arg5)

        
        READ(arg3,*)D
        READ(arg4,*)B 
        READ(arg5,*)SEED
        B = dble(B)
        fileoutk = TRIM(FILEK) // "_degrees.dat"
        filenet = TRIM(filenet)//".edge"
        FILEK = TRIM(FILEK)//".edge"
        SEED=-SEED
        write(6,*)"fileoutk"
        write(6,*)fileoutk
        
********obtaining degree list
        open(1,file=FILEK,status='unknown')
        NODOS=0
        nlist=0
        do while(.true.)
          read(1,*,END=10) i,j
          ndegree(i)=ndegree(i)+1
          ndegree(j)=ndegree(j)+1
          nlist=nlist+1
          nedge(nlist,1)=i
          nedge(nlist,2)=j
          if(i.gt.NODOS)then
            NODOS=i
          endif
          if(j.gt.NODOS)then
            NODOS=j
          endif
        enddo
10      continue
        close(1)

*****writing degrees
        INQUIRE(FILE=fileoutk, EXIST=file_exists)
        if(.NOT.file_exists)then
        open(1,file=fileoutk,status='unknown')
        indexaux=1
        ndtot=0
        nodosreal=0
        xk=0.D0
        do i=1,NODOS
         if(ndegree(i).gt.0)then
           npunteroini(i)=indexaux
           npunterofin(i)=indexaux-1
           indexaux=indexaux+ndegree(i)
           ndtot=ndtot+ndegree(i)
           nodosreal=nodosreal+1
           xk=xk+ndegree(i)
           write(1,*)ndegree(i)
         endif
        enddo
        close(1)
        endif
*****reading degrees

        do i=1,NODOS
            xd(i)=0
        enddo
        OPEN(UNIT=12,FILE=fileoutk,STATUS="UNKNOWN")
        TA=0.D0
        NODOS=0
        maxdegree=1
        do while(.true.)
          read(12,*,END=20) i
          if(i.gt.maxdegree)then
            maxdegree=i
          endif
          NODOS=NODOS+1
          xk(NODOS)=NODOS
          xd(i)=xd(i)+1
          deg(NODOS)=i
          xexpk(NODOS)=0
          TA=TA+dble(i)
        enddo
20      continue
        CLOSE(12)
        TA=TA/dble(NODOS)
        WRITE(6,*) "NODES AND AVG K:",NODOS,TA
        WRITE(6,*) "maxdegree:",maxdegree

*****************
 
        write(6,*)NODOS   
        write(6,*)"FILEOUT"
        write(6,*)filenet

******assigning coordinates
        DO I=1,NODOS
          DO J=1,D+1
            X(I,J)=DSQRT(-2.D0*DLOG(DBLE(ran2(SEED))))*
     +             DCOS(2.D0*PI*ran2(SEED))
          ENDDO
          X(I,:)=X(I,:)/DSQRT(SUM(X(I,:)**2))
        ENDDO
        R=(DBLE(NODOS)*GAMMA((DBLE(D)+1.D0)/2.D0)
     +      /(2.D0*PI**((DBLE(D)+1.D0)/2.D0)))**(1.D0/DBLE(D))
  
        xce=0.D0
        xsq=0.D0
        xpe=0.D0

        nodost=0
        ncont=0
        A=TA             !Initial Avg. degree

30      continue
        DO I=1,maxdegree
          xexpk(I)=0
        ENDDO
        write(6,*)"Calculating hidden degrees..."

        npunteroini=0
        npunterofin=0
        ndegree=0

        M=GAMMA(DBLE(D)/2.D0)*B*DSIN(DBLE(D)*PI/B)
     +      /(2.D0*PI**(1.D0+DBLE(D)/2.D0)*A)

*****************
        DO I=1,maxdegree
            DO J=I,maxdegree
                !write(6,*)"xk(I),xk(J):",xk(I),xk(J)
                !write(6,*)"xd(I),xd(J):",xd(I),xd(J)

                c=R/(((M*xk(I)*xk(J)))**(1/DBLE(D)))
                call pkk(result,D,B,c)
                !write(6,*)"c,pkk:",c,result
                if (I.ne.J)then
                    xexpk(I)=xexpk(I)+xd(J)*result
                    xexpk(J)=xexpk(J)+xd(I)*result
                else
                    xexpk(I)=xexpk(I)+(xd(J)-1)*result
                endif
            ENDDO    
        ENDDO 
        
        maxerror=0
        DO I=1,maxdegree
            IF(abs(xexpk(I)-I)>maxerror)then
                maxerror=dabs(xexpk(I)-I)
            ENDIF
        ENDDO
        WRITE(6,*) "maxerror:",maxerror
        IF(maxerror>xe)then
            DO I=1,maxdegree
                xk(I)=abs(xk(I)+(I-xexpk(I))*ran2(SEED))
            ENDDO
            goto 30
        ENDIF
        
 
*****************

        AUX=0.D0
        nlist=0
        NDS=0

        OPEN(UNIT=10,FILE=filenet,STATUS="UNKNOWN")
          DO I=1,NODOS-1  !connecting the network
            DO J=I+1,NODOS
              !write(6,*)"deg(I),xk(deg(I)):",deg(I),xk(deg(I))
              PROB=1.D0/(1.D0+(R*DACOS(SUM(X(I,:)*X(J,:)))
     +             /(M*xk(deg(I))*xk(deg(J)))**(1.D0/DBLE(D)))**B)
              IF(ran2(SEED).LE.PROB)then
                WRITE(10,*) I,J
                nlist=nlist+1
                nedge(nlist,1)=i
                nedge(nlist,2)=j
                ndegree(i)=ndegree(i)+1
                ndegree(j)=ndegree(j)+1
                if(i.gt.NDS)then
                  NDS=i
                endif
                if(j.gt.NDS)then
                  NDS=j
                endif
              endif
              AUX=AUX+2.D0*PROB/DBLE(NODOS)
            ENDDO
          ENDDO
        CLOSE(10)
        indexaux=1
        ndtot=0
        nodosreal=0
        do i=1,NDS
         if(ndegree(i).gt.0)then
           npunteroini(i)=indexaux
           npunterofin(i)=indexaux-1
           indexaux=indexaux+ndegree(i)
           ndtot=ndtot+ndegree(i)
           nodosreal=nodosreal+1
         endif
        enddo
      END PROGRAM SD



      !call pkk(result,d,beta,c)

      subroutine pkk(result,d,beta,c)
      implicit none
      double precision Intnew,pi,a,b,Int0,Int1,error,c,beta
      double precision Intfirst,Intsecond,delta,result
      integer m,d
      delta=1.d-5
      pi=4.d0*datan(1.d0)
      a=0.d0
      b=pi
      m=50
      error=2.d0*delta
      do while(error.gt.delta)
       m=2*m

       call trapecis(m,Intnew,a,b,d,c,beta)
       Int0=Intnew
       call trapecis(2*m,Intnew,a,b,d,c,beta)
       Int1=Intnew
       Intfirst=(4.d0*Int1-Int0)/3.d0

       m=2*m

       call trapecis(m,Intnew,a,b,d,c,beta)
       Int0=Intnew
       call trapecis(2*m,Intnew,a,b,d,c,beta)
       Int1=Intnew
       Intsecond=(4.d0*Int1-Int0)/3.d0

       error=dabs(Intsecond-Intfirst)/dabs(Intfirst)

      enddo

       result=Intsecond*dgamma(dble(1+d)/2.d0)
     +/(dsqrt(pi)*dgamma(dble(d)/2.d0))

      return
      end
	  

      subroutine trapecis(n,Int,a,b,d,c,beta)
      double precision Int,h,x,a,b,func,c,beta
      integer n,i,d

       x=a
       Int=func(a,d,c,beta)+func(b,d,c,beta)
       h=(b-a)/dble(n)

       do i=1,n-1
        x=x+h
        Int=Int+2.d0*func(x,d,c,beta)
       enddo
       Int=Int*h/2.d0

       return
       end


      double precision function func(x,d,c,beta)
      double precision x,c,beta
      integer d

       func=(dsin(x))**(d-1)/(1.d0+(c*x)**beta)


      return
      end


******** Uniform Random generator ***********************
       FUNCTION ran2(idum)
       INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
       DOUBLE PRECISION ran2,AM,EPS,RNMX
       PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     +   IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,
     +   IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
!Long period (> 2 x 10 18 ) random number generator of L'Ecuyer with Bays-Durham shuffle
!and added safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive
!of the endpoint values). Call with idum a negative integer to initialize; thereafter, do not
!alter idum between successive deviates in a sequence. RNMX should approximate the largest
!floating value that is less than 1.
       INTEGER idum2,j,k,iv(NTAB),iy
       SAVE iv,iy,idum2
       DATA idum2/123456789/, iv/NTAB*0/, iy/0/
       if (idum.le.0) then               !Initialize.
         idum=max(-idum,1)             !Be sure to prevent idum = 0.
         idum2=idum
         do j=NTAB+8,1,-1           !Load the shuffle table (after 8 warm-ups).
           k=idum/IQ1
           idum=IA1*(idum-k*IQ1)-k*IR1
           if (idum.lt.0) idum=idum+IM1
           if (j.le.NTAB) iv(j)=idum
         enddo
         iy=iv(1)
       endif
       k=idum/IQ1                        !Start here when not initializing.
       idum=IA1*(idum-k*IQ1)-k*IR1       !Compute idum=mod(IA1*idum,IM1) without over-
       if (idum.lt.0) idum=idum+IM1      !flows by Schrage's method.
       k=idum2/IQ2
       idum2=IA2*(idum2-k*IQ2)-k*IR2     !Compute idum2=mod(IA2*idum2,IM2) likewise.
       if (idum2.lt.0) idum2=idum2+IM2
       j=1+iy/NDIV                       !Will be in the range 1:NTAB.
       iy=iv(j)-idum2                    !Here idum is shuffled, idum and idum2 are com-
       iv(j)=idum                        !bined to generate output.
       if(iy.lt.1)iy=iy+IMM1
       ran2=dmin1(AM*dble(iy),RNMX)              !Because users don't expect endpoint values.
       return
       END



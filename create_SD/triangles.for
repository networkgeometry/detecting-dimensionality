* Calcualtes number of cycles in networks

      implicit double precision(x)
      character*90 filename,fileout
      parameter  (nodosmax=130000,Nedgesmax=2000000)

      dimension nconnectivity(1:Nedgesmax)
      dimension npunteroini(1:nodosmax)
      dimension npunterofin(1:nodosmax)
      dimension ndegree(1:nodosmax)
      dimension nedge(1:Nedgesmax,1:2)
      dimension nvec(1:nodosmax)

      data npunteroini/nodosmax*0/
      data npunterofin/nodosmax*0/
      data ndegree/nodosmax*0/

*Lee red de fichero y carga en vector de dos sentidos
*****************************************************
      CALL get_command_argument(1, filename)
      CALL get_command_argument(2, fileout)

      !fileoutk='test_edgelist_0_degrees.dat'
    
      open(1,file=filename,status='unknown')
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
10    continue
      close(1)

      !      open(1,file=fileoutk,status='unknown')
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
        !       write(1,*)ndegree(i)
             endif
            enddo
       !     close(1)
            xk=xk/dble(nodosreal)

            do l=1,nlist
             i=nedge(l,1)
             j=nedge(l,2)
             npunterofin(i)=npunterofin(i)+1
             nconnectivity(npunterofin(i))=j
             npunterofin(j)=npunterofin(j)+1
             nconnectivity(npunterofin(j))=i
            enddo


*********clustering
*            write(6,*)'Calculating node clustering ...'
            nds1=0
            xckiter=0.D0
            do i=1,NODOS
            if(ndegree(i).gt.1)then
              nds1=nds1+1
              nclust=0
              do j=npunteroini(i),npunterofin(i)-1
               do k2=j+1,npunterofin(i)
                 do l=npunteroini(nconnectivity(j)),
     +              npunterofin(nconnectivity(j))
                    if(nconnectivity(l).eq.nconnectivity(k2))then
                      nclust=nclust+1
                    endif
                  enddo
                enddo
              enddo
              xckiter=xckiter
     +        +2*dble(nclust)/(dble(ndegree(i))*(dble(ndegree(i))-1.))
            endif
            enddo
            xckiter=xckiter/dble(nds1)

*********cycles of edge
            nl1=0
            squaresP=0
            pentagonsP = 0

            xceiter=0.D0
            xsqiter=0.D0
            xpiter=0.D0
            ntt=0
            nts=0
            ntp=0
            do i=1,nlist   !do edges
              nodo1=nedge(i,1)
              nodo2=nedge(i,2)
              mledge=0
              msq=0
              mpent=0

*Calculating edge clustering ...
              if((ndegree(nodo1).gt.1).and.(ndegree(nodo2).gt.1))then
              nl1=nl1+1
              nvec=0
              do j=npunteroini(nodo1),npunterofin(nodo1)
               nodo3=nconnectivity(j)
               do k2=npunteroini(nodo2),npunterofin(nodo2)
                 nodo4=nconnectivity(k2)
                 if(nodo3.eq.nodo4)then
                   mledge=mledge+1
                   nvec(mledge)=nodo3
                 endif
               enddo
              enddo
              ntt=ntt+mledge
              xceiter=xceiter+dble(mledge)/
     +         (DMIN1(dble(ndegree(nodo1)),dble(ndegree(nodo2)))-1.0D0)
              endif !endif degrees larger than 1

            enddo   !enddo edges

            xceiter=xceiter/dble(nl1)
            open(1,file=fileout,status='unknown')
            write(1,*)xceiter
            close(1)
      end

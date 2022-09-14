* Calculates number of chordless cycles in networks
* Author: M. Ángeles Serrano Moral

      implicit double precision(x)
      character*90 filename,fileout,fileoutk
      parameter  (nodosmax=130000,Nedgesmax=2000000)

      dimension nconnectivity(1:Nedgesmax)
      dimension npunteroini(1:nodosmax)
      dimension npunterofin(1:nodosmax)
      dimension ndegree(1:nodosmax)
      dimension nedge(1:Nedgesmax,1:2)
      dimension nvec(1:nodosmax)
      dimension nvec1(1:nodosmax)
      dimension nvec2(1:nodosmax)

      data npunteroini/nodosmax*0/
      data npunterofin/nodosmax*0/
      data ndegree/nodosmax*0/

*Lee red de fichero y carga en vector de dos sentidos
*****************************************************
      CALL get_command_argument(1, filename)
      CALL get_command_argument(2, fileout)


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
               write(1,*)i,ndegree(i)
             endif
            enddo
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
            xceiter=0.D0
            xsqiter=0.D0
            xpiter=0.D0
            x2ceiter=0.D0
            x2sqiter=0.D0
            x2piter=0.D0
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
                   nvec(mledge)=nodo3   !common neighbors
                 endif
               enddo
              enddo
              ntt=ntt+mledge
              xtrans=dble(mledge)/
     +         (DMIN1(dble(ndegree(nodo1)),dble(ndegree(nodo2)))-1.0D0)
              xceiter=xceiter+xtrans
              x2ceiter=x2ceiter+xtrans**2


*Calculating edge squares ...
              if((ndegree(nodo1).gt.(mledge+1)).and.
     +           (ndegree(nodo2).gt.(mledge+1)))then !if more than multiplicity

              do j=npunteroini(nodo1),npunterofin(nodo1)
               nodo3=nconnectivity(j)
               if(nodo3.ne.nodo2)then

                 nsn=0
                 do k2=1,mledge
                   if(nodo3.eq.nvec(k2))then
                     nsn=1
                     goto 40
                   endif
                 enddo
40               continue

                 if(nsn.eq.0)then
                   do k2=npunteroini(nodo2),npunterofin(nodo2)
                     nodo4=nconnectivity(k2)
                     if(nodo4.ne.nodo1)then

                       nsn2=0
                       do k3=1,mledge
                         if(nodo4.eq.nvec(k3))then
                           nsn2=1
                           goto 50
                         endif
                       enddo
50                     continue

                       if(nsn2.eq.0)then
                         do k4=npunteroini(nodo3),npunterofin(nodo3)
                           nodo5=nconnectivity(k4)
                           if(nodo5.eq.nodo4)then
                             msq=msq+1
                             goto 60
                           endif
                         enddo
60                       continue
                       endif

                     endif
                   enddo

                 endif
               endif
              enddo

              nts=nts+msq
              if(msq.gt.0)then
               nkc1=ndegree(nodo1)-mledge-1
               nkc2=ndegree(nodo2)-mledge-1
               xledge=dble(msq)/(dble(nkc1)*dble(nkc2))
              else
                xledge=0.D0
              endif
              xsqiter=xsqiter+xledge
              x2sqiter=x2sqiter+xledge**2

              endif !endif more than multiplicity

*Calculating edge pentagones ...
              if((ndegree(nodo1).gt.(mledge+1)).and.
     +           (ndegree(nodo2).gt.(mledge+1)))then !if more than m in p
!normalization: sum[(k_n1-1)*(k_n2-1)] where n1 are neighbors of i different from j, n2 are neighbors of j different from i, n1.ne.n2 (no common neighbors),and n1 and n2 not connected (do not form a square)

               ndsv1=0 !adding neighbors of node1 which are not in triangles with node2 and with degree >1
               do j2=npunteroini(nodo1),npunterofin(nodo1)
               nodo3=nconnectivity(j2)
               if(ndegree(nodo3).gt.1)then
               if(nodo3.ne.nodo2)then
                 nsn=0
                 do k3=1,mledge
                   if(nvec(k3).eq.nodo3)then
                     nsn=1    !the neighbor is in a triangle
                     goto 70
                   endif
                 enddo
70               continue

                 if(nsn.eq.0)then
                   ndsv1=ndsv1+1
                   nvec1(ndsv1)=nodo3
                 endif
               endif
               endif
               enddo

               ndsv2=0 !adding neighbors of node2 which are not in triangles with node1 and with degree >1
               do j2=npunteroini(nodo2),npunterofin(nodo2)
               nodo4=nconnectivity(j2)
               if(ndegree(nodo4).gt.1)then
               if(nodo4.ne.nodo1)then
                 nsn=0
                 do k3=1,mledge
                   if(nvec(k3).eq.nodo4)then
                     nsn=1
                     goto 80
                   endif
                 enddo
80               continue

                 if(nsn.eq.0)then
                   ndsv2=ndsv2+1
                   nvec2(ndsv2)=nodo4
                 endif
               endif
               endif
               enddo

               nnor=0
               do k3=1,ndsv1
               do k4=1,ndsv2
                 nnei=0
                 do k6=npunteroini(nvec2(k4)),npunterofin(nvec2(k4))
                   nodo6=nconnectivity(k6)
                   if(nodo6.eq.nvec1(k3))then
                     nnei=1   ! the link is in a square
                   endif
                 enddo

                 if(nnei.eq.0)then
                   do k5=npunteroini(nvec1(k3)),npunterofin(nvec1(k3))
                   nodo5=nconnectivity(k5)
                   do k6=npunteroini(nvec2(k4)),npunterofin(nvec2(k4))
                     if(nodo5.eq.nconnectivity(k6))then

                       nin=0
                       do j2=npunteroini(nodo1),npunterofin(nodo1)
                         nodo3=nconnectivity(j2)
                         if(nodo3.eq.nodo5)then
                           nin=1  !node5 is connected to node1
                           goto 90
                         endif
                       enddo
                       do j2=npunteroini(nodo2),npunterofin(nodo2)
                         nodo3=nconnectivity(j2)
                         if(nodo3.eq.nodo5)then
                           nin=1  !node5 is connected to node2
                           goto 90
                         endif
                       enddo

                       if(nin.eq.0)then
                         mpent=mpent+1
                         goto 90
                       endif

                     endif
                   enddo  !checking neighbors of neighbors of node2
90                 continue
                   enddo  !checking neighbors of neighbors of node1
                   nnor=nnor+(ndegree(nvec1(k3))-1)*
     +                       (ndegree(nvec2(k4))-1)
                 endif    !if nnei

               enddo   !enddo ndsv1
               enddo   !enddo ndsv2
                

               ntp=ntp+mpent
               if(mpent.gt.0)then
                 xpent=dble(mpent)/dble(nnor)
               else
                 xpent=0.0D0
               endif
               xpiter=xpiter+xpent
               x2piter=x2piter+xpent**2

              endif !endif more than m in p

              endif !endif degrees larger than 1

            enddo   !enddo edges

            xceiter=xceiter/dble(nl1)
            xsqiter=xsqiter/dble(nl1)
            xpiter=xpiter/dble(nl1)

            x2ceiter=x2ceiter/dble(nl1)
            xerrorace=dsqrt(x2ceiter-xceiter**2)/dsqrt(dble(nl1))
            x2sqiter=x2sqiter/dble(nl1)
            xerrorasq=dsqrt(x2sqiter-xsqiter**2)/dsqrt(dble(nl1))
            x2piter=x2piter/dble(nl1)
            xerrorap=dsqrt(x2piter-xpiter**2)/dsqrt(dble(nl1))

            open(1,file=fileout,status='unknown')
            write(1,*)xceiter,",",xsqiter,",",xpiter,",",NODOS,",",nlist
            close(1)


      end

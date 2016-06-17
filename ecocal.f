c-----------------------------------------------------------------------
c
c     -- Physical Properties and State Variables --
c                                               Author: Daisuke Kitazawa
c                                               Last Update: 2005.7.29
c
c-----------------------------------------------------------------------
      subroutine ecocal
      implicit double precision (a-h,o-z)
      include 'param.h'
c
      if(nt.ne.1) goto 200
c
c     -- density --
c
      do 10 k=1,nz
       dens(k)=1028.14d0-.0735d0*tmp(k)-.00469d0*tmp(k)*tmp(k)
     &                  -35.d0*(.802d0-.002d0*tmp(k))
   10 continue
c
  200 continue
c
c     -- vertical eddy diffusivity --
c
	  call khcal
c
c     -- source terms --
c
      call pro
c
c     -- water temperature --
c
	  ww=0.d0
	  call adcal(ww,tmp,tmpn,qtmp,fkh)
c
c     -- state variables in the ecosystem --
c
c	if(nsw(6).eq.1) then
	  if(nsw(6).eq.1.and.mod(nt,intbal).eq.0) then
c
       do 20 m=1,np
        ww=wphy(m)
        call adcalf(m,ww,phy,phyn,qphy,fkh)
   20 continue
c
       do 22 m=1,nzp
        ww=0.d0
        call adcalf(m,ww,zoo,zoon,qzoo,fkh)
   22 continue
c
c       ww=0.d0
c       call adcal(ww,bac,bacn,qbac,fkh)
       ww=wpoc
       call adcal(ww,poc,pocn,qpoc,fkh)
       ww=0.d0
       call adcal(ww,doc,docn,qdoc,fkh)
       call adcal(ww,dip,dipn,qdip,fkh)
       call adcal(ww,dinh,dinhn,qdinh,fkh)
c      call adcal(ww,dis,disn,qdis,fkh)
       call adcal(ww,dox,doxn,qdox,fkh)
       call adcalb(ww,hb,hbn,qhb)
       call adcalb(ww,poly,polyn,qpoly)
       call adcalb(ww,other,othern,qother)
       call adcalp(ww,spfish,spfishn,qspfish)
       call adcalp(ww,sdfish,sdfishn,qsdfish)
       call adcalp(ww,pfish,pfishn,qpfish)
c
	  end if
c
c     -- update water temperature and density --
c
      do 30 k=0,nz
       tmp(k)=tmpn(k)
       dens(k)=1028.14d0-.0735d0*tmp(k)-.00469d0*tmp(k)*tmp(k)
     &                  -35.d0*(.802d0-.002d0*tmp(k))
   30 continue
c
c     -- update state variables in the ecosystem --
c
c	  if(nsw(6).eq.1) then
	  if(nsw(6).eq.1.and.mod(nt,intbal).eq.0) then
c
       do 32 k=1,nz
        do 34 m=1,np
         phy(m,k)=phyn(m,k)
         if(phy(m,k).le.1.d-8) phy(m,k)=1.d-8
   34 continue
       do 36 m=1,nzp
        zoo(m,k)=zoon(m,k)
        if(zoo(m,k).le.1.d-8) zoo(m,k)=1.d-8
   36 continue
c	  bac(k)=bacn(k)
        poc(k)=pocn(k)
        doc(k)=docn(k)
        dip(k)=dipn(k)
        dinh(k)=dinhn(k)
c        dis(k)=disn(k)
        dox(k)=doxn(k)
        hb(k)=hbn(k)
        poly(k)=polyn(k)
        other(k)=othern(k)
        spfish(k)=spfishn(k)
        sdfish(k)=sdfishn(k)
        pfish(k)=pfishn(k)
c   
c        if(bac(k).le.1.d-8) bac(i,j,k)=1.d-8
        if(poc(k).le.1.d-8) poc(k)=1.d-8
        if(doc(k).le.1.d-8) doc(k)=1.d-8
        if(dip(k).le.1.d-8) dip(k)=1.d-8
        if(dinh(k).le.1.d-8) dinh(k)=1.d-8
c        if(dis(k).le.1.d-8) dis(k)=1.d-8
        if(dox(k).le.1.d-8) dox(k)=1.d-8
        if(hb(k).le.1.d-8) hb(k)=1.d-8
        if(poly(k).le.1.d-8) poly(k)=1.d-8
        if(other(k).le.1.d-8) other(k)=1.d-8
        if(spfish(k).le.1.d-8) spfish(k)=1.d-8
        if(sdfish(k).le.1.d-8) sdfish(k)=1.d-8
        if(pfish(k).le.1.d-8) pfish(k)=1.d-8
   32 continue
c
	  end if
c
c     -- vertical mixing --
c
      if(nsw(7).ne.2) call mixcal
c
      return
      end
c-----------------------------------------------------------------------
c
c     -- Vertical Eddy Diffusivity Coefficients --
c                                               Author: Daisuke Kitazawa
c                                               Last Update: 2005.7.29
c
c-----------------------------------------------------------------------
      subroutine khcal
      implicit double precision (a-h,o-z)
      include 'param.h'
c
c     -- Kh is constant --
c
	  if(nsw(7).eq.0) then
c
	  if(nsw(2).eq.1) then
	   do 10 k=0,nz
	    fkh(k)=fkh0
   10 continue
      end if
c
c	-- Kh depends on stratification function --
c
	  elseif(nsw(7).eq.1) then
       if(nsw(2).eq.1) then
        do 20 k=1,nz-1
         rhodz=(dens(k+1)-dens(k))/dzh(k)
	   if(rhodz.lt.0.d0) rhodz=0.d0
        rhom=(dens(k)*ddz(k+1)+dens(k+1)*ddz(k))*.5d0/dzh(k)
        udz=(ui(k+1)-ui(k))*.5d0/dzh(k)
        vdz=(vj(k+1)-vj(k))*.5d0/dzh(k)
        uvmdz2=udz*udz+vdz*vdz
        epsi=btkc*g*rhodz/(rhom*((fkhmin/fkh0)**(1.d0/alkc)-1.d0))
c
       if(rhodz.gt.1.0d-10) then
        if(uvmdz2.le.epsi) then
         fkh(k)=fkhmin
        else
         fkh(k)=fkh0*(1.d0+btkc*g*rhodz/rhom/uvmdz2)**alkc
        end if
       end if
c
       if(fkh(k).lt.fkhmin) fkh(k)=fkhmin
c
   20 continue
c
	   end if
	  end if
c
      return
      end
c-----------------------------------------------------------------------
c
c     -- Production Terms --
c                                               Author: Daisuke Kitazawa
c                                               Last Update: 2005.7.29
c
c-----------------------------------------------------------------------
      subroutine pro
      implicit double precision (a-h,o-z)
      include 'param.h'
      dimension b1(npmax),b2(npmax),b3(npmax),b4(npmax)
      dimension b6(npmax),b7(npmax),b8(npmax),b9(npmax)
      dimension d1(npmax),d2(npmax),d3(npmax)
c    エラー対策↓-----------param.h に追加-------------------------------------------------------*解決？
c      common /pol/qpoly(nzmax),poly(nzmax),bother(nzmax),qbother(nzmax)
c
      do 10 k=nz,1,-1
c
	  qtmp(k)=0.d0
c
c       if(nsw(6).eq.0) ec=0.2d0
       ec=0.2d0
c
c     -- bulk formula --
c
      if(nsw(4).eq.0) then 
c
       q1=qs0*(1.d0-.71d0*cloud)*(1.d0-ref)
c
c        q2=s_s*sigma*(tmpaq+273.15d0)**4.d0
c     &     *(.39d0-5.8d-2*sqrt(e_a))*(1.d0-clat*(cldq**2.d0))
c     &     +4.d0*s_s*sigma*(tmpaq+273.15d0)**3.d0
c     &     *(tmp(i,j,k)-tmpaq)
       q2=s_s*sigma*(tmp(k)+273.15d0)**4.d0
     &     *(.49d0-.066d0*sqrt(e_a))*(1.d0-clat*(cldq**2.d0))
     &     +4.d0*s_s*sigma*(tmp(k)+273.15d0)**3.d0
     &     *(tmp(k)-tmpaq)
c
       ww=sqrt(wx*wx+wy*wy)
       if(ww.ge.1.d0)then
         q3=f_i*densa*ce*(q_s-q_e)*ww
         q4=c_a*densa*ch*(tmp(k)-tmpa)*ww
        else
         tmpfc=(tmp(k)-tmpa)+.61d0*(tmpa+273.15d0)*(q_s-q_e)
       if(tmpfc.lt.0.d0)then
         q3=-f_i*densa*ce2*(q_s-q_e)*(abs(tmpfc))**(1.d0/3.d0)
         q4=-c_a*densa*ch2*(tmp(k)-tmpa)
        elseif(tmpfc.ge.0.d0)then
         q3=f_i*densa*ce2*(q_s-q_e)*tmpfc**(1.d0/3.d0)
         q4=c_a*densa*ch2*(tmp(k)-tmpa)*tmpfc**(1.d0/3.d0)
       end if
       end if
c
       elseif(nsw(4).eq.1) then
c
        e_s=6.1078d0*(10.d0**(7.5d0*tmpaq/(237.3d0+tmpaq)))
        e_a=e_s*humq
        q_s=.622d0*e_s/(presq-.378d0*e_s)
        q_e=.622d0*e_a/(presq-.378d0*e_a)
        densa=1.293d0*273.15d0*(presq-.378d0*e_a)
     &       /(273.15d0+tmpaq)/1013.25d0
c
        q1=sunnq*(1.d0-ref)
c
c        q2=s_s*sigma*(tmpaq+273.15d0)**4.d0
c     &     *(.39d0-5.8d-2*sqrt(e_a))*(1.d0-clat*(cldq**2.d0))
c     &     +4.d0*s_s*sigma*(tmpaq+273.15d0)**3.d0
c     &     *(tmp(i,j,k)-tmpaq)
        q2=s_s*sigma*(tmp(k)+273.15d0)**4.d0
     &     *(.49d0-.066d0*sqrt(e_a))*(1.d0-clat*(cldq**2.d0))
     &     +4.d0*s_s*sigma*(tmp(k)+273.15d0)**3.d0
     &     *(tmp(k)-tmpaq)
c
       ww=sqrt(wx*wx+wy*wy)
        if(ww.ge.1.d0)then
         q3=f_i*densa*ce*(q_s-q_e)*ww
         q4=c_a*densa*ch*(tmp(k)-tmpaq)*ww
        else
         tmpfc=(tmp(k)-tmpaq)
     &         +.61d0*(tmpaq+273.15d0)*(q_s-q_e)
        if(tmpfc.lt.0.d0)then
         q3=-f_i*densa*ce2*(q_s-q_e)*(abs(tmpfc))**(1.d0/3.d0)
         q4=-c_a*densa*ch2*(tmp(k)-tmpaq)
     &       *(abs(tmpfc))**(1.d0/3.d0)
         elseif(tmpfc.ge.0.d0)then
          q3=f_i*densa*ce2*(q_s-q_e)*tmpfc**(1.d0/3.d0)
          q4=c_a*densa*ch2*(tmp(k)-tmpaq)
     &           *tmpfc**(1.d0/3.d0)
        end if
       end if
	  end if
c
      if(k.eq.1) then
       fu=q1*(exp(-ec*zb(k-1))-exp(-ec*zb(k)))/(dens0*c_0)
     &    -(q2+q3+q4)/(dens0*c_0)
       else
      if(k.eq.nz) then
       fu=q1*(exp(-ec*zb(k-1))-0.d0)/(dens0*c_0)
       else
       fu=q1*(exp(-ec*zb(k-1))-exp(-ec*zb(k)))/(dens0*c_0)
      end if
      end if
c
       qtmp(k)=fu/ddz(k)
c
c     -- river --
c
      if(k.eq.1) then
       do 20 n=1,nriv
        qtmp(k)=qtmp(k)+qriv(n)*(triv(n)-tmp(k))
     &                 /(ddx*ddy*ddz(k)+qriv(n))
       goto 22
   20 continue
      end if
   22 continue
c
c	  if(nsw(6).eq.1) then
	  if(nsw(6).eq.1.and.mod(nt,intbal).eq.0) then
c
c     -- extinction coefficient --
c
       physum=0.d0
       phys1=0.d0
       chlaave=0.d0
        do 30 m=1,np
         physum=physum+phy(m,k)
         phys1=phys1+phy(m,1)
         chlaave=chlaave+chlac(m)
   30 continue
       chlaave=chlaave/np
       ec=ec0+ec1*chlaave*phys1
c
c     -- phytoplankton --
c
       do 40 m=1,np
        b1(m)=gp(m)*thep(m)**(tmp(k)-20.d0)
     &        *dmin1(dip(k)/(dip(k)+hsp(m)),
     &               dinh(k)/(dinh(k)+hsn(m)))
     &        *(dexp(1.d0-q1*dexp(-ec*zb(k))/aip(m))
     &         -dexp(1.d0-q1*dexp(-ec*zb(k-1))/aip(m)))
     &         /ec*ddz(k)*phy(m,k)
	  if(k.eq.nz) then
       b1(m)=gp(m)*thep(m)**(tmp(k)-20.d0)
     &       *dmin1(dip(k)/(dip(k)+hsp(m)),
     &        dinh(k)/(dinh(k)+hsn(m)))
     &        *(dexp(1.d0-q1*0.d0/aip(m))
     &       -dexp(1.d0-q1*dexp(-ec*zb(k-1))/aip(m)))
     &          /ec*ddz(k)*phy(m,k)
	  end if
       b2(m)=rp(m)*thep(m)**(tmp(k)-20.d0)*phy(m,k)
       b3(m)=erp(m)*dexp(gamp(m)*chlac(m)*phy(m,k))*b1(m)
       b4(m)=dmp(m)*phy(m,k)*phy(m,k)
   40 continue
c
c     -- zooplankton --
c   
      zoosum=0.d0
       do m2=1,nzp
        zoosum=zoosum+zoo(m,k)
       end do 
       do 42 m=1,nzp
c         b6(m)=az(m)*cz(m)*hsz(m)*(physum+poc(i,j,k))
c     &        /(hsz(m)+physum+poc(i,j,k))*zoo(m,i,j,k)
        b6(m)=cz(m)*thez(m)**(tmp(k)-20.d0)*(1.d0-
     &        dexp(.007d0*(hsz(m)-physum-poc(k))))*zoo(m,k)
        b7(m)=rz(m)*thez(m)**(tmp(k)-20.d0)*zoo(m,k)
        b8(m)=(1.d0-az(m))*b6(m)
        b9(m)=dmz(m)*zoo(m,k)*zoo(m,k)
   42 continue
c
c     -- particulate organic carbon --
c
       b10=rpoc*thepoc**(tmp(k)-20.d0)*poc(k)
c         c11=rpoc*thepoc
c         c12=(tmp(k)-20.d0)*poc(k)
       b11=ft*b13
c
c     -- dissolved organic carbon --
c
       b13=rdoc*thedoc**(tmp(k)-20.d0)*doc(k)
c
c     -- dissolved oxygen --
c
      if(k.eq.1) then
        dosu=14.161d0-.3943d0*tmp(k)+7.714d-3*tmp(k)**2.d0
     &       -6.46d-5*tmp(k)**3.d0
     &       -0.d0*(.1519d0-4.62d-3*tmp(k)+6.76d-5*tmp(k)**2.d0)
       b17=rea*(dosu-dox(k))/ddz(k)
       else
       b17=0.d0
      end if
c   -- hervibous benthos -- 
c feeding = (P/Q)*(search rate*vulnerability*B prey*B predator/2*vul+search rate*B predator)
        pq4hb=.101d0
        hbsr2phy=4.228d-9
        vulhb2phy=1.d0
       b18=pq4hb*(hbsr2phy*vulhb2phy*phy(1,k)*hb(k)
     &    /(2.d0*vulhb2phy+hbsr2phy*hb(k)))
c
c  other mortality = (1-EE)*PB*B
        ee4hb=.383d0
        pb4hb=4.499d-8
       b20=(1-ee4hb)*pb4hb*hb(k)
c
c   -- polychaeta -- 
c         
c  feeding = (P/Q)*(search rate*vulnerability*B prey*B predator/2*vul+search rate*B predator)
        pq4pol=.205d0
        polsr2det=7.045d-10
        vulpol2det=1.03d0
       b31=pq4pol*polsr2det*vulpol2det*poc(k)*poly(k)
     &    /(2.d0*vulpol2det+polsr2det*poly(k))
c
c  other mortality = (1-EE)*PB*B
      ee4ply=.01d0
      pb4ply=2.708d-09
c     non-predation mortality
      nonpm=1.d0-ee4ply
      b32=nonpm*pb4ply*poly(k)
c
c   -- other benthos -- 
c  feeding = (P/Q)*(search rate*vulnerability*B prey*B predator/2*vul+search rate*B predator)
c   other benthos VS detritus         
          pq4oth=.185d0
          othsr2det=2.818d-10
          vuloth2det=1.76d0
         bofd1=pq4oth*othsr2det*vuloth2det*poc(k)*other(k)
         bofd3=2.d0*vuloth2det+othsr2det*other(k)
        obvsdet=bofd1/bofd3
c   other benthos VS  zooplankton         
          othsr2zp=1.114d-09
          vuloth2zp=1.9d+12
         bofz1=pq4oth*othsr2zp*vuloth2zp*zoo(1,k)*other(k)
         bofz2=2.d0*vuloth2zp+othsr2zp*other(k)
        obvszp=bofz1/bofz2
c       
       b41=obvsdet+obvszp
c
c   other mortality = (1-EE)*PB*B
          ee4oth=.125d0
          pb4oth=1.957d-09
       b42=(1.d0-ee4oth)* pb4oth * other(k)
c
c   -- small pelagic fish -- 
c  feeding = (P/Q)*(search rate*vulnerability*B prey*B predator/2*vul+search rate*B predator)
c   small pelagic fish VS phtoplankton        
          pq4spf=.113d0
          spfsr2phy=1.818d-09
          vulspf2phy=3.41d0
         spffphy1=pq4spf*spfsr2phy*vulspf2phy
     &          *phy(1,k)*spfish(k)
         spffphy2=(2.d0*vulspf2phy+spfsr2phy)*spfish(k)
        spfvsphy=spffphy1/spffphy2
c   small pelagic fish VS  zooplankton  
          spfsr2zp=7.183d-10
          vulspf2zp=1.d0
         spffz1=pq4spf*spfsr2zp*vulspf2zp*zoo(1,k)*spfish(k)
         spffz2=(2.d0*vulspf2zp+spfsr2zp)*spfish(k)
        spfvszp=spffz1/spffz2
       b51=spfvsphy+spfvszp
c
c  other mortality = (1-EE)*PB*B
          ee4spf=.294d0
          pb4spf=1.283d-09
       b52=(1.d0-ee4spf)*pb4spf*spfish(k)
c
c   -- small demersal fish -- 
c  feeding = (P/Q)*(search rate*vulnerability*B prey*B predator/2*vul+search rate*B predator)
c  small demersal fish VS   hervibous benthos      
        pq4sdf=.106d0
        sdfsr2hb=1.756d-10
        vulsdf2hb=1.d0
       sdffhb1=pq4sdf*sdfsr2hb*vulsdf2hb*hb(k)*sdfish(k)
       sdffhb2=2.d0*vulsdf2hb+sdfsr2hb*sdfish(k)
      sdfvshb=sdffhb1/sdffhb2
c     
c  small demersal fish VS   polychaeta      
        pq4sdf=.106d0
        sdfsr2ply=2.026d-11
        vulsdf2ply=2.06d0
       sdffply1=pq4sdf*sdfsr2ply*vulsdf2ply*poly(k)*sdfish(k)
       sdffply2=2.d0*vulsdf2poly+sdfsr2ply*sdfish(k)
      sdfvsply=sdffply1/sdffply2
c
c  small demersal fish VS   other benthos      
        pq4sdf=.106d0
        sdfsr2other=2.d-10
        vulsdf2other=3.76d+1
       sdffoth1=pq4sdf*sdfsr2other*vulsdf2other
     &        *other(k)*sdfish(k)
       sdffoth2=2.d0*vulsdf2other+sdfsr2other*sdfish(k)
      sdfvsoth=sdffoth1/sdffoth2
c        
c   small demersal fish VS   detritus      
        pq4sdf=.106d0
        sdfsr2det=2.37d-10
        vulsdfanddet=1.02d0
       sdffdet1=pq4sdf*sdfsr2det*vulsdfanddet*poc(k)*sdfish(k)
       sdffdet3=2.d0*vulsdfanddet+sdfsr2det*sdfish(k)
      sdfvsdet=sdffdet1/sdffdet3
c
      b61=sdfvshb+sdfvsply+sdfvsoth+sdfvsdet 
c
c  other mortality = (1-EE)*PB*B
        ee4sdf=.248d0
        pb4sdf=1.030d-09
       npm4sdf=1.d0-ee4sdf
      b62=npm4sdf*pb4sdf*sdfish(k)
      write (*,*)  b61,b62
c
c   -- piscivorous fish -- 
c  feeding = (P/Q)*(search rate*vulnerability*B prey*B predator/2*vul+search rate*B predator)
c  piscivorous fish VS  hervibous benthos      
        pq4pfish=.091d0
        pfsr2hb=1.065d-10
        vulpf2hb=1.d0
       pffhb1=pq4pfish*pfsr2hb*vulpf2hb*hb(k)*pfish(k)
       pffhb3=2.d0*vulpf2hb+pfsr2hb*pfish(k)
      pfishvshb=pffhb1/pffhb3
c  piscivorous fish VS   polychaeta      
       pq4pfish=.091d0
       pfishsr2ply=.178d-10
       vulpf2ply=.109d+1
      pffply1=pq4pfish*pfishsr2ply*vulpf2ply*poly(k)*pfish(k)
      pffply3=2.d0*vulpf2ply+pfishsr2ply*pfish(k)
c     
c      write (*,*)  pffply3
      pfishvspoly=pffply1/pffply3
c  piscivorous fish VS   other benthos      
        pq4pfish=.091d0
        pfishsr2oth=1.21d-10
        vulpfish2oth=1.01d0
       pffoth1=pq4pfish*pfishsr2oth*vulpfish2oth*other(k)*pfish(k)
       pffoth3=2.d0*vulpfish2oth+pfishsr2oth*pfish(k)
      pfishvsoth=pffoth1/pffoth3
c  piscivorous fis h VS  small pelagic fish      
        pq4pfish=.091d0
        pfishsr2spf=1.31d-09
        vulpfish2spf=1.65d0
       pffspf1=pq4pfish*pfishsr2spf*vulpfish2spf*spfish(k)*pfish(k)
       pffspf3=2.d0*vulpfish2spf+pfishsr2spf*pfish(k)
      pfishvsspf=pffspf1/pffspf3
c  piscivorous fish VS  small demersal fish      
        pq4pfish=.091d0
        pfishsr2sdf=6.00d-10
        vulpfish2sdf=1.15d0
       pffsdf1=pq4pfish*pfishsr2sdf*vulpfish2sdf*sdfish(k)*pfish(k)
       pffsdf2=2.d0*vulpfish2sdf+pfishsr2sdf*pfish(k)
      pfishvssdf=pffsdf1/pffsdf2
c  piscivorous fish VS  zooplankton      
        pq4pfish=.091d0
        pfishsr2zp=1.08d-10
        vulpfish2zp=1.7d0
       pffzp1=pq4pfish*pfishsr2zp*vulpfish2zp*zoo(1,k)*pfish(k)
       pffzp3=2.d0*vulpfish2zp+pfishsr2zp*pfish(k)
      pfishvszp=pffzp1/pffzp3
c          
       b71=pfishvshb+pfishvspoly+pfishvsoth
     &     +pfishvsspf+pfishvssdf
c
c  other mortality = (1-EE)*PB*B
        ee4pfish=0.373d0
        pb4pfish=5.798d-10
       morpfish=(1.d0-ee4pfish)*pb4pfish
      b72=morpfish*pfish(k)
c
c     -- river --
c
       if(k.eq.1) then
        if(nriv.ge.1) then
         do 50 n=1,nriv
          qph=qriv(n)*0.d0/ddx/ddy/ddz(k)
          qzo=qriv(n)*0.d0/ddx/ddy/ddz(k)
          qpc=qriv(n)*(pcriv(n)-poc(k))/ddx/ddy/ddz(k)
          qdc=qriv(n)*(dcriv(n)-doc(k))/ddx/ddy/ddz(k)
          qdp=qriv(n)*(dpriv(n)-dip(k))/ddx/ddy/ddz(k)
          qdnh=qriv(n)*(dnhriv(n)-dinh(k))/ddx/ddy/ddz(k)
          qdo=qriv(n)*(dxriv(n)-dox(k))/ddx/ddy/ddz(k)
         goto 52
   50 continue
        end if
       else
          qph=0.d0
          qzo=0.d0
          qpc=0.d0
          qdc=0.d0
          qdp=0.d0
          qdnh=0.d0
          qds=0.d0
          qdo=0.d0
         end if 
   52 continue
c
c     -- prediction of state variables --
c
       b1psum=0.d0
       b1nsum=0.d0
       b1ssum=0.d0
       b2psum=0.d0
       b2nsum=0.d0
       b2ssum=0.d0
       b3sum=0.d0
       b4sum=0.d0
       b6sum=0.d0
       b7psum=0.d0
       b7nsum=0.d0
       b7ssum=0.d0
       b8sum=0.d0
       b9sum=0.d0
       d1sum=0.d0
       d2sum=0.d0
       d3sum=0.d0
c
      do 60 m=1,np
       b1psum=b1psum+dipph(m)*b1(m)
       b1nsum=b1nsum+dinph(m)*b1(m)
       b1ssum=b1ssum+disph(m)*b1(m)
       b2psum=b2psum+dipph(m)*b2(m)
       b2nsum=b2nsum+dinph(m)*b2(m)
       b2ssum=b2ssum+disph(m)*b2(m)
       b3sum=b3sum+b3(m)
       b4sum=b4sum+b4(m)
   60 continue
       do 62 m=1,nzp
        b6sum=b6sum+b6(m) 
        b7psum=b7psum+dipzo(m)*b7(m) 
        b7nsum=b7nsum+dinzo(m)*b7(m) 
        b7ssum=b7ssum+diszo(m)*b7(m) 
        b8sum=b8sum+b8(m) 
        b9sum=b9sum+b9(m) 
   62 continue
       do 64 m=1,np
        d1sum=d1sum+todph(m)*b1(m)
        d2sum=d2sum+todph(m)*b2(m)
   64 continue
       do 66 m=1,nzp
        d3sum=d3sum+todzo(m)*b7(m)
   66 continue
c
       do 70 m=1,np
        qphy(m,k)=b1(m)-b2(m)-b3(m)-b4(m)
     &            -b6sum*phy(m,k)/(physum+poc(k))
     &            -b18-spfvsphy+qph
   70 continue
       do 72 m=1,nzp
        qzoo(m,k)=b6(m)-b7(m)-b8(m)-b9(m)+qzo
     &         -obvszp-spfvszp-pfishvszp  
   72 continue
        qpoc(k)=b4sum+b8sum+b9sum+qpc
     &         -(b6sum*poc(k))/(physum+poc(k))
     &         -b10-b11-b31-obvsdet-sdfvsdet
c        pocpos=b4sum+b8sum+b9sum+qpc
c        pocneg=b6sum*(poc(k)/(physum+poc(k)))
c     &        +b10+b11+b31+obvsdet+sdfvsdet
c        qpoc(k)=pocpos-pocneg
c        
c        write (*,*)  qpoc(k)
        qdoc(k)=b3sum+b11-b13+qdc
c
        dipzod=(dipph(1)*(b7(1)+b8(1)+b9(1))-dipzo(1)*b9(1))
     &        /(b7(1)+b8(1))
        dinzod=(dinph(1)*(b7(1)+b8(1)+b9(1))-dinzo(1)*b9(1))
     &        /(b7(1)+b8(1))
        b8psud=dipzod*b8(1)
        b8nsud=dinzod*b8(1)
c
        dppomd=(dipph(1)*b4(1)+dipzod*b8(1)+dipzo(1)*b9(1))
     &        /(b4(1)+b8(1)+b9(1))
	    dnpomd=(dinph(1)*b4(1)+dinzod*b8(1)+dinzo(1)*b9(1))
     &        /(b4(1)+b8(1)+b9(1))
c
        dpdomd=(dipph(1)*b3(1)+dppomd*b11)
     &        /(b3(1)+b11)
	    dndomd=(dinph(1)*b3(1)+dnpomd*b11)
     &        /(b3(1)+b11)
c
        qdip(k)=-b1psum+b2psum+b7psum+dppomd*b10+dpdomd*b13+qdp
        qdinh(k)=-b1nsum+b2nsum+b7nsum+dnpomd*b10+dndomd*b13+qdnh
        qdox(k)=d1sum-d2sum-d3sum-topom*b10-todom*b13+b17+qdo
c
        qhb(k)=b18-b20-sdfvshb-pfishvshb
        qpoly(k)=b31-b32-sdfvsply-pfishvspoly
c         qpolyneg=b32+sdfvsply+pfishvspoly
        qother(k)=b41-b42-sdfvsoth-pfishvsoth
        qspfish(k)=b51-b52-pfishvsspf
        qsdfish(k)=0
c        b61-b62-pfishvssdf
        qpfish(k)=b71-b72
c       
c
c     -- exchange of nutrients and oxygen --
c
       if(k.eq.nz) then
c
        wphsum=0.
         do 78 m=1,np
        wphsum=wphsum+wphy(m)*phy(m,k)
   78 continue
        wpcsum=wpoc*poc(k)
c
        qdip(k)=qdip(k)+.0d0*(dipph(1)*wphsum+dppomd*wpcsum)/ddz(k)
        qdinh(k)=qdinh(k)+.4d0*(dinph(1)*wphsum+dnpomd*wpcsum)/ddz(k)
	    qdox(k)=qdox(k)-topom*(wphsum+wpcsum)/ddz(k)
c
       end if
c
	  end if
c
   10 continue
c
      return
      end
c-----------------------------------------------------------------------
c
c     -- Vertical Mixing under Unstable Condition --
c                                               Author: Daisuke Kitazawa
c                                               Last Update: 2005.7.29
c
c-----------------------------------------------------------------------
      subroutine mixcal
      implicit double precision (a-h,o-z)
      include 'param.h'
c
c     -- mean values of physical properties and state variables --
c
      do 10 k=2,nz-1
c
       if(dens(k-1).ge.dens(k)) then
c
        dzu=ddz(k-1)
        dzl=ddz(k)
c
        tmpm=(tmp(k)*dzl+tmp(k-1)*dzu)/(dzu+dzl)
        tmp(k)=tmpm
        tmp(k-1)=tmpm
c
        densm=1028.14-.0735d0*tmpm-.00469*tmpm*tmpm
     &               -35.*(.802-0.002*tmpm)
        dens(k)=densm
        dens(k-1)=densm
c
        do 12 m=1,np
         phym=(phy(m,k)*dzl+phy(m,k-1)*dzu)/(dzu+dzl)
         phy(m,k)=phym
         phy(m,k-1)=phym
   12 continue
c
        do 14 m=1,nzp
         zoom=(zoo(m,k)*dzl+zoo(m,k-1)*dzu)/(dzu+dzl)
         zoo(m,k)=zoom
         zoo(m,k-1)=zoom
   14 continue
c
c        bacm=(bac(k)*dzl+bac(k-1)*dzu)/(dzu+dzl)
c        bac(k)=bacm
c        bac(k-1)=bacm
c
        pocm=(poc(k)*dzl+poc(k-1)*dzu)/(dzu+dzl)
        poc(k)=pocm
        poc(k-1)=pocm
c
        docm=(doc(k)*dzl+doc(k-1)*dzu)/(dzu+dzl)
        doc(k)=docm
        doc(k-1)=docm
c
        dipm=(dip(k)*dzl+dip(k-1)*dzu)/(dzu+dzl)
        dip(k)=dipm
        dip(k-1)=dipm
c
        dinhm=(dinh(k)*dzl+dinh(k-1)*dzu)/(dzu+dzl)
        dinh(k)=dinhm
        dinh(k-1)=dinhm
c
c        dism=(dis(k)*dzl+dis(k-1)*dzu)/(dzu+dzl)
c        dis(k)=dism
c        dis(k-1)=dism
c
        doxm=(dox(k)*dzl+dox(k-1)*dzu)/(dzu+dzl)
        dox(k)=doxm
        dox(k-1)=doxm
c
c        以下cut(20160408)
c        hbm=(hb(k)*dzl+hb(k-1)*dzu)/(dzu+dzl)
c        hb(k)=hbm
c        hb(k-1)=hbm
c        
c        polym=(poly(k)*dzl+poly(k-1)*dzu)/(dzu+dzl)
c        poly(k)=polym
c        poly(k-1)=polym
c        ここまで
       end if
c
   10 continue
c
      return
      end

;calculate velocity fields for fernando's data
;
;USES:      .compile read_fernandos_data.pro
;MKD 2017, June 6
;
;Example:
;           
;           .r calce4klaus

;pro calce4klaus

    ;read in the data
    read_fernandos_data,bx,by,bz,vx,vy,vz,tim,xx,yy
    nt=n_elements(bx[0,0,*])
    
    ds=yy[1]-yy[0]
    carr=dblarr(nt-1)
    
    ;time step
    it=1  
    wpil=5.0
    dsx=xx-shift(xx,1)
    dsx=dsx[1:*]
    dsy=yy-shift(yy,1)
    dsy=dsy[1:*]
    mx=moment(dsx,sdev=sdevx)
    my=moment(dsy,sdev=sdevy)
    print,'X grid, mean and stdev [km] ',mx[0],sdevx
    print,'Y grid, mean and stdev [km] ',my[0],sdevy
    print,'X grid, mean and stdev [Mm] ',mx[0]/1d3,sdevx/1d3
    print,'Y grid, mean and stdev [Mm] ',my[0]/1d3,sdevy/1d3
    whitewin
    !p.multi=[0,2,2]    
    window,0,xs=800,ys=800,tit='Inductivity and E-field quality: CLEM'
    
    for i=0,nt-it-1 do begin
        i=15
        dt=tim[i+it]-tim[i] 
        dbx=(bx[*,*,i+it]-bx[*,*,i])/dt
        dby=(by[*,*,i+it]-by[*,*,i])/dt
        dbz=(bz[*,*,i+it]-bz[*,*,i])/dt
        dt=tim[i+it]-tim[i]
        vx_hlf=(vx[*,*,i+it]+vx[*,*,i])/2.
        vy_hlf=(vy[*,*,i+it]+vy[*,*,i])/2.
        vz_hlf=(vz[*,*,i+it]+vz[*,*,i])/2.
        bx_hlf=(bx[*,*,i+it]+bx[*,*,i])/2.
        by_hlf=(by[*,*,i+it]+by[*,*,i])/2.
        bz_hlf=(bz[*,*,i+it]+bz[*,*,i])/2.
        
        ;SOLVING VANILLA PTD EQUATIONS
        ptdsolve,dbx,dby,dbz,ds,dt,scriptbt,dscriptbdzt,scriptjt

        ;CALCULATE ELECTRIC FIELD VECTOR
        scriptb2evec,scriptbt,dscriptbdzt,scriptjt,ds,eptd
        
        eflct=e_flct(bx_hlf,by_hlf,bz_hlf,vx_hlf,vy_hlf,ds,sigma=wpil)
        edopp=e_doppler(bx_hlf,by_hlf,bz_hlf,vz_hlf,ds,sigma=wpil)
        bvec=dblarr(n_elements(bx(*,0,0)),n_elements(bx(0,*,0)),3)
        bvec(*,*,0)=bx_hlf
        bvec(*,*,1)=by_hlf
        bvec(*,*,2)=bz_hlf

        mk_relax_psi_3d,bvec,eptd+edopp+eflct, psi, dpsi_dz,dx=ds,max_iter=20,etot=etot

        ;CALCULATE ACTUAL ELECTRIC FIELD VECTOR EX EY EZ
        ;E-actual at half time: t[i+1/2]=(t[i+1]+t[i])/2
        axb,vx_hlf,vy_hlf,vz_hlf,bx_hlf,by_hlf,bz_hlf,$
            axbx_hlf,axby_hlf,axbz_hlf
        ex_hlf=-axbx_hlf
        ey_hlf=-axby_hlf
        ez_hlf=-axbz_hlf

        ;E-actual at t=t[i]
         axb,vx[*,*,i],vy[*,*,i],vz[*,*,i],$
             bx[*,*,i],by[*,*,i],bz[*,*,i],$
             axbx_i,axby_i,axbz_i
         ex_i=-axbx_i
         ey_i=-axby_i
         ez_i=-axbz_i
        
         sz_hlf=-10.^5*(ey_hlf*bx_hlf-ex_hlf*by_hlf)/(4*!Pi)
         sz=-10.^5*(etot(*,*,1)*bx_hlf-etot(*,*,0)*by_hlf)/(4*!Pi)

         ;Test inductivity of the simulation by comparing curlh(E-actual_hlf) with dbz/dt
         ind_test,bz[*,*,i],bz[*,*,i+it],ex_hlf,ey_hlf,ds,dt,curl=curl,dbz=dbz,/quiet

         ;ind_test,bz[*,*,i],bz[*,*,i+it],eptd[*,*,0],eptd[*,*,1],ds,dt,curl=curl,dbz=dbz,/quiet
         
         carr[i]=correlate(curl,-dbz)      
         xyplot,ex_hlf,etot[*,*,0],tit='Ex'
         xyplot,ey_hlf,etot[*,*,1],tit='Ey'
         xyplot,ez_hlf,etot[*,*,2],tit='Ez'
        ;if i eq 10 then stop 
        stop
    endfor

    result=moment(carr)
    
    print, 'Mean of corr. b/w dbz(t)/dt and curlh(E(t)) : ', result[0] & PRINT, 'Variance: ', result[1] & $
    PRINT, 'Skewness: ', result[2] & PRINT, 'Kurtosis: ', result[3]
    whitewin
    
end



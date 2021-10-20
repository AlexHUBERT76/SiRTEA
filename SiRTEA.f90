program Sirtea
 
implicit none

    real*8 :: t, ti, dz, zref, Mexo, Rexo,Tstar,Rstar,Rexostar, g, H, mu, l, Te, dp, nref, Ratm
    real*8 :: zp, sp, ap, ai, da, ds, D, ray,api,Lumi
    integer :: j, i, nloop,iray,imin
    integer, parameter :: nc = 1000 !nombre de couches atmosphérique modélisés
    real*8, parameter :: Gr = 6.67408E-11, Kb = 1.38064852E-23, pi = 3.141592653589793238462643383279 ,mh = 1.6737236E-27!constantes utiles
    real*8, parameter :: Mt = 5.972E24, Rt = 6371000.0, UA = 149597870700.0 ! Masse et Rayon terrestre, distance à Terre Soleil
    real*8, parameter :: Rsol = 6.957E8 !Rayon solaire
    real*8, dimension (nc+1) :: n


    real*8 ::p0,nol,Lum,planck,p,k,nz,raymin,raymax,raypas,nolmin,nolmax



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!data!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




    real*8 , allocatable  :: tab(:),tq(:),tMM(:),tC(:),tn0(:),tS0(:),tA(:),tgamair(:),tgamself(:),tEn(:),tnair(:),tshift(:)
    integer, allocatable :: tmolid(:)
    integer  :: imax 
    integer :: molid 
    real*8  :: no0,S0,A,gamair,gamself,E,expT,shift,S0moy


    real*8 , dimension (200,6) :: Mmc 
    integer :: gj
    real*8  :: ab,Q,MM,C
    
    
    logical,allocatable :: mask(:),maskmed(:)
    integer :: itrue,itruemed
    real*8,allocatable,dimension(:) :: tab_short, tq_short,tC_short,tn0_short,tS0_short
    real*8,allocatable,dimension(:) :: tA_short,tgamair_short,tgamself_short,tnair_short
    real*8,allocatable,dimension(:) :: tab_med,tq_med,tC_med,tn0_med,tS0_med,tA_med
    real*8,allocatable,dimension(:) :: tgamair_med,tgamself_med,tnair_med
    logical :: in_atmo
    integer :: phase
    
    
    
    open(200,status = 'old',form = 'formatted',file = 'Parameter.dat')
    open(100,status = 'old',form = 'formatted',file = 'hitran.dat')
    open (90,file='spectre.dat', status = "replace")
    open (91,file='spectretaux.dat', status = "replace")
    open (92,file='spectrecorpnoir.dat', status = "replace")
    open (10,file='deviation.dat', status = "replace")
    open (42,file='chemin.dat', status = "replace")
    Mmc = 0
    i = 1

    10000 continue


    do 
        read(200,1000,err = 10000,end = 20000)molid,ab,Q,gj,MM,C

        Mmc(i,1) = molid
        Mmc(i,2) = ab
        Mmc(i,3) = q
        Mmc(i,4) = gj
        Mmc(i,5) = MM
        Mmc(i,6) = C


        i = i+1
    enddo

    20000 continue
    1000 format (i12,e13.5,e14.4,i5,f14.6,f11.5)




    imax = 10000000 !!!!!
    i = 0
    
    30000 continue
    do  !determination de la taille de la matrice des données
        read(100,*,end=300,err = 30000)

        i = i+1


    end do
    300 continue

    close (100)
    open(100,status = 'old',form = 'formatted',file = 'hitran.dat')



    imax = i

    allocate (tmolid(imax))
    allocate (tab(imax))
    allocate (tq(imax))
    allocate (tmm(imax))
    allocate (tC(imax))
    allocate (tn0(imax))
    allocate (tS0(imax))
    allocate (tA(imax))
    allocate (tgamair(imax))
    allocate (tgamself(imax))
    allocate (tEn(imax))
    allocate (tnair(imax))
    allocate (tshift(imax))
    allocate (maskmed(imax))


    i = 1

    50000 continue
    do while (i  <= imax) !chargement du ficher data dans la ram
        read(100,*,end = 200,err = 50000)molid,no0,S0,A,gamair,gamself,E,expT,shift

        do j = 1,200
            if (int(molid) == int(Mmc(j,1))) then !recoupage des information des deux plages de données




            tmolid(i) = int(Mmc(j,1))
            tab(i) = Mmc(j,2)
            tq(i) = Mmc(j,3)
            tmm(i) = Mmc(j,5)
            tC(i) = Mmc(j,6)
            tn0(i) = no0
            tS0(i) = S0
            tA(i) = A
            tgamair(i) = gamair
            tgamself(i) = gamself
            tEn(i) = E
            tnair(i) = expt
            tshift(i) = shift




            endif








        end do

        !900 format(i2,i1,f12.6,f10.3,e10.3,f5.4,f5.3,f10.4,f4.2,f8.6,a15,a15,a15,a15,6i1,6i2,a1,f7.1,f7.1)
        i = i +1
    enddo
    200 continue



    close (100)
    close (200)








    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!data!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




    Mexo = 1.0* Mt
    Rexo = 1.0* Rt
    Ratm = 1.05*Rexo
    !Masse et rayon du corp du systeme étudier
    
    Rstar = 1.0*Rsol
    Rexostar = 1.0*UA

    Tstar = 5778 !temperature de la surface de l'étoile d'ou provient les rayons en kelvin 

    dz = (Ratm-Rexo)/nc
    dp = dz/1000


    p0 = 1 !pression au niveau du sol (z = 0) en atm

    Te = 15 + 273.15 ! Température d'équilibre de l'atmoshpère (isothermie vertical et azimutal)
    call masse_reduite(mu) ! masse reduite des composants de l'atmosphere
    zref = 0

    g = Gr*Mexo/Rexo**2 !pesanteur à la surface. On considere g constante, approximation legitime pour une planete tellurique
    H = kb*Te/(mu*g*mh) !H est une echelle de pression, on consière que l'atmosphère est isotherme 

    raymin = 1000 !etendue de la longueur d'onde balayé en nm
    raymax = 2500
    raypas = 0.1
    
    
    nolmax = (1/raymin)*1E7
    nolmin = (1/raymax)*1E7


        !maskmed = .true.
        maskmed = ((((tn0 < nolmax + 20*(tgamair+tgamself)) .and.&
        ( tn0 > nolmin - 20*(tgamair+tgamself))) .or. (tS0 > 1000*S0moy)) .and. (tab > 1E-4)&
        .and. (tS0 /= 0.0) .and. (tC /= 0) ) !sensibilité de detection des raies
        itruemed = count(maskmed)
        allocate(tab_med(itruemed), tq_med(itruemed),tC_med(itruemed),tn0_med(itruemed),tS0_med(itruemed))
        allocate(tA_med(itruemed),tgamair_med(itruemed),tgamself_med(itruemed),tnair_med(itruemed),mask(itruemed))

        tab_med = pack(tab,maskmed)
        tq_med = pack(tq,maskmed)
        tC_med = pack(tc,maskmed)
        tn0_med = pack(tn0,maskmed)
        tS0_med = pack(ts0,maskmed)
        tA_med = pack(ta,maskmed)
        tgamair_med = pack(tgamair,maskmed)
        tgamself_med = pack(tgamself,maskmed)
        tnair_med = pack(tnair,maskmed)
        !print*, tab_med    


        S0moy = sum(tS0_med)/itruemed
    
    deallocate (tmolid)
    deallocate (tab)
    deallocate (tq)
    deallocate (tmm)
    deallocate (tC)
    deallocate (tn0)
    deallocate (tS0)
    deallocate (tA)
    deallocate (tgamair)
    deallocate (tgamself)
    deallocate (tEn)
    deallocate (tnair)
    deallocate (tshift)
    deallocate (maskmed) 
    nloop = int((raymax-raymin)/raypas)
    
    imin = 0
    api = 0
    zp = 0
    ti = 0
    t = 0
    sp = 0 
    ray = 0
    phase = 0
    p = 0
    nz = 0
    nref = 0
    nol = 0
    lum = 0
    l = 0
    k = 0
    itrue = 0
    in_atmo = .true.
    ds = 0
    da = 0
    D = 0
    ap = 0
    ai = 0
    Lumi = 0
   
    !parallélisation, num_treads à ajuster suivant la machine
    !$OMP PARALLEL DO SCHEDULE(STATIC) num_threads(8) DEFAULT(FIRSTPRIVATE) &
    !$OMP SHARED(tab_med,tq_med,tC_med,tn0_med,tS0_med,tA_med,tgamair_med,tgamself_med,S0moy) &
    !$OMP SHARED(tnair_med,tab,tq,tMM,tC,tn0,tS0,tA,tgamair,tgamself,tEn,tnair,tshift,tmolid)
    
    do iray = 0, nloop

        ray = raymin + iray * raypas






        t = 72.310 !angle d'impacte au zenith

        t = t*2*pi/360
        ti = t

        ap = 0 
        api = ap
        zp = Ratm

        l = ray

        nol = (1/l)*1E7 !nombre d'onde en cm-1


        Lum = planck(nol,Tstar)*(Rstar**2)/(Rexostar**2)
        Lumi = Lum
        

        mask = ((abs(nol-tn0_med) < 20*(tgamair_med+tgamself_med)) .or. (tS0_med > 1000*S0moy)) !sensibilité de detection des raies


        itrue = count(mask)
        allocate(tab_short(itrue), tq_short(itrue),tC_short(itrue),tn0_short(itrue),tS0_short(itrue))
        allocate(tA_short(itrue),tgamair_short(itrue),tgamself_short(itrue),tnair_short(itrue))

        tab_short = pack(tab_med,mask)
        tq_short = pack(tq_med,mask)
        tC_short = pack(tc_med,mask)
        tn0_short = pack(tn0_med,mask)
        tS0_short = pack(ts0_med,mask)
        tA_short = pack(ta_med,mask)
        tgamair_short = pack(tgamair_med,mask)
        tgamself_short = pack(tgamself_med,mask)
        tnair_short = pack(tnair_med,mask)
        

        l = l*1E-3

        nref = 0.05792105/(238.0185-l**(-2))+0.00167917/(57.362-l**(-2)) + 1 ! dispersion de l'air, 15°C, 1bar

        do j = 1, nc

            n(j) =(nref-1)*exp((zref-(j-1)*dz)/H)+1 !changement d'indice selon la couche atmoshérique

        enddo
        n(nc+1) = 1.0 !indice du vide
        n = 1.000






        i = nc
        imin  = i
        


        sp = 0
        ds = 0
        da = 0
        
        in_atmo = .TRUE.
        phase = 1
        
        do while (in_atmo ) !entrée dans l'atmosphère
        

        

            if  (zp <=  Rexo )  then

                t = t + da
                imin = 0
                D = ap+ti-t
                
                p = p0*exp(-(zp-rexo)/H)

                nz = (p*101315/(kb*Te))*1E-6 !molec/cm**3

                call k2_short(p,nz,te,tC_short,tab_short,tS0_short,tgamair_short,tgamself_short,tnair_short,tn0_short,nol,k,itrue)

                Lum = Lum*exp(-k*ds*100)+planck(nol,Te)*(1-(exp(-k*ds*100)))
    
                in_atmo = .FALSE.
                exit
            else if (zp > Ratm + dp) then !sortie de l'atmosphère


                D = t + ti + api -pi
                
                
                
                p = p0*exp(-(zp-rexo)/H)

                nz = (p*101315/(kb*Te))*1E-6 !molec/cm**3

                call k2_short(p,nz,te,tC_short,tab_short,tS0_short,tgamair_short,tgamself_short,tnair_short,tn0_short,nol,k,itrue)

                Lum = Lum*exp(-k*ds*100)+planck(nol,Te)*(1-(exp(-k*ds*100)))
                
                in_atmo = .FALSE.
                exit
             else if ((t+da >= pi/2) .and. (phase == 1)) then
             

                imin = i
                phase = 2
                continue
            endif
        
            select case(phase) !phase d'entrée
            case(1)
                if (zp <= (i-1)*dz + Rexo) then ! refraction sur une couche inférieur
                    
                    api = ap
                    
                    call descente(p0, zp, rexo,H,kb,Te, tC_short,tab_short,tS0_short,tgamair_short,tgamself_short,&
                                    tnair_short,tn0_short, nol, itrue, lum,ds,n,t,da,i,nc,dz)
               ! call xy(zp,ap)

                endif
                zp = zp - dp*cos(t+da)           
                ap = ap + dp*sin(t)/zp
                !le photon avance d'un pas dp
                da = da + dp*sin(t)/zp

                sp = sp + dp*n(i)

                ds = ds + dp


            case(2) !phase de transition
                



                zp = (i) * dz + Rexo
                sp = sp - ds*n(i)
                
                ap = api
            
                da = pi -2*t

                ap = ap + da
                
                api = ap
                
                ds = 2*zp*tan(da/2)
                

                call transition(p0, zp, rexo,H,kb,Te, tC_short,tab_short,tS0_short,tgamair_short,tgamself_short,&
                                tnair_short,tn0_short, nol, itrue, lum,ds)


                t = asin(n(i)/n(i+1)*sin(t)) !changement d'angle
                da = 0
                ds = 0

                i = i + 1  ! passage dans la couche supérieur
               ! call xy(zp,ap)


               
                
                phase = 3
            case(3) !phase de sortie
                if (zp >= (i)*dz + Rexo) then 
                
                    api = ap
                    
                    call montee(p0, zp, rexo,H,kb,Te, tC_short,tab_short,tS0_short,tgamair_short,tgamself_short,&
                                    tnair_short,tn0_short, nol, itrue, lum,ds,n,t,da,i,nc,dz)
              !  call xy(zp,ap)

                endif            
                zp = zp + dp*cos(t-da)
                ap = ap + dp*sin(t)/zp
                !le photon avance d'un pas dp
                da = da + dp*sin(t)/zp

                sp = sp + dp*n(i)

                ds = ds + dp



            

            end select

    

        enddo

        deallocate(tab_short, tq_short,tC_short,tn0_short,tS0_short,tA_short,tgamair_short,tgamself_short,tnair_short)



        D = D*360/(2*pi)
        t = t*360/(2*pi)
        ap = ap*360/(2*pi)

        !$OMP CRITICAL
        !print*,'ray=',ray,'nm'
        !print*, 'no=',nol,'cm-1'
        !print*, 'couche',i
        !print*, 'couche de transition',imin
        !print*, 'nmax=',n(imin)
        !print*, 'nref=',nref
        !print*, 'tf =',t, '° (Zenith parcouru)'
        !print*, 'apf =',ap, '° (azimut parcourue)'
        !print*, 'D =',D,'° (Angle de déviation)'
        !print*, 'sp =',sp,'m (chemin optique)'
        !print*, 'L =',Lum,'w.m-2.sr/cm-1 (luminance final)'
        !print*, 'fraction de luminance transmise =',Lum/Lumi,'sans unité'


       ! write (91,*) ray,D,imin
        !write (90,*) ray,lum
        write (91,*) ray,Lum/Lumi,imin,D,sp !!!!tri des données necessaire (parallélisation)
        !write (92,*) ray,planck(nol,tstar)
        !$OMP END critical
        
    enddo
    !$OMP END PARALLEL DO 
    
    deallocate (mask) 
    deallocate(tab_med,tq_med,tC_med,tn0_med,tS0_med,tA_med,tgamair_med,tgamself_med,tnair_med)
    
    close(10)
    close(90)
    close(91)
    close(92)
    close(42)
end program SiRTEA




subroutine montee(p0, zp, rexo,H,kb,Te, tC_short,tab_short,tS0_short,tgamair_short,tgamself_short,&
                    tnair_short,tn0_short, nol, itrue, lum,ds,n,t,da,i,nc,dz)
implicit none                    
integer, intent(in) :: itrue, nc
real*8, intent(in) :: p0, rexo, H, kb,te,tC_short(itrue),tab_short(itrue),tS0_short(itrue),tgamair_short(itrue)
real*8, intent(in) :: tgamself_short(itrue),tnair_short(itrue),tn0_short(itrue),nol,n(nc+1),dz
real*8, intent(inout) :: t,lum,zp, da,ds
real*8  ::  p, nz
integer, intent(inout) :: i
real*8 :: k
real*8 :: planck

   

    p = p0*exp(-(zp-rexo)/H)

    nz = (p*101315/(kb*Te))*1E-6 !molec/cm**3

    call k2_short(p,nz,te,tC_short,tab_short,tS0_short,tgamair_short,tgamself_short,tnair_short,tn0_short,nol,k,itrue)

    Lum = Lum*exp(-k*ds*100)+planck(nol,Te)*(1-(exp(-k*ds*100)))

    t = asin(n(i)/n(i+1)*sin(t-da)) !changement d'angle
    da = 0
    ds = 0

    i = i + 1 ! passage dans la couche superieur

    zp = (i-1) * dz + Rexo !fixation du nouveau point



end subroutine montee




subroutine descente(p0, zp, rexo,H,kb,Te, tC_short,tab_short,tS0_short,tgamair_short,tgamself_short,&
                    tnair_short,tn0_short, nol, itrue, lum,ds,n,t,da,i,nc,dz)
implicit none                    
integer, intent(in) :: itrue, nc
real*8, intent(in) :: p0, rexo, H, kb,te,tC_short(itrue),tab_short(itrue),tS0_short(itrue),tgamair_short(itrue)
real*8, intent(in) :: tgamself_short(itrue),tnair_short(itrue),tn0_short(itrue),nol,n(nc+1),dz
real*8, intent(inout) :: t,lum,zp, da,ds
real*8  ::  p, nz
integer, intent(inout) :: i
real*8 :: k
real*8 :: planck

   

    p = p0*exp(-(zp-rexo)/H)

    nz = (p*101315/(kb*Te))*1E-6 !molec/cm**3

    call k2_short(p,nz,te,tC_short,tab_short,tS0_short,tgamair_short,tgamself_short,tnair_short,tn0_short,nol,k,itrue)

    Lum = Lum*exp(-k*ds*100)+planck(nol,Te)*(1-(exp(-k*ds*100)))

    t = asin(n(i)/n(i-1)*sin(t+da)) !changement d'angle
    da = 0
    ds = 0

    i = i - 1 ! passage dans la couche inferieur

    zp = (i) * dz + Rexo !fixation du nouveau point



end subroutine descente



subroutine transition(p0, zp, rexo,H,kb,Te, tC_short,tab_short,tS0_short,tgamair_short,tgamself_short,&
                    tnair_short,tn0_short, nol, itrue, lum,ds)
implicit none                    
integer, intent(in) :: itrue
real*8, intent(in) :: p0, rexo, H, kb,te,tC_short(itrue),tab_short(itrue),tS0_short(itrue),tgamair_short(itrue)
real*8, intent(in) :: tgamself_short(itrue),tnair_short(itrue),tn0_short(itrue),nol,zp,ds
real*8, intent(inout) :: lum
real*8  :: p, nz
real*8 :: k
real*8 :: planck

   

    p = p0*exp(-(zp-rexo)/H)

    nz = (p*101315/(kb*Te))*1E-6 !molec/cm**3

    call k2_short(p,nz,te,tC_short,tab_short,tS0_short,tgamair_short,tgamself_short,tnair_short,tn0_short,nol,k,itrue)

    Lum = Lum*exp(-k*ds*100)+planck(nol,Te)*(1-(exp(-k*ds*100)))





end subroutine transition

subroutine masse_reduite(mu)
    implicit none
    real*8 , intent(out) :: mu
    integer :: i,molid,gj
    real*8  :: ab,Q,MM,C
    open(200,status = 'old',form = 'formatted',file = 'Parameter.dat')
    mu = 0
    i = 1

    10000 continue


    do while (i <= 2000)



        read(200,1000,err = 10000,end = 20000)molid,ab,Q,gj,MM,C
        mu = mu + ab*MM*C


        i = i+1
    enddo

    20000 continue
    1000 format (i12,e13.5,e14.4,i5,f14.6,f11.5)



    close (200)

end subroutine masse_reduite







subroutine k2_short(p,nz,T,tC,tab,tS0,tgamair,tgamself,tnair,tn0,nol,k,imax)

    implicit none

    integer,intent(in) :: imax
    real*8,dimension(imax),intent(in) :: tS0,tC,tgamair,tgamself,tnair,tn0,tab
    real*8,intent(in) :: nz,nol,T,p
    real*8,intent(out) :: k
    real*8,parameter :: tref = 296.0,kb = 1.3806488E-16,pi = 3.141592653589793238462643383279 !cgs
    real*8, allocatable :: tgama(:),tflor(:),tk(:), tref2,tmp(:)
    real*8,parameter :: Rpi = 1/pi

    allocate(tgama(imax),tflor(imax),tk(imax), tmp(imax))



    tref2 = tref/T
    tmp = (nol-tn0) ** 2



    tgama = p*(tref2**tnair)*(tgamair*(1-tC*tab)+tgamself*tc*tab)



    tflor = Rpi*tgama/((tgama**2)+tmp)



    tk  = tS0*tflor*tC

    k = sum(tk)*nz

    deallocate(tgama,tflor,tk,tmp)

end subroutine k2_short


function planck(no,T) result(L0)

    real*8 :: no,T,L0


    real*8,parameter :: K = 1.38064852E-23,c = 299792458,h = 6.62607015E-34


    L0 = 2E8*h*(c**2)*(no**3)/(exp(100*h*c*no/(k*T))-1) !loi de planck en nombre d'onde





end function planck

subroutine xy(zp,ap)

real*8,intent(in) :: zp,ap
real*8 :: xp,yp

xp = -zp*cos(ap)
yp = zp*sin(ap)
!print*,'xp,yp',xp,yp
write(42,*)xp,yp

endsubroutine xy

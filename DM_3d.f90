! Programa que simula a evolucao temporal de um sistema de NAT atomos interagentes
! Pela tecnica de dinamica molecular. Os atomos interagem via potencial de Lennard-
! Jones e as equacoes de movimento sao integradas via algoritmo de Verlet.

program dinamol
!use ifport
implicit none

! -------------- Declaracao das variaveis do programa -------------------------

integer              :: nx, ny, nz, nat, i, j, k, m, nsamp
integer              :: npassos, npasseq
real*8               :: par, temp, dt, vol, den, tini, virial
real*8               :: Ecin, Epot, Emec, Emmed, Epmed, Ecmed, prmed, dqm
real*8               :: fcx, fcy, fcz, xit, yit, zit, tinst, press
real*8, allocatable  :: rx(:), ry(:), rz(:), vx(:), vy(:), vz(:)
real*8, allocatable  :: fx(:), fy(:), fz(:), fax(:), fay(:), faz(:)
real*8, allocatable  :: rx0(:), ry0(:), rz0(:)
! -----------------------------------------------------------------------------
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Modificacao para aceitar argumentos de entrada
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

integer :: l
character(len=32) :: arg

call getarg(1, arg)
read (arg,*) den

call getarg(2, arg)
read (arg,*) temp


! Arquivoa manuseados pelo programa:

open(7,file='input1.dat',status='unknown')
open(8,file='energia.dat',status='unknown')
open(9,file='dqm.dat',status='unknown')

! Lendo o arquivo de entrada:

read(7,*) nx, ny, nz  
read(7,*) npassos, npasseq
!read(7,*) den
!read(7,*) temp
read(7,*) dt

! Calculando NAT, VOL e DEN:

nat = nx*ny*nz
par = (1/den)**(1.d0/3.d0)
vol = nx*ny*nz*par**3

! Alocando os vetores:

allocate(rx(nat),ry(nat),rz(nat),vx(nat),vy(nat),vz(nat))
allocate(fx(nat),fy(nat),fz(nat),fax(nat),fay(nat),faz(nat))
allocate(rx0(nat),ry0(nat),rz0(nat))

! Gerando as posicoes e as velocidades iniciais: atomos numa 
! Estrutura quadrada e velocidades iniciais aleatorias

m = 1
do i = 0, nx - 1
  do j = 0, ny - 1
    do k = 0, nz - 1
  
  
   rx(m) = i*par
   ry(m) = j*par
   rz(m) = k*par
   vx(m) = 2*rand() - 1
   vy(m) = 2*rand() - 1
   vz(m) = 2*rand() - 1
 
   m = m + 1
   
    enddo
  enddo
enddo

! Definindo RX0, RY0 e RZ0 e as dimensoes laterais da celula de simulacao:

rx0 = rx
ry0 = ry
rz0 = rz
fcx = nx*par
fcy = ny*par
fcz = nz*par

! Rescalonando as velocidades: seus valores serao multiplicados por uma consante
! De modo a resultar numa energia termica KT igual a da temperatura de entrada.

Ecin = 0.
do i = 1, nat
  Ecin = Ecin + 0.5*(vx(i)**2 + vy(i)**2 + vz(i)**2)
enddo

tini = 2.d0*Ecin/dfloat(3*nat)

vx(:) = vx(:)*dsqrt(temp/tini)
vy(:) = vy(:)*dsqrt(temp/tini)
vz(:) = vz(:)*dsqrt(temp/tini)

! Chamdando a subrotina que calcula a energia potencial e as forcas:

fx(:) = 0.0d0
fy(:) = 0.0d0
fz(:) = 0.0d0

call lenjon(nat,nx,ny,nz,par,rx,ry,rz,fx,fy,fz,Epot,den,virial)

fax(:) = fx(:)
fay(:) = fy(:)
faz(:) = fz(:)

! ------------------- ESCREVENDO ALGUMAS INFORMACOES -------------------------

print*, "--------------------- PROGRAMA DinMol2D -------------------------"
print*, "Dinamica  molecular 2D de um sistema de  NAT atomos posicionados "
print*, "Inicialmente numa  rede quadrada e temperatura  TEMP. O Potencial"
print*, "Interacao interatomico que modela a interacao entre os  atomos eh"
print*, "O Lennard-Jones e condicoes periodicas sobre as posicoes atomicas"
print*, "Sao empregadas"
print*, "-----------------------------------------------------------------"
print*
print*
print*, "(i) Parametros de input da simulacao:"
print*
write(*,'(a,i6)')   ' Numero de Atomos do sistema (NAT)           =', nat
write(*,'(a,i8)')   ' Numero de passos de Integracao (NPASSOS)    =', npassos
write(*,'(a,i8)')   ' Numero de passos de equilibracao (NPASSEQ)  =', npasseq
write(*,'(a,f12.4)')' Densidade = Nat/Vol (DEN)                   =', den
write(*,'(a,f12.4)')' Temperatura de Entrada (TEMP)               =', temp
write(*,'(a,f12.4)')' Passo de integracao (DT)                    =', dt
print*

! ---------------- INTEGRACOES DAS EQUACOES DE MOVIMENTO ---------------------

Emmed = 0.0d0
Epmed = 0.0d0
Ecmed = 0.0d0
prmed = 0.0d0
nsamp = 0

do i = 0, npassos

call verletmod(nat,nx,ny,nz,dt,par,rx,ry,rz,vx,vy,vz,fx,fy,fz,fax,fay,faz,Epot,den,virial)

! Atualizando FAX e FAY:

fax(:) = fx(:)
fay(:) = fy(:)
faz(:) = fz(:)

! *************** Escrevendo os resultados a cada 10 passos *******************

if(mod(i,10).eq.0) then

! Calculo do desvio quadratico medio:

dqm = 0.
do j = 1, nat
   xit = rx(j) - rx0(j)
   yit = ry(j) - ry0(j)
   zit = rz(j) - rz0(j)
   xit = xit - dnint(xit/fcx)*fcx
   yit = yit - dnint(yit/fcy)*fcy
   zit = zit - dnint(zit/fcz)*fcz
  dqm = dqm + xit**2 + yit**2 + zit**2
enddo

write(9,*) i, dqm/dfloat(nat)

! Calculo da energia cinetica:

Ecin = 0.
do j = 1, nat
  Ecin = Ecin + 0.5*(vx(j)**2 + vy(j)**2 + vz(j)**2)
enddo

Emec = Ecin + Epot

! Calculando a temperatura instantanea e a pressao:

tinst = (2.d0*Ecin)/(3.d0*nat)
press = den*tinst + virial/dfloat(nat)

write(8,'(i10,4f12.5)') i, Ecin/dfloat(nat), Epot/dfloat(nat), Emec/dfloat(nat), press

! Acumulando as energias (I > NPASSEQ) e Ajustando a temperatura (I < NPASSEQ)

if(i.ge.npasseq) then
nsamp = nsamp + 1
 Epmed = Epmed + Epot
 Ecmed = Ecmed + Ecin
 Emmed = Emmed + Emec
 prmed = prmed + press
 else
 tini = 2.d0*Ecin/dfloat(3*nat)
vx(:) = vx(:)*dsqrt(temp/tini)
vy(:) = vy(:)*dsqrt(temp/tini)
vz(:) = vz(:)*dsqrt(temp/tini)
endif

endif

enddo ! Fim do loop sobre as integracoes
! *****************************************************************************

! Tomando a media das grandezas acumuladas:

Ecmed = Ecmed/dfloat(nsamp)
Epmed = Epmed/dfloat(nsamp)
Emmed = Emmed/dfloat(nsamp)
prmed = prmed/dfloat(nsamp)

! Escrevendo as medias das energias ao longo da simulacao:

print*, "(ii) Medias das grandezas fisicas apos a simulacao:"
print*
write(*,'(a,f12.4)')   ' Temperatura Media da Simulacao <T>    =', 2.d0*Ecmed/dfloat(3*nat)
write(*,'(a,f12.4)')   ' Energia Cinetica Media  <K>/nat       =', Ecmed/dfloat(nat)
write(*,'(a,f12.4)')   ' Energia Potencial Media <P>/nat       =', Epmed/dfloat(nat)
write(*,'(a,f12.4)')   ' Energia Mecanica Media  <E>/nat       =', Emmed/dfloat(nat)
write(*,'(a,f12.4)')   ' Pressao Media  <P>                    =', prmed
print*

! -----------------------------------------------------------------------------
end program dinamol
! -----------------------------------------------------------------------------

subroutine lenjon(nat,nx,ny,nz,par,rx,ry,rz,fx,fy,fz,Epot,den,virial)
implicit none

integer, intent(in) :: nat, nx, ny, nz
real*8, intent(in)  :: rx(nat), ry(nat), rz(nat), par, den
real*8, intent(out) :: fx(nat), fy(nat), fz(nat), Epot, virial
real*8  :: xij, yij, zij, fcx, fcy, fcz
real*8  :: rij14, rij12, rij8, rij6, rij, pref
integer :: i, j

fcx = nx*par
fcy = ny*par
fcz = nz*par
fx(:) = 0.0d0
fy(:) = 0.0d0
fz(:) = 0.0d0
  Epot = 0.0d0
virial = 0.0d0

do i = 1, nat - 1
  do j = i+1, nat
 
! Posicao relativa com CCP:

   xij = rx(i) - rx(j)
   yij = ry(i) - ry(j)
   zij = rz(i) - rz(j)
   xij = xij - dnint(xij/fcx)*fcx
   yij = yij - dnint(yij/fcy)*fcy
   zij = zij - dnint(zij/fcz)*fcz
   rij = dsqrt(xij**2 + yij**2 + zij**2)
   
! Potencias de RIJ:

   rij6 = rij**6
   rij8 = rij**8
  rij12 = rij**12
  rij14 = rij**14

  Epot = Epot + 4.*((1./rij12) - (1./rij6))
  pref = 24.*((2./rij14) - (1./rij8))

! Forcas:

  fx(i) = fx(i) + pref*xij
  fy(i) = fy(i) + pref*yij
  fz(i) = fz(i) + pref*zij
  fx(j) = fx(j) - pref*xij
  fy(j) = fy(j) - pref*yij
  fz(j) = fz(j) - pref*zij

! Termo de virial:

virial = virial + pref*rij*rij

  enddo
enddo

! Completando o termo de virial:

virial = den*virial/3.d0

! ----------------- FIM DA SUBROTINA LENJON ------------------------------
end subroutine lenjon
! ------------------------------------------------------------------------

subroutine verletmod(nat,nx,ny,nz,dt,par,rx,ry,rz,vx,vy,vz,fx,fy,fz,fax,fay,faz,Epot,den,virial)
implicit none

! ----------------- Variaveis da Subrotina -------------------------------
! Externas
integer,intent(in)   :: nat, nx, ny, nz
real*8, intent(in)   :: fax(nat), fay(nat), faz(nat), par, dt, den
real*8, intent(inout):: rx(nat), ry(nat), rz(nat), vx(nat), vy(nat), vz(nat)
real*8, intent(inout):: fx(nat), fy(nat), fz(nat), virial
real*8, intent(inout):: Epot
! Internas e inteiros para loops
integer :: i, j, k
real*8  :: fcx, fcy, fcz
! ------------------------------------------------------------------------

! Dimensoes laterais da celula de simulacao:

fcx = par*nx
fcy = par*ny
fcz = par*nz

! Posicoes em t + Dt:

do j = 1, nat

rx(j) = rx(j) + vx(j)*dt + (fx(j)/2.d0)*(dt*dt)
ry(j) = ry(j) + vy(j)*dt + (fy(j)/2.d0)*(dt*dt)
rz(j) = rz(j) + vz(j)*dt + (fz(j)/2.d0)*(dt*dt)

! Condicoes de contorno periodicas para as posicoes:

rx(j) = rx(j) - dnint(rx(j)/fcx)*fcx
ry(j) = ry(j) - dnint(ry(j)/fcy)*fcy
rz(j) = rz(j) - dnint(rz(j)/fcz)*fcz

enddo

! Chamando a subrotina LENJON para calcular as novas forcas:

call lenjon(nat,nx,ny,nz,par,rx,ry,rz,fx,fy,fz,Epot,den,virial)

! Velocidades em t + dt:

do k = 1, nat
vx(k) = vx(k) + (fx(k) + fax(k))*(dt/2.d0)
vy(k) = vy(k) + (fy(k) + fay(k))*(dt/2.d0)
vz(k) = vz(k) + (fz(k) + faz(k))*(dt/2.d0)
enddo

! ---------------- FIM DA SUBROTINA VERLETMOD ----------------------------
end subroutine verletmod
! ------------------------------------------------------------------------


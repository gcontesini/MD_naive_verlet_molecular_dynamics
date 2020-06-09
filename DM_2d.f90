! Programa que simula a evolucao temporal de um sistema de NAT atomos interagentes
! Pela tecnica de dinamica molecular. Os atomos interagem via potencial de Lennard-
! Jones e as equacoes de movimento sao integradas via algoritmo de Verlet.

program dinamol
!use ifport
implicit none

! -------------- Declaracao das variaveis do programa -------------------------

integer              :: nx, ny, nat, i, j, k, nsamp
integer              :: npassos, npasseq
real*8               :: par, temp, dt, vol, den, tini
real*8               :: Ecin, Epot, Emec, Emmed, Epmed, Ecmed
real*8, allocatable  :: rx(:), ry(:), vx(:), vy(:), fx(:), fy(:), fax(:), fay(:)
! -----------------------------------------------------------------------------

! Arquivoa manuseados pelo programa:

open(7,file='input.dat',status='unknown')
open(8,file='energia.dat',status='unknown')

! Lendo o arquivo de entrada:

read(7,*) nx, ny  
read(7,*) npassos, npasseq
read(7,*) den
read(7,*) temp
read(7,*) dt

! Calculando NAT, VOL e DEN:

nat = nx*ny
par = dsqrt(1/den)
vol = nx*ny*par**2

! Alocando os vetores:

allocate(rx(nat),ry(nat),vx(nat),vy(nat),fx(nat),fy(nat),fax(nat),fay(nat))

! Gerando as posicoes e as velocidades iniciais: atomos numa 
! Estrutura quadrada e velocidades iniciais aleatorias

k = 1
do i = 0, nx - 1
  do j = 0, ny - 1
  
   rx(k) = i*par
   ry(k) = j*par
   vx(k) = 2*rand() - 1
   vy(k) = 2*rand() - 1
 
   k = k + 1

  enddo
enddo

! Rescalonando as velocidades: seus valores serao multiplicados por uma consante
! De modo a resultar numa energia termica KT igual a da temperatura de entrada.

Ecin = 0.
do i = 1, nat
  Ecin = Ecin + 0.5*(vx(i)**2 + vy(i)**2)
enddo

tini = Ecin/dfloat(nat)
vx(:) = vx(:)*dsqrt(temp/tini)
vy(:) = vy(:)*dsqrt(temp/tini)

! Chamdando a subrotina que calcula a energia potencial e as forcas:

fx(:) = 0.
fy(:) = 0.

call lenjon(nat,nx,ny,par,rx,ry,fx,fy,Epot)

fax(:) = fx(:)
fay(:) = fy(:)
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
nsamp = 0

do i = 0, npassos

call verletmod(nat,nx,ny,dt,par,rx,ry,vx,vy,fx,fy,fax,fay,Epot)

! Atualizando FAX e FAY:

fax(:) = fx(:)
fay(:) = fy(:)

! *************** Escrevendo os resultados a cada 10 passos *******************

if(mod(i,10).eq.0) then

! Calculo da energia cinetica:

Ecin = 0.
do j = 1, nat
  Ecin = Ecin + 0.5*(vx(j)**2 + vy(j)**2)
enddo

Emec = Ecin + Epot

write(8,'(i10,3f12.5)') i, Ecin/dfloat(nat), Epot/dfloat(nat), Emec/dfloat(nat) 

! Acumulando as energias

if(i.ge.npasseq) then
nsamp = nsamp + 1
 Epmed = Epmed + Epot
  Ecmed = Ecmed + Ecin
 Emmed = Emmed + Emec
endif

endif

! *****************************************************************************

enddo ! Fim do loop sobre as integracoes

! Escrevendo as medias das energias ao longo da simulacao:

print*, "(ii) Medias das grandezas fisicas apos a simulacao:"
print*
write(*,'(a,f12.4)')   ' Energia Cinetica Media  <K>/nat       =', Ecmed/dfloat(nat*nsamp)
write(*,'(a,f12.4)')   ' Energia Potencial Media <P>/nat       =', Epmed/dfloat(nat*nsamp)
write(*,'(a,f12.4)')   ' Energia Mecanica Media  <E>/nat       =', Emmed/dfloat(nat*nsamp)
print*

! -----------------------------------------------------------------------------
end program dinamol
! -----------------------------------------------------------------------------

subroutine lenjon(nat,nx,ny,par,rx,ry,fx,fy,Epot)
implicit none

integer, intent(in) :: nat, nx, ny
real*8, intent(in)  :: rx(nat), ry(nat), par
real*8, intent(out) :: fx(nat), fy(nat), Epot
real*8  :: xij, yij, fcx, fcy
real*8  :: rij14, rij12, rij8, rij6, rij, pref
integer :: i, j

fcx = nx*par
fcy = ny*par
fx(:) = 0.
fy(:) = 0.
Epot = 0.

do i = 1, nat - 1
  do j = i+1, nat
 
! Posicao relativa com CCP:

   xij = rx(i) - rx(j)
   yij = ry(i) - ry(j)
   xij = xij - dnint(xij/fcx)*fcx
   yij = yij - dnint(yij/fcy)*fcy
   rij = dsqrt(xij**2 + yij**2)
   
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
  fx(j) = fx(j) - pref*xij
  fy(j) = fy(j) - pref*yij

  enddo
enddo

! ----------------- FIM DA SUBROTINA LENJON ------------------------------
end subroutine lenjon
! ------------------------------------------------------------------------

subroutine verletmod(nat,nx,ny,dt,par,rx,ry,vx,vy,fx,fy,fax,fay,Epot)
implicit none

! ----------------- Variaveis da Subrotina -------------------------------
! Externas
integer,intent(in)   :: nat, nx, ny
real*8, intent(in)   :: fax(nat), fay(nat), par, dt
real*8, intent(inout):: rx(nat), ry(nat), vx(nat), vy(nat), fx(nat), fy(nat)
real*8, intent(inout):: Epot
! Internas e inteiros para loops
integer :: i, j, k
real*8  :: fcx, fcy
! ------------------------------------------------------------------------

! Dimensoes laterais da celula de simulacao:

fcx = par*nx
fcy = par*ny

! Posicoes em t + Dt:

do j = 1, nat

rx(j) = rx(j) + vx(j)*dt + (fx(j)/2.d0)*(dt*dt)
ry(j) = ry(j) + vy(j)*dt + (fy(j)/2.d0)*(dt*dt)

! Condicoes de contorno periodicas para as posicoes:

rx(j) = rx(j) - dnint(rx(j)/fcx)*fcx
ry(j) = ry(j) - dnint(ry(j)/fcy)*fcy

enddo

! Chamando a subrotina LENJON para calcular as novas forcas:

call lenjon(nat,nx,ny,par,rx,ry,fx,fy,Epot)

! Velocidades em t + dt:

do k = 1, nat
vx(k) = vx(k) + (fx(k) + fax(k))*(dt/2.d0)
vy(k) = vy(k) + (fy(k) + fay(k))*(dt/2.d0)
enddo

! ---------------- FIM DA SUBROTINA VERLETMOD ----------------------------
end subroutine verletmod
! ------------------------------------------------------------------------


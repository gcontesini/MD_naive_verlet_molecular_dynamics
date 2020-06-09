!Autor Douglas
program Ha !Matriz simétrica
   implicit none
   integer :: k,n=4,i,lin,col , p
   real*8 ::  H
   real*8, dimension(4,4)  :: A
   real*8, dimension(3,3) :: Al
   real*8, dimension(3) :: v, g, x 
   real*8, dimension(3) :: u,w
   
   
   !real*8, allocatable  :: A(:,:) ,Al(:,:) , x(:) , v(:), g(:), u(:),w(:)
   !allocate( A(n,n) ,Al(n-1,n-1) , x(n),v(n), g(n), u(n-1), w(n-1) )
 
   A(1,:)=(/ 7. , 2. , 3., -1. /)
   A(2,:)=(/ 2. , 8. , 5., 1. /)
   A(3,:)=(/ 3. , 5. , 12. , 9. /)
   A(4,:)=(/ -1. , 1. , 9., 7. /)
   print*, "Matriz A: "
   do i=1 , n
   	print*, A(i,:)
   end do
   !------------------------------------------Comando de repetição---------------------------------------------
   do i=1, n-1

  	!-------------------------------------Criando A' ------------------------------------------------------
  	 
   	do lin=1, n-i
   		do col=1, n-i
   			Al(lin,col) = A(lin+1,col+1)
   		end do
   	end do
   	
   	!Imprimindo A'
   	print*, "Matriz A': "
   	do lin=1 , n-i
 	  	print*, Al(lin,:)
        end do
      
   	
   	
   	
       ! Encontrando X
       x(1: n-i) = A(i+1:n,i)
       
       do lin=1 , n-i
         print*, "x",lin,":",x(lin)
       end do
   
       k=dsign(dsqrt(dot_product(x,x)),x(1))
       !print*, 'k=',k
       !if ( x(1) >= 0. ) then
       ! 	k= ( dsqrt(dot_product(x,x)) )
       !end if
       !if (x(1) < 0. ) then 
       !k = - ( dsqrt(dot_product(x,x)) )   
       !end if

   	u(1) = k+x(1)
   	u(2: n-i ) = x(2:n-i)
   	do lin=1 , n-i
           print*, "u",lin,":",u(lin)
        end do
   	!H= dot_product(u,u)/2.
   	H= (dot_product(u,u))/2.
   	v= matmul(Al,u)/(H)
   	g= dot_product(u,v) / (H*2.)
   	w= v - g*u
   	!Laço
   	
   	do col= 1 , n-i
  		do lin= 1, n-i
  			
    			Al(lin,col) =  Al(lin, col)-w(lin)*u(col)-u(lin)*w(col)		
    			
   		end do
   	end do
   	!------------Imprimindo -------------
   	do col= 1 , n-i
  		do lin= 1, n-i
  			
    			print*, " Imprimir : " ,Al(lin,col)
    			
   		end do
   	end do
   	
   	A(i,i+1) = -k
   	A(i+1,i) = -k
	A(i,i+2:n)= 0.
   	A(i+2:n,i)= 0.
   	
   	
   end do
   
   
   !Imprimindo resultados---------------------------------------------------------------------
   
   print*, "Matriz A': "
   	do lin=1 , n-i
 	  	print*, Al(lin,:)
   end do
   print*, "Matriz A: "
   do lin=1 , n
   	print*, A(lin,:)
   end do
   
end program

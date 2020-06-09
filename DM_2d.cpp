#include "iostream"
#include "fstream"
#include "stdlib.h"
#include "math.h"
#include "time.h"

using namespace std;

main()
{
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Declaracao de variavel

  int i_int, j_int, k_int, numero_atomos_int, size_x_int, size_y_int,          \
    numero_passos_integracao_int, numero_passos_equilibrio_int                 ;

  double temperatura_double,temperatura_inicial_double, parametro_rede_double, \
    volume_double, variacao_tempo_double, densidade_rede_double,               \
    energia_cinetica_double, energia_potencial_double, energia_mecanica_double,\
    energia_cinetica_media_double, energia_potencial_media_double,             \
    energia_mecanica_media_double ;
  
  double *posicao_x_vector, *posicao_y_vector, *velocidade_x_vector,           \
    *velocidade_y_vector, *forca_x_vector, *forca_y_vector,                    \
    *forca_aux_x_vector, *forca_aux_y_vector ;

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Aqruivos

  ofstream energia_data_file                                                   ;
  energia_data_file.open ("energia_data.txt")                                  ;
  ifstream input("input.txt")                                                  ;

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Leitura de Aqruivo

  if (input.is_open())
  {
    while ( getline(input,line) )
    {
      cin >> size_x_int                                                        ;
      cin >> size_y_int                                                        ;
      cin >> numero_passos_integracao_int                                      ;
      cin >> numero_passos_equilibrio_int                                      ;
      cin >> densidade_rede_double                                             ;
      cin >> temperatura_inicial_double                                        ;
      cin >> variacao_tempo_double                                             ;
    }
    input.close();
  }

  else cout << "Unable to open file"                                           ; 

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Escopo de funcao

  double potencial_leonard_jones( 
    int size_x_int              ,
    int size_y_int              ,
    int numero_atomos_intm      ,
    double parametro_rede_double,

  );

  double velert_modificado(

  );

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Calculo inicial

  numero_atomos_int = size_x_int * size_y_int                                  ;
  parametro_rede_double = sqrt( 1/densidade_rede_double )                      ;
  volume_double = size_x_int * size_y_int * pow( parametro_rede_double, 2)     ;               

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Alocacao dinamica
  
  posicao_x_vector = (double*)malloc(numero_atomos_int*sizeof(double))         ;
  posicao_y_vector = (double*)malloc(numero_atomos_int*sizeof(double))         ;
  
  velocidade_x_vector = (double*)malloc(numero_atomos_int*sizeof(double))      ;
  velocidade_y_vector = (double*)malloc(numero_atomos_int*sizeof(double))      ;

  forca_x_vector = (double*)malloc(numero_atomos_int*sizeof(double))           ;
  forca_y_vector = (double*)malloc(numero_atomos_int*sizeof(double))           ;

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Geracao Posicao Velocidade
  
  k_int = 0
  for(i_int=0; i_int < size_x_int; i_int++)
  {
    for( j_int=0; j_int < size_y_int; j_int++ )
    {
      posicao_x_vector[ k_int ] = i_int * parametro_rede_double                ;
      posicao_y_vector[ k_int ] = j_int * parametro_rede_double                ;
      velocidade_x_vector[ k_int ] = 2 * srand(time(NULL)) - 1                 ; 
      velocidade_y_vector[ k_int ] = 2 * srand(time(NULL)) - 1                 ;
      k_int++                                                                  ;
    }
  }

  energia_cinetica_double = 0.0
  for( i_int=0; i_int<numero_atomos_int; i_int++)
  {
    energia_cinetica_double = energia_cinetica_double +                        \
    0.5 * ( pow( velocidade_x_vector[ i_int ], 2) +                            \
    pow( velocidade_y_vector[ i_int ], 2) )
  }

  temperatura_inicial_double = \
  energia_cinetica_double /(double)numero_atomos_int                           ;

  for( i_int=0; i_int<numero_atomos_int; i_int++ )
  {
    velocidade_x_vector[ i_int ] = velocidade_x_vector[ i_int ] *              \
    sqrt(temperatura_double / temperatura_inicial_double)                      ;
    velocidade_y_vector[ i_int ] = velocidade_y_vector[ i_int ] *              \
    sqrt(temperatura_double / temperatura_inicial_double)                      ;
    
    forca_x_vector[i_int] = 0.0                                                ;
    forca_y_vector[i_int] = 0.0                                                ;

    forca_aux_x_vector[i_int] = 0.0                                            ;
    forca_aux_y_vector[i_int] = 0.0                                            ;    
  }

  cout <<
  << "\n--------------------- PROGRAMA DinMol2D -------------------------"
  << "\nDinamica  molecular 2D de um sistema de  NAT atomos posicionados "
  << "\nInicialmente numa  rede quadrada e temperatura  TEMP. O Potencial"
  << "\nInteracao interatomico que modela a interacao entre os  atomos eh"
  << "\nO Lennard-Jones e condicoes periodicas sobre as posicoes atomicas"
  << "\nSao empregadas"
  << "\n-----------------------------------------------------------------"
  << "\n"
  << "\n(i) Parametros de input da simulacao:"
  << '\nNumero de Atomos do sistema      ='<< numero_atomos_int
  << '\nNumero de passos de Integracao   ='<< numero_passos_integracao_int
  << '\nNumero de passos de equilibracao ='<< numero_passos_equilibrio_int
  << '\nDensidade                        ='<< densidade_rede_double
  << '\nTemperatura de Entrada           ='<< temperatura_double
  << '\nPasso de integracao              ='<< variacao_tempo_double
  << endl                                                                      ;

  energia_mecanica_media_double = 0.0                                          ;
  energia_potencial_media_double = 0.0                                         ;
  energia_cinetica_media_double = 0.0                                          ;

  // nsamp = 0;

  for( i_int=0; i_int<numero_passos_integracao_int; i_int++)
  {
    velert_modificado(
      numero_atomos_int         ,
      size_x_int                ,
      size_y_int                ,
      variacao_tempo_double     ,
      parametro_rede_double     ,
      posicao_x_vector          ,
      posicao_y_vector          ,
      velocidade_x_vector       ,
      velocidade_y_vector       ,
      forca_x_vector            ,
      forca_y_vector            ,
      forca_aux_x_vector        ,
      forca_aux_y_vector        ,
      energia_potencial_double
    )

    for(k_int=0; k_int<numero_atomos_int; k_int++)
    {
      forca_aux_x_vector[ k_int ] = forca_x_vector[ k_int ]                    ;
      forca_aux_y_vector[ k_int ] = forca_y_vector[ k_int ]                    ;
    }

    if( (double)i_int % 10 == 0 )
    {
      energia_cinetica_double = 0.0
      for( j_int=0; j_int<numero_atomos_int; j_int++ )
      {
        energia_cinetica_double = energia_cinetica_double + 0.5*               \
        (pow( velocidade_x_vector[j_int], 2 ) +                                \
        pow( velocidade_y_vector[j_int], 2 ))                                  ;        
      }

      energia_mecanica_double =                                                \
      energia_cinetica_double + energia_potencial_double                       ;


      energia_data_file <<
      << i_int 
      << energia_cinetica_double/(double)numero_atomos_int
      << energia_potencial_double/(double)numero_atomos_int
      << energia_mecanica_double/(double)numero_atomos_int
      << endl                                                                  ;

      if( i_int > numero_passos_equilibrio_int )
      {
        // nusamp = nusamp +1
        energia_potencial_media_double =                                       \
        energia_potencial_media_double + energia_potencial_double              ;
        
        energia_cinetica_media_double =                                        \
        energia_cinetica_media_double + energia_cinetica_double                ;

        energia_mecanica_media_double =                                        \
        energia_mecanica_media_double + energia_mecanica_double                ;
      }
    }
  }

  cout <<
  <<  "(ii) Medias das grandezas fisicas apos a simulacao:"
  <<  '\nEnergia Cinetica Media  =', energia_mecanica_media_double/            \
  (double)(numero_atomos_int * nsamp)
  <<  '\nEnergia Potencial Media =', energia_potencial_media_double/           \
  (double)(numero_atomos_int* nsamp)
  <<  '\nEnergia Mecanica Media  =', energia_mecanica_media_double/            \
  (double)(numero_atomos_int * nsamp)
  << endl                                                                      ;
  energia_data_file.close();
}

// nsamp = 
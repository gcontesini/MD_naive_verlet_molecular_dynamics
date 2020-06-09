// ConsoleApplication2.cpp : Defines the entry point for the console application.
//
// #include "stdio.h"
#include "math.h"
#include "iostream"
// #include "time.h"
#include "stdlib.h"

using namespace std;

FILE *arq;

double lenjon(int nat, double *rx, double *ry, double *fx, double *fy)
{
  //Declarando variaveis
  int i, j;
  double epot, rij, rij12, rij6, rij14, rij8, pref;
  epot = 0.;
  //fij força em i devido a j
  for (i = 1; i <= nat ; i++)
  {
    for (j = i + 1; i <= nat ; i++)
    {
      rij = sqrt( pow( ( rx[i]- rx[j] ) ,2) + pow( ( ry[i] - ry[j] ) , 2 ) );
      rij6 = pow(rij,6);
      rij8 = pow(rij, 8);
      rij12 = pow(rij, 12);
      rij14 = pow(rij, 14);
      //Fórmula de Lennar Jones
      epot = epot + 4 * ((1.0 / rij12) - (1.0 / rij6));
      pref = - 24. * ((2.0 / rij14) - (1.0 / rij8)); // prefator
      //forças
      fx[i] = fx[i] + pref* (rx[j] - rx[i]);
      fy[i] = fy[i] + pref* (ry[j] - ry[i]);
      fx[j] = fx[j] - pref* (rx[j] - rx[i]);
      fy[j] = fy[j] - pref* (ry[j] - ry[i]);
      
    }
  }
  return epot;
}

int main(/*int argc, _TCHAR* argv[]*/)
{
  //Declarando as variáveis 
  int nx=4, ny=4, nat, i, j, k;
  double *rx, *ry, *fx, *fy, par=2.0, temp=300., dt=0.005, vol, den, ecin, tini, epot,emec;
  double *vx,*vy;
  double med, var;
  //char url[] = "input.txt";
  //abrindo arquivos

  /*
  arq = fopen(url, "w"); // criar outro projeto com gerenciamento de arquivo do pdf
  
  
  //Leitura do arquivo de entrada
  
  if ( fopen("input.txt", "r") == NULL)
  {
    cout << "Erro: Arquivo vazio!";
    exit(1);
  }
  fscanf(arq, "%d %d %df %df %df", &nx, &ny,&par, &temp, &dt);
  */
  
  //Calculando NAT, VOL e DEN :

  nat = nx*ny;
  vol = nx*ny*pow(par,2);
  den = ((double) nat) / vol;


  // leitura do número de valores 
  //cout << "Digite o numero de elementos: " << endl;
  //scanf("%d", &n);

  //Alocando dinâmicamente os vetores
  rx = (double*)malloc(nat*sizeof(double));
  ry = (double*)malloc(nat*sizeof(double));
  vx = (double*)malloc(nat*sizeof(double));
  vy = (double*)malloc(nat*sizeof(double));
  fx = (double*)malloc(nat*sizeof(double));
  fy = (double*)malloc(nat*sizeof(double));

  k = 1;
  srand(rand());
  for (i = 1; i <= nx ; i++)
  {
    for (j = 1; j <= ny ; j++)
    {
      rx[k] = i*par;
      ry[k] = j*par;
      vx[k] = 2 * ( rand()  % 100)/100. ;
      vy[k] = 2 * ( rand()  % 100)/100. ;
      cout << " " << rx[k] << " " << ry[k] << " " << vx[k] << " " << " " << vy[k] << endl;
      k += 1;
    }
  }
  //Rescalonando as velocidades
  ecin = 0.;
  for (i =1; i <= nat; i++)
  {

    ecin = ecin + 0.5*(pow(vx[i],2) + pow(vy[i],2));
    cout << ecin << endl;
  }
  tini = ecin /((double) nat);
  for (i = 1; i <= nat; i++)
  {
    vx[i] = vx[i] * sqrt((double) (temp*tini));
    vy[i] = vy[i] * sqrt((double) (temp*tini));
    cout << "Vx: " << vx[i] << " Vy: " << vy[i] <<endl;
    fx[i] = 0.;
    fy[i] = 0.;
  }
  //Chamandi a subrotina que calcula a energia potencial e as forças
  cout << "Energias do sistema" << endl;
  epot = lenjon(nat, rx, ry, fx, fy);
  cout << "\nEnergia Cinetica: "<< ecin << "\nEnergia Potencial: " << epot << "\nEnergia Total: " <<  (ecin + epot) <<endl;
  //free(rx);
  //free(ry);
  //free(fx);
  //free(fy);
  //free(vx);
  //free(vy);
  // system("pause");
  return 0;
}


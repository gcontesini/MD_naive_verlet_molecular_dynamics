#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys

def main():

  densidade_inicial_float   = float(input("Densidade inicial:")) 
  densidade_final_float     = float(input("Densidade final:"))  
  temperatura_inicial_float = float(input("Temperatura inicial:"))
  temperatura_final_float   = float(input("Temperatura final:"))

  delta_density_float     = ( densidade_final_float - densidade_inicial_float )
  delta_temperatura_float = ( temperatura_final_float - \
  	temperatura_inicial_float )

  dt_float = 0.1

  for j_int in range( int( delta_density_float*10 )+1 ):

    for i_int in range( int( delta_temperatura_float*10 )+1 ):

      update_progress( float( float(j_int)*float(i_int) )/\
      	( 10.*delta_density_float*10.*delta_temperatura_float ) ) 

      os.system(
        "./DM_3d "+str(densidade_inicial_float+(j_int*dt_float))+" "+\
        str(temperatura_inicial_float+(i_int*dt_float))+\
        ' > '+str("saida_"+str(densidade_inicial_float+(j_int*dt_float))+"_"+\
        	str(temperatura_inicial_float+(i_int*dt_float))+'.txt')
      )

  print '\n'

def update_progress( progress ):
  
  barLength = 10 # Modify this to change the length of the progress bar
  status = ""
  if isinstance(progress, int):
    progress = float(progress)
  if not isinstance(progress, float):
    progress = 0
    status = "error: progress var must be float\r\n"
  if progress < 0:
    progress = 0
    status = "Halt...\r\n"
  if progress >= 1:
    progress = 1
    status = "Done...\r\n"
  block = int(round(barLength*progress))
  text = "\rPercent: [{0}] {1}% {2}".format( "="*block + " "*\
  	(barLength-block), int(progress*100), status)
  sys.stdout.write(text)
  sys.stdout.flush()

if __name__ == "__main__":

  main()
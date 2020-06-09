#!/bin/sh

grep 'Temperatura Media da Simulacao <T>' saida_0.75_*.txt | awk -F "=" '{ print $2 }' | sed 's/      //g' > status_07a.txt
grep 'Energia Total Media <T>/nat' 		    saida_0.75_*.txt | awk -F "=" '{ print $2 }' | sed 's/      //g' > status_07b.txt

paste status_07a.txt status_07b.txt | awk '{ print $1,$2 }' > status_07.txt
rm status_07a.txt
rm status_07b.txt

grep 'Temperatura Media da Simulacao <T>' saida_0.85_*.txt | awk -F "=" '{ print $2 }' | sed 's/      //g' > status_08a.txt
grep 'Energia Total Media <T>/nat' 		    saida_0.85_*.txt | awk -F "=" '{ print $2 }' | sed 's/      //g' > status_08b.txt

paste status_08a.txt status_08b.txt | awk '{ print $1,$2 }' > status_08.txt
rm status_08a.txt
rm status_08b.txt

grep 'Temperatura Media da Simulacao <T>' saida_0.95_*.txt | awk -F "=" '{ print $2 }' | sed 's/      //g' > status_09a.txt
grep 'Energia Total Media <T>/nat' 		    saida_0.95_*.txt | awk -F "=" '{ print $2 }' | sed 's/      //g' > status_09b.txt

paste status_09a.txt status_09b.txt | awk '{ print $1,$2 }' > status_09.txt
rm status_09a.txt
rm status_09b.txt

rm -rf saida_*.txt

gnuplot gnuscript.gp

reset                                           
set terminal png
set output 'status.png'         
set xlabel "Temperature (K)"
set ylabel "Total Energy (meV)"
set title "Temperatura x Energia Total"
set multiplot

plot "status_07.txt" u 1:2 w lp title "densidade 0.7"  lc rgb "blue" ;
rep "status_08.txt" u 1:2 w lp title "densidade 0.8"  lc rgb "red" ;
rep "status_09.txt" u 1:2 w lp title "densidade 0.9"  lc rgb "black"
reset
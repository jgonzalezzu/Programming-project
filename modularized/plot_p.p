set xlabel 't [year]'
set ylabel 'concentration of CO_2 (t)'
set xtics 1000
set ytics 1
plot 'Datos.txt' using 1:2 w l t'P(t)' ,\


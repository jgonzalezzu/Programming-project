set xlabel 't [year]'
set ylabel 'emission of CO_2 f(t)'
set xtics 1000
set ytics 2
plot 'Datos.txt' using 1:5 w l t'f(t)'
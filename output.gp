set   autoscale                        # scale axes automatically
set title "How The Temperature Distribution of Nuclear Fuel Rod Varies Over Time"
set xlabel "Radius(cm)"
set ylabel "Temperature(K)"
set xr [0:100]
plot "output.dat" using 1:2 title 'Time = 1 year' with lines,  "output.dat" using 1:3 title 'Time = 10 years' with lines,  "output.dat" using 1:4 title 'Time = 50 years' with lines,  "output.dat" using 1:5 title 'Time = 100 years' with lines ,
pause -1

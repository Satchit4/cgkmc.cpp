set datafile separator comma
set terminal pngcairo size 1280,900 enhanced font "Arial,24"
set border linewidth 2
set grid xtics ytics linewidth 1 linecolor rgb "#bdbdbd"
set key top right box opaque width 1 height 1
set xlabel "time (units of kappa^-1)"
set ylabel "surface energy (mJ/m^2)"
set ytics 80,20,210
set yrange [80:210]

set output "results/figure2_petn_surface_energy_full.png"
set xrange [0:75]
set xtics 0,15,75
plot "results/petn_surface_energy_ovito.csv" every ::1 using 4:9 with lines linewidth 3 linecolor rgb "#ff7f0e" title "KMC simulation", \
     88 with lines dashtype 2 linewidth 3 linecolor rgb "#111111" title "AE model"

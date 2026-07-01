set datafile separator comma
set terminal pngcairo size 1280,900 enhanced font "Arial,24"
set border linewidth 2
set grid xtics ytics linewidth 1 linecolor rgb "#bdbdbd"
set key top right box opaque width 1 height 1
set xlabel "time (units of kappa^-1)"
set ylabel "surface energy (mJ/m^2)"
set autoscale xy

set output "results/figure2_trp_l_tryptophan_surface_energy.png"
plot "results/trp_l_tryptophan_surface_energy_ovito.csv" every ::1 using 4:9 with lines linewidth 3 linecolor rgb "#1f77b4" title "TRP surrogate KMC"

set output "results/figure2_trp_l_tryptophan_surface_energy_time_s.png"
plot "results/trp_l_tryptophan_surface_energy_ovito.csv" every ::1 using 3:9 with lines linewidth 3 linecolor rgb "#1f77b4" title "TRP surrogate KMC"

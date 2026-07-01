set datafile separator whitespace
set terminal pngcairo size 1500,1250 enhanced font "Arial,18"
set output "results/figure3_petn_final_morphology.png"
set multiplot layout 2,2 margins 0.03,0.97,0.05,0.97 spacing 0.04,0.16

unset key
unset colorbox
unset border
unset xtics
unset ytics
unset ztics
unset xlabel
unset ylabel
unset zlabel
set xyplane 0
set view equal xyz
set xrange [-165:165]
set yrange [-165:165]
set zrange [-235:235]

set view 62,38
splot "results/petn_final_points.tsv" every ::1 using 2:3:4 with points pointtype 7 pointsize 0.13 linecolor rgb "#ef6b6b"

set view 66,126
splot "results/petn_final_points.tsv" every ::1 using 2:3:4 with points pointtype 7 pointsize 0.13 linecolor rgb "#ef6b6b"

set view 62,38
splot "results/petn_final_surface_points.tsv" every ::1 using 2:3:4 with points pointtype 7 pointsize 0.42 linecolor rgb "#9b45c6"

set view 66,126
splot "results/petn_final_surface_points.tsv" every ::1 using 2:3:4 with points pointtype 7 pointsize 0.42 linecolor rgb "#9b45c6"

unset multiplot

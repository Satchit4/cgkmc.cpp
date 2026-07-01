set datafile separator whitespace
set terminal pngcairo size 1500,1250 enhanced font "Arial,18"
set output "results/figure3_trp_l_tryptophan_final_morphology.png"
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
set autoscale x
set autoscale y
set autoscale z

set view 62,38
splot "results/trp_l_tryptophan_final_points.tsv" every ::1 using 2:3:4 with points pointtype 7 pointsize 0.13 linecolor rgb "#4e79a7"

set view 66,126
splot "results/trp_l_tryptophan_final_points.tsv" every ::1 using 2:3:4 with points pointtype 7 pointsize 0.13 linecolor rgb "#4e79a7"

set view 62,38
splot "results/trp_l_tryptophan_final_surface_points.tsv" every ::1 using 2:3:4 with points pointtype 7 pointsize 0.42 linecolor rgb "#59a14f"

set view 66,126
splot "results/trp_l_tryptophan_final_surface_points.tsv" every ::1 using 2:3:4 with points pointtype 7 pointsize 0.42 linecolor rgb "#59a14f"

unset multiplot

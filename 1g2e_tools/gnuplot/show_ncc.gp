# Parameters
id = 2

# set terminal png transparent nocrop enhanced font arial 8 size 420,320 
# set output 'fillbetween.3.png'
set term x11 enhanced font "Times-New-Roman,16" 

#set style fill solid 0.4 noborder
set style data lines
#unset title 

file_dic = "dic_out_".id.".txt"

max(x,y) = (x>y)?x:y


# Note that NCC_subpix should not be compared with NCC since the pattern is not the same.
#plot file_dic u 0:13:(max($11,$12)) w filledcurves above title 'SubPix. Enhanced' lc rgb "green" , \
#     file_dic u 0:13:(max($11,$12)) w filledcurves below title 'Below' lc rgb "red" 

#plot file_dic u 0:(max($12,$11)):11 w filledcurves above title 'Rescue Enhanced' lc rgb "green" , \
#     file_dic u 0:(max($12,$11)):11 w filledcurves below title 'Below' lc rgb "red", \
#     file_dic u 0:11 w l lc rgb "gray"

set key bottom left
set xlabel "Tracked point identifier"
set ylabel "NCC"
set title "File number ".id
plot [][0:1] file_dic u 13 w l title "Subpix." lw 1.2 lc rgb "red" 
replot file_dic u (max($12,$11)) w l title "Rescued" lw 1.2 lc rgb "blue" 
replot file_dic u 11 w l title "Followed" lw 1.2 lc rgb "gray" 

pause -1

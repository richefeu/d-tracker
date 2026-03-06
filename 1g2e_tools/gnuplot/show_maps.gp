# usage:

# > load "path/to/show_raw.gp"
# > id = 2
# > fact = 100
# > plot @NCC_rescue, @disp

# Default Parameters
id = 2        # File number
fact = 200.0  # Factor for the size of vectors

#_________________________________________________
set term wxt size 1000,750
set size ratio -1
unset key

## MACROS
set macros

ZOOMCOL0 = 'set palette defined (0.0 "black", 0.5 "red", 1.0 "green"); set cbrange [0:1]'
ZOOMCOL1 = 'set palette defined (0.5 "black", 0.75 "red", 1.0 "green"); set cbrange [0.5:1]'
ZOOMCOL2 = 'set palette defined (0.75 "black", 0.875 "red", 1.0 "green"); set cbrange [0.75:1]'
RODS = '"dic_out_".id.".txt" u ($1+$5):(-$2-$6):4 w circles lc rgb "blue" fs transparent solid 0.15'
NCC_follow = '"dic_out_".id.".txt" u ($1+$5):(-$2-$6):4:11 w circles lc palette fs transparent solid 0.15'
NCC_rescue = '"dic_out_".id.".txt" u ($1+$5):(-$2-$6):4:12 w circles lc palette fs transparent solid 0.15'
NCC = '"dic_out_".id.".txt" u ($1+$5):(-$2-$6):4:13 w circles lc palette fs transparent solid 0.15'
disp = '"dic_out_".id.".txt" u ($1+$5):(-$2-$6):($8*fact):(-$9*fact) w vectors lw 1.2 lc rgb "black"'

@ZOOMCOL2



# Title 
set title 'Temperature'

set dgrid3d 30,30
set hidden3d
splot "Temperature.dat" u 1:2:3 with lines
pause -1 "Hit any key to continue"

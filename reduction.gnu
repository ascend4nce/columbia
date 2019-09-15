set tics font ",30"
set terminal png size 1920, 1080
set pm3d map
set pm3d scansforward
set palette defined (0 "blue", 1 "cyan", 2 "green", 3 "yellow", 4 "red")
set nokey

set bmargin at screen 0.07
set tmargin at screen 0.95
set lmargin at screen 0.07
set rmargin at screen 0.85

filenames = "rdvalue rdfluxx rdfluxy rdtruevalue rdtruefluxx rdtruefluxy rdvalueerror rdfluxxerror rdfluxyerror"
do for [file in filenames] {
	outfile = "reduction/".file.".png"
	set output outfile
	splot "reduction/".file.".dat" u 1:2:3 with pm3d
}

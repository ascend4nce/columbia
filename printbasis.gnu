set terminal png size 1920, 1080
set pm3d map
set pm3d scansforward
set palette defined (0 "blue", 1 "cyan", 2 "green", 3 "yellow", 4 "red")

filenames = system('ls baseprints/*fluxx*.dat baseprints/*fluxy*.dat baseprints/*value*.dat')

do for [file in filenames] {
	outfile = file
	commands = "echo ".outfile." | sed s=.dat=.png=g"
	outfile = system(commands) 
	set output outfile
	splot file."" u 1:2:3 with pm3d
}

filenames = system('ls baseprints/*nodes.dat')

set nokey 

do for [file in filenames] {
	outfile = file
	commands = "echo ".outfile." | sed s=.dat=.png=g"
	outfile = system(commands) 
	set output outfile
	splot file."" u 1:2:3 with points palette pt 7 ps 2.5
}

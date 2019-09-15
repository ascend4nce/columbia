set terminal png size 1920, 1080
set pm3d map
set pm3d scansforward
#set pm3d interpolate 0,0
set palette defined (0 "blue", 1 "cyan", 2 "green", 3 "yellow", 4 "red")

filenames = system('ls stest/*.dat')

do for [file in filenames] {
	outfile = file
	commands = "echo ".outfile." | sed s=.dat=.png=g"
	outfile = system(commands) 
#	outfile[strlen(outfile) - 3] = 'p'
#	outfile[strlen(outfile) - 2] = 'n'
#	outfile[strlen(outfile) - 1] = 'g'
	set output outfile
	splot file."" u 1:2:3 with pm3d
}

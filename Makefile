all:
	gcc -c src/dcdplugin.c -lm -Iinclude
	gfortran -c src/trajectory.f90 -Iinclude
	g++ dcdplugin.o OLD/main.cpp -o test -Iinclude
	gfortran dcdplugin.o trajectory.o main.f90 -Jinclude -Iinclude

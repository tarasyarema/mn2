run: compile
	./main	

compile:
	gcc -ansi -Wall -pedantic -g YaremaTaras.c -o main -lm

plot:
	gnuplot -e "set terminal png size 635,635; set output 'got.png'; set grid; set size square; set size ratio -1;  plot 'results.txt' w l"

docs:
	latexmk -quiet -pdf YaremaTaras.tex
	latexmk -quiet -c YaremaTaras.tex


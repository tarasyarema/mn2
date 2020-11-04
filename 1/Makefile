run: compile
	./main	

test: setup_python compile
	./main
	python test.py
	
setup_python:
	python -m pip install numpy==1.19.3

compile:
	gcc -ansi -Wall -pedantic YaremaTaras.c -o main -lm

docs:
	latexmk -quiet -pdf YaremaTaras.tex
	latexmk -quiet -c YaremaTaras.tex


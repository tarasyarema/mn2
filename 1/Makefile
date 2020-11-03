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
	pandoc YaremaTaras.md -s --pdf-engine=xelatex -o YaremaTaras.pdf


PREFIX ?= $(HOME)/.local/

pcor: pcor.cpp
	g++ -O3 -std=c++0x -fopenmp -o $@ $^ -larmadillo

install: pcor
	cp pcor $(PREFIX)/bin/

TexOptions = --pdf --interaction=nonstopmode --halt-on-error --output-directory=build

Target = deflection_article

all: build/$(Target).pdf

build/$(Target).pdf: FORCE | build
	latexmk $(TexOptions) $(Target).tex 

build:
	mkdir -p build 

clean:
	rm -rf build 

FORCE:

.PHONY: all clean FORCE

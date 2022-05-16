TexOptions = --pdf --interaction=nonstopmode --halt-on-error --output-directory=build

Target = elsarticle-template-num

all: build/$(Target).pdf

build/$(Target).pdf: FORCE | build
	latexmk $(TexOptions) $(Target).tex 

build:
	mkdir -p build 

clean:
	rm -rf build 

FORCE:

.PHONY: all clean FORCE

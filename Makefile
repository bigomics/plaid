build: doc
	R -e "devtools::build()"

doc: vignettes/*.Rmd R/*.R
	R -e "devtools::document()"
	R -e "devtools::build_vignettes()"

install: 
	R CMD INSTALL .

check: clean
	R -e 'devtools::check()'

clean:
	rm -f `find . -name '.#*' -o -name '#*' -o -name '*~' -printf '"%p" '`

FORCE: ;

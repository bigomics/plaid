build: doc
	R -e "devtools::build()"

doc: vignettes/plaid-vignette.Rmd
	R -e "devtools::document()"
	R -e "devtools::build_vignettes()"

install: 
	R CMD INSTALL .

check: clean
	R -e 'devtools::check()'

clean:
	rm -f `find . -name '.#*' -o -name '#*' -o -name '*~' -printf '"%p" '`

FORCE: ;

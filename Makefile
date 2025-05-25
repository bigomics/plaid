install: 
	R CMD INSTALL .

check: clean
	Rscript -e 'devtools::check()'

clean: 
	rm -f `find . -name '.#*' -o -name '*~' -o -name '#*#' -o -name '*.save' -o -name '*.bak' -o -name '*.SAVE' -o -name '*.BAK'`

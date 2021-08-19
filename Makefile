# (C) 2020-2021 Marek Gagolewski, https://www.gagolewski.com

.PHONY: all check test

all:
	Rscript -e 'Rcpp::compileAttributes()'
	R CMD INSTALL .
	# AVOID ADDING THE -O0 flag!!!
	Rscript -e 'roxygen2::roxygenise(roclets=c("rd", "collate", "namespace", "vignette"), load_code=roxygen2::load_installed)'
	R CMD INSTALL .

check: all
	Rscript -e 'devtools::check()'

test: clean all
	Rscript -e 'options(width=120); testthat::test_dir("tests/")'

test-quick:
	Rscript -e 'options(width=120); testthat::test_dir("tests/")'

clean:
	rm -f src/*.o src/*.so

# Copyleft (C) 2025-2025, Marek Gagolewski <https://www.gagolewski.com/>

# pip3 install tabulate hdbscan pytest rpy2 scanpy
# R: rworldmap [gdal-devel geos-devel proj-devel udunits2-devel gdal proj]

.PHONY: python py-test py-check r r-check r r-build html docs clean

.NOTPARALLEL: python py-test py-check r r-check r r-build clean purge html docs

PKGNAME="quitefastmst"

all: r python

################################################################################

python:
	# python3 setup.py install --user  # deprecated
	# python3 -m build  # requires the 'build' package, creates an isolated envir
	# -DQUITEFASTMST_PROFILER
	CPPFLAGS="${CPPFLAGS} -UNDEBUG -DDEBUG -Wpedantic -Wextra -Wall" python3 -m pip install . --no-build-isolation

py-test: python
	pytest
	cd .devel/sphinx && make doctest && cd ../../

################################################################################

r-autoconf:
	Rscript -e 'Rcpp::compileAttributes()'
	#CXX_DEFS="${CPPFLAGS} -UNDEBUG -DDEBUG -Wall -Wextra -Wpedantic" R CMD INSTALL . --preclean

r: r-autoconf
	CXX_DEFS="${CPPFLAGS} -UNDEBUG -DDEBUG -Wall -Wextra -Wpedantic" R CMD INSTALL . --html
	# Roxygen2 adds the -O0 flag if load_installed is not passed!
	Rscript -e "\
	    source('.devel/roxygen2-patch.R');\
	    roxygenise(\
	        roclets=c('rd', 'collate', 'namespace', 'vignette'), \
	        load_code=roxygen2::load_installed\
	)"
	CXX_DEFS="${CPPFLAGS} -UNDEBUG -DDEBUG -Wall -Wextra -Wpedantic" R CMD INSTALL . --html

r-test: r
	Rscript -e 'source(".devel/tinytest.R")'

r-build: r-autoconf
	cd .. && R CMD build ${PKGNAME}

r-check: r-build
	cd .. && R_DEFAULT_INTERNET_TIMEOUT=240 \
	    _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE \
	    _R_CHECK_FORCE_SUGGESTS_=0 \
	    R CMD check `ls -t ${PKGNAME}*.tar.gz | head -1` --as-cran --no-manual

################################################################################

rd2myst:
	# https://github.com/gagolews/Rd2rst
	cd .devel/sphinx && Rscript -e "Rd2rst::Rd2myst('${PKGNAME}')"

weave-examples:
	cd .devel/sphinx/rapi && Rscript -e "Rd2rst::weave_examples('${PKGNAME}', '.')"

weave:
	cd .devel/sphinx/weave && make && cd ../../../

news:
	cd .devel/sphinx && cp ../../NEWS news.md

html: r python news weave rd2myst weave-examples
	rm -rf .devel/sphinx/_build/
	cd .devel/sphinx && make html
	.devel/sphinx/fix-html.sh .devel/sphinx/_build/html/rapi/
	.devel/sphinx/fix-html.sh .devel/sphinx/_build/html/weave/
	rm -rf .devel/sphinx/_build/html/_sources
	@echo "*** Browse the generated documentation at"\
	    "file://`pwd`/.devel/sphinx/_build/html/index.html"

docs: html
	@echo "*** Making 'docs' is only recommended when publishing the "\
	    "official release, because it updates the package homepage."
	@echo "*** Therefore, we check if the package version is like 1.2.3 "\
	    "and not 1.2.2.9007."
	Rscript --vanilla -e "stopifnot(length(unclass(packageVersion('${PKGNAME}'))[[1]]) < 4)"
	rm -rf docs/
	mkdir docs/
	cp -rf .devel/sphinx/_build/html/* docs/
	cp .devel/CNAME.tpl docs/CNAME
	touch docs/.nojekyll
	touch .nojekyll

################################################################################

clean:
	python3 setup.py clean
	rm -rf quitefastmst/__pycache__/
	rm -rf quitefastmst.egg-info/
	rm -rf ../quitefastmst/
	rm -rf dist/
	rm -rf build/
	rm -rf .devel/pytest/__pycache__/
	rm -f quitefastmst/py_*.cpp
	find src -name '*.o' -exec rm {} \;
	find src -name '*.so' -exec rm {} \;
	rm -rf .devel/sphinx/_build/
	rm -rf .devel/sphinx/rapi/
	rm -rf revdep/

purge: clean
	#rm -f man/*.Rd

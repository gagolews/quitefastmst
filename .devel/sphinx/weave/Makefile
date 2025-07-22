# Copyright (c) 2013-2025, Marek Gagolewski <https://www.gagolewski.com>

FILES_RMD = \
	definitions.Rmd \
	benchmarks.Rmd \
	python.Rmd \
	r.Rmd

RMD_MD_OUTPUTS=$(patsubst %.Rmd,%.md,$(FILES_RMD))

%.md: %.Rmd
	./Rmd2md.sh "$<"

all : rmd

rmd : $(RMD_MD_OUTPUTS)

clean:
	rm -f $(RMD_MD_OUTPUTS)

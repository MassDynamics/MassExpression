installr:
	R -e 'renv::restore()'

check:
	R -e "renv::restore();sessionInfo();devtools::check()"
	
document:
	Rscript -e "devtools::document()"

build:
	Rscript -e "devtools::build()"

install:
	Rscript -e "devtools::install(build_vignettes = TRUE, upgrade = TRUE)"

readme:
	Rscript -e "rmarkdown::render('README.Rmd')"


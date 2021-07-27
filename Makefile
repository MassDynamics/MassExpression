installr:
	Rscript install.R

check:
	Rscript -e "devtools::check()"
	
document:
	Rscript -e "devtools::document()"

build:
	Rscript -e "devtools::build()"

install:
	Rscript -e "devtools::install(build_vignettes = TRUE, upgrade = TRUE)"

readme:
	Rscript -e "rmarkdown::render('README.Rmd')"


FROM rocker/r-ver:4.1.0

#RUN yum -y install libjpeg-turbo-devel
#RUN yum --setopt=skip_missing_names_on_install=False -y install openssl-devel libcurl-devel libxml2-devel
RUN apt-get update
RUN apt-get install libssl-dev libcurl4-openssl-dev libz-dev libxml2-dev -y


ENV RENV_VERSION 0.14.0
RUN R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))"
RUN R -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"

COPY . /usr/local/src/MassExpression
WORKDIR /usr/local/src/MassExpression

RUN R -e 'renv::restore()'


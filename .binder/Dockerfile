FROM rocker/binder:3.6.0
LABEL maintainer='Ben Marwick'
USER root
COPY . ${HOME}
RUN chown -R ${NB_USER} ${HOME}
USER ${NB_USER}



RUN wget https://github.com/benmarwick/signatselect/raw/master/DESCRIPTION && R -e "options(repos = list(CRAN = 'http://mran.revolutionanalytics.com/snapshot/2019-12-03/')); devtools::install_deps(); devtools::install_github('benmarwick/signatselect')"

RUN rm DESCRIPTION.1; exit 0

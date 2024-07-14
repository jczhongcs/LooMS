FROM nvcr.io/nvidia/tensorflow:22.11-tf2-py3


LABEL maintainer="j58zhong@uwaterloo.ca,jczhongcs@gmail.com,jr2wu@uwaterloo.ca"

# FOR PROD Copy PYTHON code into the container
# MUST run mvn clean install package first, otherwise the COPY ./target... stage will fail
COPY ./src/main/python /waterlooms/src/main/python
COPY ./src/main/R /waterlooms/src/main/R

COPY ./DeBruijn/deBruijn.jar /waterlooms/deBruijn.jar
COPY ./autort /autort


COPY ./target/waterlooms-1.0-SNAPSHOT-jar-with-dependencies.jar /waterlooms/waterlooms.jar
VOLUME /waterlooms
# FOR PROD


# Compile Dependencies in order to run Pythonic Code
RUN apt-get update && \
	apt-get install -y maven && \
	pip install python-dotenv h5py tensorflow-gpu==2.9.0 scikit-learn==0.22.2.post1


RUN pip install pandas
RUN pip install --no-cache-dir numpy matplotlib psutil



# Export Python executables to the PATH
ENV PATH=/waterlooms/src/main/python:cd$PATH
ENV RUN_WITHIN_CONTAINER=TRUE

# Prevent asking for interactive status when installing the following package
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get -y install postgresql

#RUN apt install dirmngr gnupg apt-transport-https ca-certificates software-properties-common
# note the proxy for gpg
#for R
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9


# Now install R and littler, and create a link for littler in /usr/local/bin
# Also set a default CRAN repo, and make sure littler knows about it too
RUN apt-get update \
	&& apt-get install -y --no-install-recommends \
		littler \
                r-cran-littler \
		r-base \
                r-base-core \
		r-base-dev \
		r-recommended \
        && echo 'options(repos = c(CRAN = "https://cloud.r-project.org/"), download.file.method = "libcurl")' >> /etc/R/Rprofile.site \
        && echo 'source("/etc/R/Rprofile.site")' >> /etc/littler.r \
	&& ln -s /usr/share/doc/littler/examples/install.r /usr/local/bin/install.r \
	&& ln -s /usr/share/doc/littler/examples/install2.r /usr/local/bin/install2.r \
	&& ln -s /usr/share/doc/littler/examples/installGithub.r /usr/local/bin/installGithub.r \
	&& ln -s /usr/share/doc/littler/examples/testInstalled.r /usr/local/bin/testInstalled.r \
	&& install.r docopt \
	&& rm -rf /tmp/downloaded_packages/ /tmp/*.rds \
	&& rm -rf /var/lib/apt/lists/*



RUN Rscript -e 'install.packages(c("Rcpp", "reticulate", "data.table", "tensorflow", "tfdatasets", "keras", "stringi","stringr"))'

#should install python virtual environment first

RUN : \
    && apt-get update && apt install -y dbus\
    && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        software-properties-common \
    && add-apt-repository -y ppa:deadsnakes \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        python3.8-venv \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && :

RUN Rscript /waterlooms/src/main/R/initial.R



# Revert to root for execution
USER root

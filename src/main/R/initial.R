install.packages("reticulate")
library(reticulate)
virtualenv_install(envname="/opt/venv/greta",
                   packages = c("tensorflow==2.9.0",
                                "tensorflow_probability==0.7.0") )

use_virtualenv("/opt/venv/greta")

install.packages("tensorflow")
tensorflow::install_tensorflow(version="2.9.0-gpu", extra_packages="tensorflow-probability==0.7.0")
tensorflow::tf_version()
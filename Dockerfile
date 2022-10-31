# For finding latest versions of the base image see
# https://github.com/SwissDataScienceCenter/renkulab-docker
FROM carloferrigno/nustar-pipeline:0.1.0

# Uncomment and adapt if code is to be included in the image
# COPY src /code/src

# Uncomment and adapt if your R or python packages require extra linux (ubuntu) software
# e.g. the following installs apt-utils and vim; each pkg on its own line, all lines
# except for the last end with backslash '\' to continue the RUN line
#

# RUN apt-get update && \
#    apt-get install -y --no-install-recommends \
#    apt-utils \
#    vim

# install the python dependencies
USER root
RUN pip install pipupgrade && \
    pipupgrade --verbose --latest --yes

COPY requirements.txt  /tmp/
RUN pip install -r /tmp/requirements.txt

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    apt-utils nodejs python3-novnc tightvncserver libvncserver1 

#ADD js9.ipynb /home/jovyan/js9.ipynb

RUN mv /home/heasoft/js* /home/jovyan && \
    mv /home/heasoft/.jupyter /home/jovyan && \
    cat /home/jovyan/.jupyter/jupyter_notebook_config.py | sed "s/heasoft/jovyan/g" > /tmp/jupyter_notebook_config.py && \
    mv /tmp/jupyter_notebook_config.py /home/jovyan/.jupyter/jupyter_notebook_config.py && \
    chown -R jovyan /home/jovyan

RUN  usermod  --uid 1000 jovyan; usermod -g heasoft jovyan
USER heasoft
RUN cd /home/heasoft;git clone https://gitlab.astro.unige.ch/ferrigno/xspec_lmod.git;cd xspec_lmod;initpackage bwmod model.dat `pwd`;hmake;cd /home/jovyan
USER root
ADD .xspec/Xspec.init /home/jovyan/.xspec/Xspec.init
RUN chown -R jovyan /home/jovyan/.xspec
USER jovyan
# RENKU_VERSION determines the version of the renku CLI
# that will be used in this image. To find the latest version,
# visit https://pypi.org/project/renku/#history.
ARG RENKU_VERSION=1.7.1

########################################################
# Do not edit this section and do not add anything below

# Install renku from pypi or from github if it's a dev version
RUN if [ -n "$RENKU_VERSION" ] ; then \
        source .renku/venv/bin/activate ; \
        currentversion=$(renku --version) ; \
        if [ "$RENKU_VERSION" != "$currentversion" ] ; then \
            pip uninstall renku -y ; \
            gitversion=$(echo "$RENKU_VERSION" | sed -n "s/^[[:digit:]]\+\.[[:digit:]]\+\.[[:digit:]]\+\(rc[[:digit:]]\+\)*\(\.dev[[:digit:]]\+\)*\(+g\([a-f0-9]\+\)\)*\(+dirty\)*$/\4/p") ; \
            if [ -n "$gitversion" ] ; then \
                pip install --force "git+https://github.com/SwissDataScienceCenter/renku-python.git@$gitversion" ;\
            else \
                pip install renku==${RENKU_VERSION} ;\
            fi \
        fi \
    fi

#########################################################

USER root
RUN pip install mistune --upgrade
RUN pip uninstall -y jupyter_nbextensions_configurator jupyter_contrib_nbextensions
USER jovyan

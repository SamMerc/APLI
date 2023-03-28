# AP-Lab-I-2022 : Spectral and timing variability of magnetized neutron stars
Supervisor: Carlo Ferrigno

Magnetized neutron stars accreting from a companion star emit X-ray radiation that is strongly modulated with the spin phase. They are also variables on time scales varying from seconds to months.
By analyzing the pulsed emission, it is possible to study the emission mechanism and its geometry in the vicinity of the neutron star, where extreme gravity and radiation conditions are present.
By analyzing the variability on different time scales, it is possible to study how the neutron star interact with the surrounding medium to capture the material to be accreted.
The student will be asked to familiarize with one or two X-ray facilities and their data format as well as data reduction pipelines. Then, specific tools used for spectral and timing analysis 
will be applied to the reduced data sets and results visualized using python notebooks. Note that many methods are common to several research fields and can be used beyond the 
high-energy-astrophysics domain. Even if the student will be asked to contribute to the development of parts of the analysis and visualization tools using standard software versioning, 
most of the technical details will be mitigated by using software containers and self-developed python packages. 
The student will be asked to read and summarize research papers and provide a written report of their work that includes a literature review.


## Use with git
```bash
git clone git@renkulab.io:carlo.ferrigno/ap-lab-i-2022.git 
make build
JUPYTER_PORT=4444 make notebook
```
or if you need the command line
```bash
make run
```

Then
```bash
git commit "MYFILE" -m "MYMESSAGE"
git push
```

To use a graphical interface see 
`https://gist.github.com/sorny/969fe55d85c9b0035b0109a31cbcb088`
(The following has a wrong command)
`https://www.isdc.unige.ch/integral/download/osa/doc/11.2/osa_inst_guide/node9.html#SECTION00061100000000000000`

## Use in Renku (need testing)
Run it with enough storage and 8 Gb memory, better at least 1 CPU

-------------------

## Introduction

This is a Renku project - basically a git repository with some
bells and whistles. You'll find we have already created some
useful things like `data` and `notebooks` directories and
a `Dockerfile`.

## Working with the project

The simplest way to start your project is right from the Renku
platform - just click on the `Environments` tab and start a new session.
This will start an interactive environment right in your browser.

To work with the project anywhere outside the Renku platform,
click the `Settings` tab where you will find the
git repo URLs - use `git` to clone the project on whichever machine you want.

### Changing interactive environment dependencies

Initially we install a very minimal set of packages to keep the images small.
However, you can add python and conda packages in `requirements.txt` and
`environment.yml` to your heart's content. If you need more fine-grained
control over your environment, please see [the documentation](https://renku.readthedocs.io/en/latest/user/advanced_interfaces.html#dockerfile-modifications).

## Project configuration

Project options can be found in `.renku/renku.ini`. In this
project there is currently only one option, which specifies
the default type of environment to open, in this case `/lab` for
JupyterLab. You may also choose `/tree` to get to the "classic" Jupyter
interface.

## Moving forward

Once you feel at home with your project, we recommend that you replace
this README file with your own project documentation! Happy data wrangling!

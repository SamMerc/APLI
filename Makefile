IMAGE?=carloferrigno/aplab:0.1.0
IMAGE_LATEST?=carloferrigno/aplab:latest
DUSER := $(shell id -u)
JUPYTER_PORT?=1235
NBARGS_PYDICT?=''
time=$(shell date --utc +'%Y-%M-%dT%H:%m:%S')

push: build
	docker push $(IMAGE) 
	docker push $(IMAGE_LATEST) 

build: Dockerfile
	docker build . -t $(IMAGE) 
	docker build . -t $(IMAGE_LATEST) 

pull:
	docker pull $(IMAGE) 
	docker pull $(IMAGE_LATEST) 

test: build
	docker run -e HOME_OVERRRIDE=/tmp-home -v $(PWD):/tmp-home -it $(IMAGE_LATEST) bash -c 'source /usr/local/bin/init.sh; bash $$HOME/self-test.sh'
 
notebook: 
	docker run -e ODA_TOKEN=${ODA_TOKEN} -v $(PWD):/home/jovyan -it --entrypoint='' -v ${HOME}/.Xauthority:/home/jovyan/.Xauthority:rw --net=host --user $(DUSER)  $(IMAGE_LATEST) bash -c "source /init.sh;jupyter notebook --ip 0.0.0.0 --no-browser --port=$(JUPYTER_PORT)"

run: 
	docker run -e ODA_TOKEN=${ODA_TOKEN} -e DISPLAY=${DISPLAY} -v /etc/passwd:/etc/passwd -v /tmp/.X11-unix:/tmp/.X11-unix --net=host -v ${HOME}/.Xauthority:/home/jovyan/.Xauthority:rw -v $(PWD):/home/jovyan -it --entrypoint='' --user $(DUSER)  $(IMAGE_LATEST) bash 
# -p $(JUPYTER_PORT):$(JUPYTER_PORT) 
squash: build
	docker build . -t $(IMAGE)-noentry -f Dockerfile-noentry
	docker-squash -t $(IMAGE)-squashed $(IMAGE)-noentry

singularity: squash
	docker run -v /var/run/docker.sock:/var/run/docker.sock -v ${PWD}:/output --privileged -t --rm quay.io/singularity/docker2singularity $(IMAGE)-squashed


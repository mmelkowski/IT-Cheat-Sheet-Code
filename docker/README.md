# Quick-Docker

As Wikipedia say:

> *Docker is a set of platform as a service (PaaS) products that use OS-level virtualization to deliver software in packages called containers.*

TLDR: An environnement in a can

This folder aim to give a basic knowledge of Docker usage for environement virtualisation

## Docker: How to use it

The principle of Docker can be resume as writting a dockerfile which serve as an instruction to build an image.

This image can then be used to start a container which will contain all information and program installed by the `dockerfile`.

## 1) Dockerfile : The container's DNA

In order to build an image the user will have to write a `Dockerfile` as instruction for the `DockerEngine`.

This file is a plain text using the keyword in the table below

| Keyword | Definition |
|---|---|
| FROM | Base the image on a previous image or OS |
| COPY | Copy files or folder inside the container |
| RUN | Execute a command during image building |
| CMD | Define the default command that is executed once you start the container |
| ENTRYPOINT | Set the default program to execute, upon running the containner argument will be passed to the program set as entrypoint |

These are the major keyword use. A complete list can be found in the [official Docker documentation](https://docs.docker.com/).

### Example of Dockerfile

```dockerfile
FROM mambaorg/micromamba:latest

COPY environment.yml /tmp/conda-tmp/environment.yml

RUN /bin/micromamba install -f /tmp/conda-tmp/environment.yml
```

Here the docker file will be based on the tiny micromamba release and will install all conda package found in the yml file.

### Docker image creation (docker build)

```bash
docker build --file <dockerfile> --tag <image>:<version> <context_path>
```

The context path is the location from where the image will be build.

## 2) From image to container

Now that an image is build we can start a container using the command `docker run`.

```bash
docker run <image>:<version>
```

Many more options are available such as:

| option | Definition |
|---|---|
| --name | Container's name |
| -it | Use interactive mode (-it is short for --interactive + --tty) |
| --rm | Remove the container after use (when it's done, it's gone) |
| -v | Easy way to mount volume from outside to the container (see --mount for more mounting option) |
| --entrypoint | Overide or the the Dockerfile ENTRYPOINT value |
| -c | Execute a command inside the contianer using the entrypoint program |

Example:

```bash
docker run --name "myContainer" -it -v /path/to/workfolder:/workfolder --entrypoint "/bin/bash" <image>:<version>
(in container) $ ls -lh /workfolder

docker stop "myContainer"
docker rm -v "myContainer"
```

As you can starting a container is great but stopping it is better. The command `docker stop` and `docker rm` will help you get rid of those pesky background container.

If the container is suppose to have a one time use its even better to directly use the `--rm` option.

If you just need to run a single command you can directly pass it to the containner:

```bash
docker run --name "myContainer" --user 1000:1000 -v /path/to/workfolder:/workfolder --entrypoint "/bin/bash" <image>:<version> -c "ls -lh /workfolder"
```

## 4) Unusual building (docker commit)

Assuming a container was edited you might want to save it as a new image by using the `docker commit` command:

```bash
docker commit <container_id> <image>:<version>
```

You can even edit the docker image entrypoint:

```bash
docker commit --change='ENTRYPOINT ["/src/foo.sh"]' <container_id> <image>:<version>
```

## 5) Image management: sharing is caring

An image can be save locally to be transfered on remote server or friendly co-worker.

### Save

```bash
docker save --output image_docker.tar <image>:<version>
```

### Load

```bash
docker load --input image_docker.tar
```

### Remove

If you generated a wrong image or its no longer in use it can be removed using `docker remove`

```bash
docker image remove <image>:<version>
```

## 6) Testing: practice makes perfect

in the folder `test` you can find the dockerfile presented above:

### → Build-it

```bash
cd /path/to/docker/test
docker image build --file Dockerfile --tag mydockertest:1.0.0 .
```

### → Run-it

```bash
docker run --rm mydockertest:1.0.0
```

### → Delete-it

```bash
docker image remove mydockertest:1.0.0 
```

## 7) For more on youtube

### A short explenation by fireship.io

https://www.youtube.com/watch?v=gAkwW2tuIqE

### A longer explenation by fireship.io

https://www.youtube.com/watch?v=Gjnup-PuquQ
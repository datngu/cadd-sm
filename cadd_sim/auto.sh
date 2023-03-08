## BUILD DOCKER AUTOMATICLY
docker build -t cadd:v0.0.0 -f Dockerfile .

# v0.0.0
docker tag cadd:v0.0.0 ndatth/cadd:v0.0.0
docker push ndatth/cadd:v0.0.0
echo DONE


### test docker

docker run -it --rm -v /sigma4:/sigma4 --name cadd ndatth/cadd:v0.0.0
docker start cadd
docker attach cadd
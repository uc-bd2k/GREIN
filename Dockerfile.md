How to build
============
docker build . -t shiny.ilincs.org:5000/grin:v5
#docker build . -t shiny.ilincs.org:5000/grin:dev

#38384: dev version; 38385: production version

How to run
==========

docker run -d -t -p 38385:3838  -v /mnt/raid/tmp/iLincs/signatures/:/mnt/raid/tmp/iLincs/signatures/ -v /opt/raid10/:/opt/raid10/ -v /opt/raid10/genomics/naim/Thesis/geo/:/opt/raid10/genomics/naim/Thesis/geo/ -v /opt/raid10/genomics/naim/myGithub/GRIN/:/srv/shiny-server shiny.ilincs.org:5000/grin:v5

#docker run -d -t -p 38384:3838  -v /mnt/raid/tmp/iLincs/signatures/:/mnt/raid/tmp/iLincs/signatures/ -v /opt/raid10/:/opt/raid10/ -v /opt/raid10/genomics/naim/Thesis/geo/:/opt/raid10/genomics/naim/Thesis/geo/ -v /opt/raid10/genomics/naim/myGithub/GRIN_dev/:/srv/shiny-server shiny.ilincs.org:5000/grin:dev

How to push to repository
==========================

docker login -u ***** -p ***** shiny.ilincs.org:5000
docker tag grin shiny.ilincs.org:5000/ilincs/grin
docker push shiny.ilincs.org:5000/ilincs/grin

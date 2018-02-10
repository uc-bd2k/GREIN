How to build
============
docker build . -t shiny.ilincs.org:5000/grsn:v2


How to run
==========

docker run -d -t -p 38384:3838 -v /opt/raid10/:/opt/raid10/ -v /opt/raid10/genomics/naim/Thesis/geo/:/opt/raid10/genomics/naim/Thesis/geo/ -v /opt/raid10/genomics/naim/myGithub/GRSN/:/srv/shiny-server shiny.ilincs.org:5000/grsn:v2


How to push to repository
==========================

docker login -u ***** -p ***** shiny.ilincs.org:5000
docker tag grsn shiny.ilincs.org:5000/ilincs/grsn
docker push shiny.ilincs.org:5000/ilincs/grsn

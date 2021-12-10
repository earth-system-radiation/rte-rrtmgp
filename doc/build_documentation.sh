#!/bin/sh

set -e

rm -rf ../public
ford ford_site/extensions-front-matter.md
ford ford_site/kernels-front-matter.md
cd jekyll_site/
bundle exec jekyll build
cd _site/
if [ "$1" != "-ci" ]; then grep -rl --include \*.html /assets/main.css | xargs sed -i 's#/rte-rrtmgp/assets/main.css#assets/main.css#g'; fi
cd ..
mv _site/* ../../public

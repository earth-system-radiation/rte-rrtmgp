#!/bin/sh
rm -rf ../../public
ford ford_site/index.md
cd jekyll_site/
bundle exec jekyll build
cd _site/
grep -rl --include \*.html /assets/main.css | xargs sed -i 's#/rte-rrtmgp/assets/main.css#assets/main.css#g'
cd ..
cp -a _site/. ../../public

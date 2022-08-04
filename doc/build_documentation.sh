#!/bin/sh
rm -rf ../../public
mkdir ../../public
ford ford_site/rrtmgp-fortran-interface.md
ford ford_site/rrtmgp-kernels.md
ford ford_site/rte-fortran-interface.md
ford ford_site/rte-kernels.md
cd jekyll_site/
bundle exec jekyll build
cd _site/
if [ `uname` = "Linux" ]; then
    if [ "$1" != "-ci" ]; then grep -rl --include \*.html /assets/main.css | xargs sed -i 's#/rte-rrtmgp/assets/main.css#assets/main.css#g'; fi
else
    if [ "$1" != "-ci" ]; then grep -rl --include \*.html /assets/main.css | xargs sed -i '' 's#/rte-rrtmgp/assets/main.css#assets/main.css#g'; fi
fi
cd ..
cp -a _site/* ../../public

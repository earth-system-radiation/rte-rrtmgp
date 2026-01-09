#!/bin/sh
rm -rf ../../public
mkdir ../../public
#
# Ford auto-documentation file
#
ford ford_templates/rrtmgp-fortran-interface.md
ford ford_templates/rrtmgp-kernels.md
ford ford_templates/rte-fortran-interface.md
ford ford_templates/rte-kernels.md
ford ford_templates/ssm.md

#
#
#
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

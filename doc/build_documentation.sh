#!/bin/sh

# ford ford_site/extensions-front-matter.md 
# ford ford_site/kernels-front-matter.md 
cd jekyll_site/
bundle exec jekyll build
cd _site/
grep -rl --include \*.html /assets/main.css | xargs sed -i 's#/rte-rrtmgp/assets/main.css#assets/main.css#g'
grep -rl --include \*.html /assets/main.css | xargs sed -i 's#/assets/main.css#assets/main.css#g'
cd ..
mkdir ../../public/jekyll_site
cp -a _site/. ../../public/jekyll_site
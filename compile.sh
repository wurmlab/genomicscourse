#!/bin/bash
# This is used by Travis-CI to compile Markdown to html for gh-pages

shopt -s extglob # extend bash regex

echo "Copying repo content to out folder"
mkdir -p out > /dev/null 2>&1
cp -rf !(out) out/

cd out

# Add this when practical has submodules (e.g. MSMC 2016-SIB)
# echo "Removing submodules git"
# find . -mindepth 2 -name .git | xargs rm -rf
# rm .gitmodules
# sed -ie '/submodule/,+1d' .git/config

echo "Converting markdown to html"
find . -name "*.md" -type f -print0 | \
  # pv -0 | \
  xargs -0 -I{} \
  sh -c '
    repodir=$(pwd)
    dir=$(dirname $1)
    base=$(basename $1)
    name=${base%.*}
    ext=${base##*.}
    cd $dir
    pandoc -s \
      -f markdown_github+yaml_metadata_block \
      -c ${repodir}/css/github-pandoc.css \
      -A ${repodir}/footer.html \
      --self-contained \
      -i "${base}" \
      -o "${name}.html"
    cd $repodir
  ' -- {}
cd ..

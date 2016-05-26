#!/bin/bash
# This is used by Travis-CI to compile Markdown to html for gh-pages

shopt -s extglob # extend bash regex

echo "Copying genomicscourse content to website"
mkdir -p out > /dev/null 2>&1
cp -r !(out|*.sh) out/

echo "Converting markdown to html"
cd out
find . -name "*.md" -type f -print0 | \
  # pv -0 | \
  xargs -0 -I{} \
  sh -c 'dir=$(dirname $1); base=$(basename $1); name=${base%.*}; ext=${base##*.}; \
    pandoc -i "$1" -o "${dir}/${name}.html" --from markdown-yaml_metadata_block' -- {}
cd ..

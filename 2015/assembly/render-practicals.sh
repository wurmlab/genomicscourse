#!/bin/bash
pandoc -s README.md -o README.html
pandoc -s assembly-practical-part2.md -o assembly-practical-part2.html
pandoc -s assembly-practical-part3.md -o assembly-practical-part3.html
pandoc -s assembly-practical-extra-assemblers.md -o assembly-practical-extra-assemblers.html
pandoc -s assembly-practical-extra-qa.md -o assembly-practical-extra-qa.html

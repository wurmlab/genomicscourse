# Genome-Bioinformatics Practicals

<!-- Updated by Paolo Inglese, 2022 -->

Course material for the practicals organised
   * for the MSc Bioinformatics at QMUL.
   * for the 2012 SIB summer school in bioinformatics and population genomics (Adelboden)
   * for the 2016 SIB spring school in bioinformatics and population genomics (Leukerbad)

Practicals include:
* Genome assembly and annotation
* Read mapping and variant calling
* Simple gene expression analysis

Previous years courses are backed up in the `previous-courses` directory, while
the current course is saved in the directory `current-year`. Remember to backup
it as `202x` in the `previous-courses` when starting a new year.

# Website generation

Website is generated from Mardowns using `Jekyll`. Be sure that the file
`_config.yml` is updated to point to the current year course materials.  
The links in the Markdowns must point to the `html`, not the `md` files.

# Deploying

I compile, then cp out to gh-pages


int3dInterfaces
===============

### Description

Pipeline to extract protein interfaces for all the available Interactome3d structures in a given organism. The project automatically runs the following steps:

1. Downloads all the necessary files for a given organism
2. Splits the interacting partners into individual pdbs
3. Runs NACCESS to calculate the accessbilities for the individual proteins as well as for the protein complexes.
4. Calculates the difference on relative accessibility from complexes and monomers to detect residues changing their accessibiliy.

### Output files

The pipeline prints the following files on the `results` directory:

* `accessibilities.tab`. Contains all the relative accesibilities for all the residues present on the pdbs.
* `interfaces.tab`. Contains all the interface residues for each chain in each complex.

The remaining columns in the files contain additional information regarding the residues or the reliability of the protein structure (i.e. Structure, Model or Domain-Domain Model).    

### Dependencies

It requires [naccess](http://www.bioinf.manchester.ac.uk/naccess/).

### How to run it

To run the pipeline, these are the two necessary commands.

      make prepare
      make all

### Clean the project

If you want to run everything from scratch, clean the whole project running:

    make clean

##Other organisms

So far, the pipeline has not been adapted to run different organisms. By default, it runs yeast. You can select other organisms by updating the content of the variables `INT3URL` and `FILELISTPAGE` on the `Makefile`. It's as easy as substite `yeast` by the organism of interest. Remember to `make clean` before running a new organism.

The list of current organisms includes:

* Arabidopsis thaliana (athaliana)
* Bacillus subtilis (bsubtilis)
* Canorhabditis elegans (celegans)
* Campylobacter jejuni (cjejuni)
* Drosophila melanogaster (fly)
* Escherichia coli (ecoli)
* Helicobacter pylori (hpylori)
* Homo sapiens (human)
* Mus musculus (mouse)
* Mycobacterium tuberculosis (mtuberculosis)
* Mycoplasma pneumoniae (mpneumoniae)
* Plasmodium falciparum (pfalciparum)
* Rattus norvegicus (rat)
* Saccharomyces cerecisiae (yeast)
* Schizosaccharomyces pombe (spombe)
* Treponema pallidum (tpallidum)

# ActiveDriverWGS 1.1.2

* Bug fixed: datasets with very few unmutated elements (<10) previously failed to include these elements in results.
* New option: mitochondrial mutations can be analyzed (chrM). 

# ActiveDriverWGS 1.1.1

* Bug fixed: a rare set of very large elements with depletion of mutations would sometimes appear as enriched in mutations. 

# ActiveDriverWGS 1.1.0

## Major Changes

* Alternative reference genomes for human (hg19, hg38) and mouse (mm9, mm10) are now supported. Default option is hg19.
* Full rewrite of ADWGS_rest for memory and speed efficiency.
* Improved handling on indels on boundaries of elements is now default. Removed parameter element_bias that controlled this behaviour earlier. 
* README.rd and front page now include a simple example of running ActiveDriverWGS, and the vignette has been updated.
* Site IDs need to match element IDs. Example databases with PTM sites have been updated. 

# ActiveDriverWGS 1.2.0

* Updated trinucleotide quantification in sequence and mutations. Trinucleotides on the inverse strand are now reverse complemented instead of complemented to match COSMIC signatures.
* More accurate and stable numeric estimates of expected mutations are provided as output instead of the Poisson-sampled estimates provided previously. 
* The package now allows identification of genomic elements with significantly fewer mutations than expected, potentially reflecting negative selection of elements. This is enabled using an optional command line parameter detect_depleted_mutations. In the analysis depleted mutations, elements with enriched mutations are assigned non-significant P-values.

# ActiveDriverWGS 1.1.2

* Bug fixed: datasets with very few unmutated elements (<10) previously failed to include these elements in results.

# ActiveDriverWGS 1.1.1

* Bug fixed: a rare set of very large elements with depletion of mutations would sometimes appear as enriched in mutations. 

# ActiveDriverWGS 1.1.0

## Major Changes

* Alternative reference genomes for human (hg19, hg38) and mouse (mm9, mm10) are now supported. Default option is hg19.
* Full rewrite of ADWGS_rest for memory and speed efficiency.
* Improved handling on indels on boundaries of elements is now default. Removed parameter element_bias that controlled this behaviour earlier. 
* README.rd and front page now include a simple example of running ActiveDriverWGS, and the vignette has been updated.
* Site IDs need to match element IDs. Example databases with PTM sites have been updated. 
# ActiveDriverWGS 1.1.0

## Major Changes

* Alternative reference genomes for human (hg19, hg38) and mouse (mm9, mm10) are now supported. Default option is hg19.
* Full rewrite of ADWGS_rest for memory and speed efficiency.
* Improved handling on indels on boundaries of elements is now default. Removed parameter element_bias that controlled this behaviour earlier. 
* README.rd and front page now include a simple example of running ActiveDriverWGS, and the vignette has been updated.
* Site IDs need to match element IDs. Example databases with PTM sites have been updated. 
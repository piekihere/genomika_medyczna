![Logo](images/logo.png)
usage: Varianinator [-h] [-i INPUT] [--id ID] [-o OUTPUT] [--show-na] [--rare] [--pathogenic] \
 \
Search ClinVar, DBSNP and SnpEff for described variant in VCF file or by variant ID. \
 \
options: \
  -h, --help &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp; &emsp; Show this help message and exit \
  -i INPUT, --input INPUT &emsp;&emsp;&emsp; Path to VCF file \
  --id ID&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;Search for variants by ID. Possible formats: \
  &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp; &emsp; &emsp;*RSid e.g. rs58991260 \
  &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp; &emsp; &emsp;*specific SNP e.g. chr1:g.35367G>A \
  &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp; &emsp; &emsp;*ENSBML gene ID e.g. ENSG00000113368 \
  -o OUTPUT, --output OUTPUT&ensp;Path to output file. \
  --show-na&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp; &emsp;Show variants even if there is no information for them in the databases \
  --rare &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;Save only rare variants \
  --pathogenic &emsp;&emsp;&emsp;&emsp; &emsp; &emsp; &emsp;Save only pathogenic variants. 

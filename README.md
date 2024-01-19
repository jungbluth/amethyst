![image](https://github.com/rckarns8/amethyst/assets/42095420/f9da1bd9-250e-4607-8512-a125088eabed)


Written 2023, by Rachael Storo. Contributions by Giles Goetz. Based on a workflow written by Nastassia Patin.
v. 1.0
# Important notes:
- Amethyst requires python >3.6 
- Amethyst requires data files to be located in amethyst/00_data/fastq, separated in R1 for forward reads, R2 for reverse.

To get started, please see the [Quick Start Guide](https://github.com/rckarns8/amethyst/wiki/0.0-Quick-Start-Guide)


# Roadmap:
Note- this is a prioritized list of features to add to Amethyst. This is subject to change.
| Program      | Function                          | Status      |
|--------------|-----------------------------------|-------------|
| fastqc       | read QC                           | complete    |
| multiqc      | read QC                           | complete    |
| multitrim    | read trimming                     | complete    |
| megahit      | assembly                          | complete    |
| bbnorm       | normalization                     | ON HOLD     |
| bowtie2      | assembly coverage analysis        | complete    |
| prodigal     | gene prediction                   | complete    |
| prokka       | gene annotation                   | complete    |
| sourmash     | diversity analysis- assembly based| complete    |
| maxbin2      | binning                           | complete    |
| checkm       | bin/MAG quality                   | complete    |
| dRep         | dereplication of bins/MAGs        | tested      |
| GTDB-Tk      | taxonomic assignment of bins/MAGs | written     |
| nonpareil    | read-based coverage analysis      | backlog     |
| sourmash     | diversity analysis- read based    | backlog     |
| humann       | functional analysis               | backlog     |


*BBnorm is on hold as normalization is not ideal for binning, and it was having errors.





Features to add in new versions: 
- co-assembly capabilities
- overarching rules to combine pseudo rules
- benchmarking/logs- coded in the current version but not writing
- init file to make all conda environments needed to run amethyst
- yaml files for conda environments needed to run amethyst
- choice of assembler
- Docker Containers 

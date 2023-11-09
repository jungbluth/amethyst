Written 2023, by Rachael Storo. Based on a workflow written by Nastassia Patin.
v. 1.0
# Important notes:
- Amethyst requires python 3.6> 
- Amethyst requires data files to be located in amethyst/00_data/fastq, separated in R1 for forward reads, R2 for reverse.
- to create the multitrim environment, navigate to amethyst directory, which contains the relevant yaml file, before running the commands found in the [Quick Start Guide](https://github.com/rckarns8/amethyst/wiki/0.-Quick-Start-Guide)



To get started, please see the [Quick Start Guide](https://github.com/rckarns8/amethyst/wiki/0.-Quick-Start-Guide)



Features to add in new versions: 
- co-assembly capabilities
- overarching rules to combine pseudo rules
- benchmarking
- init file to make all conda environments needed to run amethyst
- yaml files for conda environments needed to run amethyst
- fix issue where adapt_trim is overwriting temp files for each sample
- choice of assembler
- read-based analyses and gene annotation
- Binning for MAGS
- Docker Containers 

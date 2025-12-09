# oliver-yanco-et-al-2026

Data analysis for publication: "Interacting effects of human presence and landscape modification on birds and mammals"

### Abstract

Sustainable human-wildlife coexistence requires a mechanistic understanding of the many ways that humans affect animals. However, progress is hampered by the lack of accessible data measuring the dynamic presence of people. Here, we leverage mobile-device data to disentangle how human presence and landscape modification differentially influence the use of geographic and environmental space for 37 mammal and bird species across the United States. Human presence affected over 65% of species, with substantial variation across species. For ~60% of species that responded to human activities, the effects were interdependent – animals tended to react more strongly to human presence in less modified habitats. Our results demonstrate that human presence and landscape modification have complex combined effects on wildlife which need to be considered for effective management.


### Repository Structure

```
/oliver-yanco-et-al-2026 
  |  
  +--/src                 # Source code directory
  |   |
  |   +--config1.env      # config for repository filepaths and conda env
  |   |
  |   +--config2.env      # config for package libraries
  |   |
  |   +--/funs            # Source code for custom functions called by other scripts
  |   |   
  |   +--/hpc             # Scripts for submitting jobs to Slurm manager sequentially
  |   |
  |   +--/mosey           # Scripts for database annotation
  |   |
  |   +--/startup.R       # Source code for some basic environment configuration
  |   | 
  |   +--/workflow        # Scripts that execute elements of the workflow
  |       |
  |       +--part1_data_prep
  |       |
  |       +--part2_modeling
  |       |
  |       +--part3_figures
  |    
  +--/ctfs                # Control files for workflow scripts
  | 
  +--/conda_envs          # Stores .yml files with conda environment specifications
  |  
  +--/raw_data            # Raw data stored as initially received, including database
  |  
  +--/processed_data      # Processed data products
  |   |
  |   +--/intermediate_db_copies  # Working version of the database
  |   |   
  |   +--/safegraph       # Human mobility data used for part 1
  |   |  |
  |   |  +--/counties-dates-2-10-22-reformatted
  |   |     |
  |   |     +--/daily-data
  |   |
  +--/gee_data            # Used during environmental annotation steps in part 1
  |   |
  |   +--/annotated
  |   |
  |   +--/csvs_gee_ingest
  |   |  
  |   |
  +--/out                 # Analytical outputs, interim products
  |    |
  |    +--/single_species_models  # Single species model .rdata files
  |    |   |
  |    |   +--/niche_interactive
  |    |   |
  |    |   +--/niche_additive
  |    |   |
  |    |   +--/area-interactive
  |    |   |
  |    |   +--/area_additive
  |    |
  |   +--/single_species_models_reruns
  |    |   |
  |    |   +--/niche_interactive
  |    |   |
  |    |   +--/niche_additive
  |    |   |
  |    |   +--/area_interactive
  |    |   |
  |    |   +--/area_additive
  |    |
  |    +--/model_diagnostics  # Single species model summary PDFs
  |    |   |
  |    |   +--/area
  |    |   |
  |    |   +--/niche
  |    |
  |    +--/model_diagnostics_reruns
  |    |   |
  |    |   +--/area
  |    |   |
  |    |   +--/niche
  |    |
  |    +--/safegraph_summary
  |    |
  |    +--/dbbmms         # Utilization distributions for individual-weeks
  |    |
  |    +--/event-annotation
  |    |
  |    +--/event-cbg-intersection
  |    |
  |    +--/covid_results  # Output CSVs from scripts in part3_model_effects.sh
  |    |
  |    +--/figures        # Output figures
  |    |
  |    +--/intra-ind-models  # Models for species with sufficient data in both 2019 & 2020

```

### Data Availability

The wildlife movement data that serves as input for part 1 of the workflow are archived publicly on the [Movebank Data Repository](https://www.movebank.org/cms/movebank-content/data-repository) for reproducibility. See the manuscript's supplementary table 1 for DOIs and dataset contacts. Select species data could not be made public due to conservation concerns. Human mobility data used in part 1 for one component of the anthropogenic annotation cannot be made publically available. The secondary data products used as input to parts 2 and 3 of the workflow are publicly available on OSF: [DOI 10.17605/OSF.IO/3UA2C](https://osf.io/3ua2c/). These tabular data products contain all species and individuals in the analysis with environmental and anthropogenic annotations derived from the animal GPS locations and timestamps. Users interested in reproducing the workflow are encouraged to start at part 2 with these secondary products downloaded from OSF and stored in the `/out` directory.

### Part 1: Data Prep

**R scripts:** `src/workflow/part1_data_prep`

**SLURM scripts:** `src/hpc/part1*.sh`

- Build database of wildlife movement data for 37 bird and mammal species pulled from [Movebank](https://www.movebank.org/cms/movebank-main) studies.
  - Data are stored as a [mosey_db](https://github.com/benscarlson/mosey_db), a SQLite relational database built to store data from Movebank.
  - [This repository release](https://github.com/julietcohen/mosey_db/releases/tag/v1.0.0) includes the forked `mosey_db` code used to build the animal movement database used as input for this repository's workflow part1.
- Subset wildlife movement data to region and time period of interest.
- Annotate database with environmental layers for temperature, NDVI, and elevation.
- Annotate database with human mobility using daily mobile device counts.
- Annotate database with landscape modification based on a multi-year aggregated metric of anthropogenic modification. 
- Filter database for analysis minimum criteria.
- Clean database, including removing outlier events.
- Estimate utilization distributions via dynamic Brownian bridge movement models for each individual-week combination.
- Estimate environmental niche sizes based on the pooled variance of multidimensional hypervolumes of the environmental conditions.
- Execute niche breadth sensitivity analysis.
- Calculate fix rate per species.

### Part 2: Modeling

Use tabular niche and space use estimations for each individual-week as input to species-specific Bayesian mixed effects models across all species.

**R scripts:** `src/workflow/part2_modeling`

**SLURM scripts:** `src/hpc/part2*.sh`

- Fit space use interactive and additive models.
- Fit niche interactive and additive models.
- Fit intra-individual interactive and additive models for individuals with sufficient data in both 2019 and 2020.
- Select interactive or additive models based on significance and produce model summaries.

### Part 3: Figures

**R scripts:** `src/workflow/part3_figures`

**SLURM scripts:** `src/hpc/part3*.sh`

- Summarize model effects from single species models.
- Wildlife responses to the major components of human activity across the United States (Fig 1)
- Interacting effects of human activities on wildlife’s use of geographic and environmental space (Fig 2)
- Plastic behavioral responses to human mobility (Fig 3)
- Combined impact of human activities on wildlife use of geographic and environmental space (Fig 4) 
- Relationship between weekly area size and weekly sample size (Fig S1)
- Plot niche breadth subsample sizes (Fig S2)
- Summarize sample sizes for area and niche single species models presented in Fig 2, paired with fix rate (Table S2)
- Posterior distributions of species-specific estimates (Figs S3 and S4):
  - effect of human mobility on area size 
  - interactive effect of human modification and human mobility on area size
  - effect of human mobility on niche size 
  - interactive effect of human modification and human mobility on niche size
- Distribution of census block group sizes (Fig S5)
- Distribution of utilization distribution sizes (Fig S6)
- Species-specific distributions of census block group sizes and utilization distribution sizes

### Development Environment

#### Python

This workflow was run with conda 24.11.1 and Python 3.12.8. With these installations on your machine, run the following commands in a terminal from the root of this repository to recreate our conda environment.

```
conda env create -f conda_envs/r_spatial2_direct_dependencies_environment.yml

conda activate r_spatial2
```

#### R

This workflow was run with R/4.3.1. Necessary packages with specified versions can be installed with [renv](https://rstudio.github.io/renv/articles/renv.html).

In order to activate the R env documented in this repo, clone the repo and run `renv::restore()`

For the few packages that failed to install within the archived `renv` environment, the user can install them into their R environment afterwards as needed. This may be the case in particular for part 3. Additionally, note that the versions of packages generated by `renv` are a mix of specific versions or ranges of appropriate versions.

#### Config and running locally versus on HPC

Update files `config1.env` and `config2.env` to include local filepaths prior to running workflow on a server. These are read in by the shell scripts to activate the conda environment, activate R, set working directory to the repository root, etc. 

The workflow is set up to run on a server using sequential SLURM jobs launched via shell scripts. Alternatively, users can run through all scripts locally for parts 2 and 3 by following the order of scripts and the arguments provided in the shell scripts. Note that the cores arguments specified in shell scripts exceed the amount available for most machines, so these values will need to be reduced. Reducing the cores will increase the run times. Any modeling steps can be skipped and instead the output models can be downloaded from the OSF project and used as input into sequential scripts.

### Contributing

Contacts:

- Ruth Oliver rutholiver@ucsb.edu
- Scott Yanco yancos@si.edu

We welcome feedback and questions. Please open an issue or create a fork of this repository.

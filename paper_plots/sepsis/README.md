# Analyze sepsis data

In this folder, we provide code for sequentially estimating the effect of IV fluid caps on 30-day mortality in sepsis-3 patients. The data are available in the Medical Information Mart for Intensive Care III (MIMIC-III) database, but the process of turning the MIMIC-III data into a downstream causal analysis is involved. Here, we outline this process.

## 1. Get access to [MIMIC-III](https://mimic.mit.edu/docs/iii/) in a SQL database

MIMIC-III is freely-available online, but in order to access it, you'll need to [create a PhysioNet account](https://physionet.org/login/?next=/settings/credentialing/) and [sign the data use agreement](https://physionet.org/sign-dua/mimiciii/1.4/).

Once you have the necessary permissions to access MIMIC-III, you can access it by setting up a [local PostgreSQL database](https://mimic.mit.edu/docs/gettingstarted/local/) or through [a cloud provider](https://mimic.mit.edu/docs/gettingstarted/cloud/) like AWS or Google BigQuery. We took the local PostgreSQL approach, but either one should be fine for what follows.

## 2. Create the sepsis-3 table

Alistair Johnson and Tom Pollard have created [some fantastic SQL scripts](https://github.com/alistairewj/sepsis3-mimic) for creating tables of sepsis patients. 
Specifically, running the `make-tables.sql` script from [this folder](https://github.com/alistairewj/sepsis3-mimic/tree/master/query) will generate a table called `sepsis3` which contains patients satisfying the [sepsis-3 criteria](https://jamanetwork.com/journals/jama/fullarticle/2492881). 

## 3. Extract the relevant data

Run our SQL query [`get_sepsis_patients.sql`](/paper_plots/sepsis/get_sepsis_patients.sql) found in this folder to extract sepsis-3  patients along with their 24-hour post-admission IV fluid intake, demographic information, and other medically relevant covariates. 

Save the resulting table as a csv to `paper_plots/sepsis/data/sepsis_patients.csv`.

## 4. Compute doubly robust confidence sequences

Now that we have the required sepsis patient data, we are ready to use our confidence sequences to sequentially estimate the effect of IV fluid caps on 30-day mortality. Simply run in your terminal:

```zsh
# Enter the sepsis directory
cd sequential.causal/paper_plots/sepsis

# Run the R script. Note this may take a while...
Rscript compute_sepsis_cs.R
```

The above script will generate .RData files which can be read in by the RMarkdown file [`paper_plots/plots_for_paper.Rmd`](/paper_plots/plots_for_paper.Rmd) which generates all of the paper's plots. For more information, see the [paper_plots](/paper_plots) folder.

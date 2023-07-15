## Set up the environment

Install [(mini)conda](https://conda.io/miniconda.html) as a light-weighted package management tool. Run the following commands to initialize and setup the conda environment for the program

```bash
# add channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# create conda environment
conda create --name prog-env

# activate conda environment
conda activate prog-env

conda install -c bioconda -c conda-forge python=3 minimap2 numpy pandas seaborn matplotlib
```

After successfully setup the environment and dependencies, clone the program into your desirable place.


### Running Program
```bash
main.py build -f<a|q> <k> <out_dir> <ref> <ref> ... <ref>
main.py query -f<a|q> <db_dir> <out_dir> <query>
main.py analysis <qry_dir>
```





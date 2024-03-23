# run-pipeline

## test databases formation

Databases of fastq results from Genorobotics exppeditions can be found under `testdbs/`. The species are known and each sub-folder contains the associated fastq file as well as the reference sequences for the four barcoding genes extracted from NCBI.

The databases are formed through the sorting of the expeditions' results folders by the [`create_test_db` notebook](create_test_db.ipynb). The expedition folders are added to the `input` folder before running the notebook


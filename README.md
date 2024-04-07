# run-pipeline

## test databases formation

Databases of fastq results from Genorobotics expeditions can be found under `testdbs/`. The species are known and each sub-folder contains the associated fastq file as well as the reference sequences for the sequenced genes extracted from NCBI.

The databases are formed through the sorting of the expeditions' results folders by the [create_test_db notebook](create_test_db.ipynb). The unsorted expedition folders are added to the `input` folder before running the notebook.

Two test databases are already provided, they are the results of two Genorobotics expeditions: the summer expedition used for the publication and the Jardin botanique expedition. The raw contents of the unsorted expedition folders are not provided in the `input` folder because of their size. Members can check them on the Genorobotics google drive.

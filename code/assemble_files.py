file_name = "aggregated/"+snakemake.wildcards['phenotype']+".runs.csv"

import csv
with open(file_name, 'wt',  newline= '') as out_file:
  csv_writer = csv.writer(out_file, delimiter=' ')
  csv_writer.writerows(snakemake.input)

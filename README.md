## Identifying regional relatedness of African-Americans and European-Americans across the United States

### General info
MATHEMATICA code to
  1. read, efficiently store, analyze, and visualize genomic relatedness of African-Americans and European-Americans across the United States
  2. read, store, model, and visualize relatedness of individuals across the United States based on the data from US Census Bureau (data avilable from [IPUMS](https://usa.ipums.org/usa/))
  3. assess the correlation between genomic data and census data
Results are available in the preprint article [The Great Migration and African-American genomic diversity](http://biorxiv.org/content/early/2015/10/15/029173), to appear in PLoS Genetics.

### NOTES
  - I have removed the hard-coded paths to the data files (hosted on our clusters) from this project for security reasons. So, this code, in its current state on GitHub, will not work! However, it is shared publicly to show the steps I took to solve the specific problems discussed in our publication.
  - The code is not the most elegant, and there are some parts that intentionally have not been refactored to keep different parts of the project physically separate from each other.

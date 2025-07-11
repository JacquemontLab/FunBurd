# FunBurd
![FunBurd.png](FunBurd.png)

## Description
Compute the functional burden of a biological function for a phenotype.

## How to install:

1. **Clone the Repository**:
   ```bash
   git clone https://github.com/JacquemontLab/FunBurd.git
   cd FunBurd
   ```

2. **Create a Virtual Environment (Optional but Recommended)**:
   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows use `venv\Scripts\activate`
   ```

3. **Install the Required Dependencies**:
   ```bash
   pip install -r requirements.txt
   ```

4. **Install the Module**:
   ```bash
   pip install .
   ```

## Dependencies

The module requires the following dependencies:
- Python 3.8 or higher
- pandas
- numpy
- scipy
- statsmodels

These dependencies will be automatically installed when you run `pip install -r requirements.txt`.

## Usage

After installation, you can import the module in your Python scripts:
```python
import FunBurd
```

## Functions

### Function Documentation: `FunBurd_one_gene_set`

#### Description
The `FunBurd_one_gene_set` function performs a burden analysis for a single gene set. It takes phenotype and variant data, collapses the variants into gene sets, and runs a regression model to assess the association between the gene burden and the phenotype. The function can handle both continuous and binary phenotypes and allows for the inclusion of covariates and interaction terms.

```python
FunBurd_one_gene_set(phenotype_table: pd.DataFrame,
                         phenotype_column_name: str,
                         variants_table: pd.DataFrame,
                         gene_list: list,
                         keep_all_regression_results: bool = False,
                         list_covariantes: list = None,
                         interaction_column: list = None,
                         column_for_conditional: str = None,
                         correction_outside_gene_set: dict = None)
```

#### Parameters
   Parameter | Type | Description | Default |
 |-----------|------|-------------|---------|
 | `phenotype_table` | `pd.DataFrame` | DataFrame containing the phenotypical data of each individual (See input example for mandatory columns). | |
 | `phenotype_column_name` | `str` | Name of the column in the phenotype table that contains the phenotype scores, can be either binary or continous scores, FunBurd will run accordingly. | |
 | `variants_table` | `pd.DataFrame` | DataFrame containing the variants altering each gene for each individual of the analysis. (See input example for mandatory columns) | |
 | `gene_list` | `list` | List of genes to collapse as a unique functional group and used for the burden analysis. | |
 | `keep_all_regression_results` | `bool` | If `False`, will only keep the results of the regression model for the gene burden. If `True`, will keep all the results of the regression model (including covariates and interactions terms). | `False` |
 | `list_covariantes` | `list` | List of covariates to include in the regression model in addition to the gene burden. | `None` |
 | `interaction_column` | `list` | List of columns for interaction terms in the regression model. | `None` |
 | `column_for_conditional` | `str` | Column name for conditional analysis. To use if the individuals are part of same families for example. | `None` |
 | `correction_outside_gene_set` | `dict` | To use if you want to correct the burden results based on the rest of the genome outside of the studied functional set of genes. Dictionary containing the categories of the genes in the genome for correction outside the gene set. (See input example for more details) | `None` |

#### Returns
 | Type | Description |
 |------|-------------|
 | `pd.DataFrame` | DataFrame containing the results of the burden analysis with each line representing a covariate association with the phenotype. |

### Function Documentation: `FunBurd_one_gene_set_wih_name`

#### Description
The `FunBurd_one_gene_set_wih_name` function performs a burden analysis for a single gene set with its name. It takes phenotype and variant data, collapses the variants into gene sets, and runs a regression model to assess the association between the gene burden and the phenotype. This function is similar to `FunBurd_one_gene_set` but includes the gene set name in the results.

```python
FunBurd_one_gene_set_wih_name(phenotype_table: pd.DataFrame,
                             phenotype_column_name: str,
                             variants_table: pd.DataFrame,
                             genes_item: list,
                             keep_all_regression_results: bool = False,
                             list_covariantes: list = None,
                             interaction_column: list = None,
                             column_for_conditional: str = None,
                             correction_outside_gene_set: dict = None)
```

#### Parameters
 | Parameter | Type | Description | Default |
 |-----------|------|-------------|---------|
 | `phenotype_table` | `pd.DataFrame` | DataFrame containing the phenotype data. | |
 | `phenotype_column_name` | `str` | Name of the column in the phenotype table that contains the phenotype scores (Can be either continuous or binary, the model will adapt accordingly). | |
 | `variants_table` | `pd.DataFrame` | DataFrame containing the variants altering each gene for each individual of the analysis. | |
 | `genes_item` | `list` | List containing the name of the gene set and the list of genes to collapse as a unique functional group for the burden analysis. The first element is the name of the gene set and the second element is the list of genes. | |
 | `keep_all_regression_results` | `bool` | If False, will only keep the results of the regression model for the gene burden. If True, will keep all the results of the regression model (including covariates and interactions). | `False` |
 | `list_covariantes` | `list` | List of covariates to include in the regression model in addition to the gene burden. | `None` |
 | `interaction_column` | `list` | List of columns for interaction terms in the regression model. | `None` |
 | `column_for_conditional` | `str` | Column name for conditional analysis. | `None` |
 | `correction_outside_gene_set` | `dict` | Dictionary containing the categories of the genes in the genome for correction outside the gene set. | `None` |

#### Returns
 | Type | Description |
 |------|-------------|
 | `pd.DataFrame` | DataFrame containing the results of the burden analysis with each line representing a covariate association with the phenotype for the gene set. |

### Function Documentation: `FunBurd_multiple_gene_sets`

#### Description
The `FunBurd_multiple_gene_sets` function performs a burden analysis for multiple gene sets in parallel. It takes phenotype and variant data, collapses the variants into gene sets, and runs a regression model to assess the association between the gene burden and the phenotype for each gene set. This function allows for parallel processing to speed up the analysis.

```python
FunBurd_multiple_gene_sets(phenotype_table: pd.DataFrame,
                          phenotype_column_name: str,
                          variants_table: pd.DataFrame,
                          gene_sets_dict: dict,
                          keep_all_regression_results: bool = False,
                          list_covariantes: list = None,
                          interaction_column: list = None,
                          column_for_conditional: str = None,
                          correction_outside_gene_set: dict = None,
                          n_cpus: int = 1)
```

#### Parameters
 | Parameter | Type | Description | Default |
 |-----------|------|-------------|---------|
 | `phenotype_table` | `pd.DataFrame` | DataFrame containing the phenotype data. | |
 | `phenotype_column_name` | `str` | Name of the column in the phenotype table that contains the phenotype scores (Can be either continuous or binary, the model will adapt accordingly). | |
 | `variants_table` | `pd.DataFrame` | DataFrame containing the variants altering each gene for each individual of the analysis. | |
 | `gene_sets_dict` | `dict` | Dictionary where keys are gene set names and values are lists of genes to collapse as unique functional groups for the burden analysis. | |
 | `keep_all_regression_results` | `bool` | If False, will only keep the results of the regression model for the gene burden. If True, will keep all the results of the regression model (including covariates and interactions). | `False` |
 | `list_covariantes` | `list` | List of covariates to include in the regression model in addition to the gene burden. | `None` |
 | `interaction_column` | `list` | List of columns for interaction terms in the regression model. | `None` |
 | `column_for_conditional` | `str` | Column name for conditional analysis. | `None` |
 | `correction_outside_gene_set` | `dict` | Dictionary containing the categories of the genes in the genome for correction outside the gene set. | `None` |
 | `n_cpus` | `int` | Number of CPUs to use for parallel processing. | `1` |

#### Returns
 | Type | Description |
 |------|-------------|
 | `pd.DataFrame` | DataFrame containing the results of the burden analysis with each line representing a covariate association with the phenotype for each gene set. |

### Function Documentation: `bootstraped_FunBurd`

#### Description
The `bootstraped_FunBurd` function performs a bootstrap burden analysis for a list of genes. It creates bootstrap samples of the gene list, performs a burden analysis for each sample, and returns the results. This function is useful for assessing the robustness of the burden analysis results.

```python
bootstraped_FunBurd(phenotype_table: pd.DataFrame,
                   phenotype_column_name: str,
                   variants_table: pd.DataFrame,
                   gene_list: list,
                   n_bootstrap: int = 1000,
                   bootstrap_size: int = None,
                   pvalue_threshold: float = 0.05,
                   return_effect_size_df: bool = False,
                   n_cpus: int = 1,
                   seed: int = None,
                   list_covariantes: list = None,
                   interaction_column: list = None,
                   column_for_conditional: str = None,
                   correction_outside_gene_set: dict = None)
```

#### Parameters
 | Parameter | Type | Description | Default |
 |-----------|------|-------------|---------|
 | `phenotype_table` | `pd.DataFrame` | DataFrame containing the phenotype data. | |
 | `phenotype_column_name` | `str` | Name of the column in the phenotype table that contains the phenotype scores (Can be either continuous or binary, the model will adapt accordingly). | |
 | `variants_table` | `pd.DataFrame` | DataFrame containing the variants altering each gene for each individual of the analysis. | |
 | `gene_list` | `list` | List of genes to collapse as a unique functional group and used for the burden analysis. | |
 | `n_bootstrap` | `int` | Number of bootstrap samples to create. | `1000` |
 | `bootstrap_size` | `int` | Size of each bootstrap sample. If None, will be set to half the size of the gene list. | `None` |
 | `pvalue_threshold` | `float` | P-value threshold for significance. | `0.05` |
 | `return_effect_size_df` | `bool` | If True, will return the DataFrame with the effect sizes of the gene sets. | `False` |
 | `n_cpus` | `int` | Number of CPUs to use for parallel processing. | `1` |
 | `seed` | `int` | Random seed for reproducibility. | `None` |
 | `list_covariantes` | `list` | List of covariates to include in the regression model in addition to the gene burden. | `None` |
 | `interaction_column` | `list` | List of columns for interaction terms in the regression model. | `None` |
 | `column_for_conditional` | `str` | Column name for conditional analysis. | `None` |
 | `correction_outside_gene_set` | `dict` | Dictionary containing the categories of the genes in the genome for correction outside the gene set. | `None` |

#### Returns
 | Type | Description |
 |------|-------------|
 | `dict` | A dictionary containing the results of the bootstrap burden analysis. If `return_effect_size_df` is True, the dictionary will also contain a DataFrame with the effect sizes of the gene sets. |
## Example of input data
## Input Data Format

### Genotyping DataFrame

The genotyping dataframe contains information about genetic variants for each individual. It must include the following columns:
   Column | Type | Description | Example |
 |--------|------|-------------|---------|
 | `SampleID` | str | Unique identifier for each sample/individual | "Sample0", "Sample1" |
 | `Gene` | str | Gene identifier (ENSG format recommended) | "ENSG00000132549" |
 | `Variant_Type` | str | Type of genetic variant | "Missense", "Splice_variant", "Stop_Gained", "Frameshift" |

**Example:**
 | SampleID | Gene | Variant_Type |
 |----------|------|--------------|
 | Sample0 | ENSG00000132549 | Splice_variant |
 | Sample0 | ENSG00000167548 | Missense |
 | Sample0 | ENSG00000182359 | Missense |
 | Sample0 | ENSG00000165474 | Stop_Gained |
 | Sample1 | ENSG00000164303 | Missense |

**Notes:**
- Each row represents a unique gene-variant combination for a specific sample
- Multiple rows with the same SampleID are allowed (and expected) as each individual typically has multiple variants
- The `Variant_Type` column should use standardized variant nomenclature

### Phenotyping DataFrame

The phenotyping dataframe contains phenotypic information and covariates for each individual. It must include the following:

**Mandatory Columns:**
 | Column | Type | Description | Example |
 |--------|------|-------------|---------|
 | `SampleID` | str | Unique identifier for each sample/individual (must match genotyping dataframe) | "Sample0", "Sample1" |
 | `[Phenotype]` | numeric | Phenotypic measurement (column name can be customized) | 0, 1, 0.45, 120.5 |

**Optional Columns:**
 | Column | Type | Description | Example |
 |--------|------|-------------|---------|
 | `FID` | str/int | Family identifier (useful for family-based studies) | 16, 21 |
 | `[Covariate1]` | numeric | Continuous or categorical covariates | -0.0023, 25.4 |
 | `[Covariate2]` | numeric | Additional covariates | -0.0189, 1.2 |

**Example:**
 | SampleID | FID | Phenotype | Covariate1 | Covariate2 |
 |----------|-----|-----------|------------|------------|
 | Sample0 | 16 | 0 | -0.0023 | -0.0189 |
 | Sample1 | 16 | 1 | -0.0021 | -0.0204 |
 | Sample2 | 21 | 0 | -0.0010 | -0.0200 |

**Requirements:**
1. The `SampleID` column must exist in both dataframes and contain matching values
2. The phenotype column must contain numerical data (binary or continuous)
3. Any additional columns can be used as covariates or interaction terms if they contain numerical data
4. For family-based studies, the `FID` column can be used to specify family relationships

**Data Notes:**
- Missing values should be handled appropriately before analysis
- Covariates should be checked for collinearity
- Binary phenotypes should be coded as 0/1
- Continuous phenotypes should be normally distributed or transformed if necessary
- Sample sizes should be sufficient for the planned statistical analysis


### Single Gene Set Input

For functions analyzing individual gene sets (`FunBurd_one_gene_set` and `FunBurd_one_gene_set_wih_name`), the gene set should be provided as:

#### Format
- A simple Python list of gene identifiers
- Gene names must match exactly with those in the `Gene` column of the genotyping dataframe

#### Example
```python
gene_list = [
    "ENSG00000132549",
    "ENSG00000167548",
    "ENSG00000182359",
    "ENSG00000165474"
]
```

### Multiple Gene Sets Input

For functions analyzing multiple gene sets (FunBurd_multiple_gene_sets), the input should be structured as:

#### Format

A Python dictionary where:
Keys are descriptive names for each gene set
Values are lists of gene identifiers
Gene names must match exactly with those in the Gene column of the genotyping dataframe

#### Example
```python
gene_sets_dict = {
    "cel_type 1": [
        "ENSG00000164303",
        "ENSG00000010626",
        "ENSG00000112559"
    ],
    "cell_type 2": [
        "ENSG00000139144",
        "ENSG00000131016",
        "ENSG00000157800"
    ],
    "cell type 3": [
        "ENSG00000134313",
        "ENSG00000102763",
        "ENSG00000129467"
    ]
}
```

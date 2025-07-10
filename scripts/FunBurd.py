import pandas as pd

import scripts.Burden_fcts as Burden_fcts
import scripts.data_preparation as data_preparation
from joblib import Parallel, delayed
from tqdm import tqdm
from multiprocessing import Pool
import time


def FunBurd_one_gene_set(phenotype_table: pd.DataFrame,
                         phenotype_column_name: str,
                         variants_table: pd.DataFrame,
                         gene_list: list,
                         keep_all_regression_results: bool = False,
                         list_covariantes: list = None,
                         interaction_column: list = None,
                         column_for_conditional: str = None,
                         correction_outside_gene_set: dict = None):
    """
    Perform a burden analysis for a single gene set.
    Args:
        phenotype_table (pd.DataFrame): DataFrame containing the phenotype data.
        phenotype_column_name (str): Name of the column in the phenotype table that contains the phenotype scores (Can be either continuous or binary, the model will adapt accordingly).
        variants_table (pd.DataFrame): DataFrame containing the variants altering each gene for each individual of the analysis.
        gene_list (list): List of genes to collapse as a unique functional group and used for the burden analysis. The names of the genes have to match the 'Gene' column in the variants table.
        keep_all_regression_results (bool): If False, will only keep the results of the regression model for the gene burden. If True, will keep all the results of the regression model (including covariates and interactions). Default is False.
        list_covariantes (list): List of covariates to include in the regression model in addition of the gene burden (Can be ancestry, sex, age, etc.). If not None the column names have to be in the phenotype table. Default is None.
        interaction_column (list): List of columns for interaction terms in the regression model. If not None the column name has to be in the phenotype table. Default is None.
        column_for_conditional (str): Column name for conditional analysis. If not None, the column name has to be in the phenotype table. Default is None. 
        correction_outside_gene_set (dict): Dictionary containing the categories of the genes in the genome. Could be LOEUF categories for example (LOEUF below or Above 1), if only one category included in the dictionnary, will corrrect for the average of the burden analysis outside the gene set. If None, no correction will be applied. This correction is mostly usefull for multigenic variants such as CNVs. Default is None.
    Returns:
        pd.DataFrame: DataFrame containing the results of the burden analysis with each line representing a covariate association with the phenotype.
    """

    # Check if all the mandatory columns are present in the phenotype and variants tables
    data_preparation.check_mandatory_columns(
        phenotype_table, variants_table, phenotype_column_name, list_covariantes, interaction_column, column_for_conditional)

    # Get the type of phenotype distribution to adapt the model (logisitic regression for binary, linear regression for continuous)
    pheno_type = data_preparation.get_phenotype_type(
        phenotype_table, phenotype_column_name)

    options = {
        'pheno_score': phenotype_column_name,
        'pheno_type': pheno_type,
        'correction_outside_gene_set': correction_outside_gene_set,
    }

    # Collapse the variants table to get the burden of the gene set
    annotated_collapsed_genes = data_preparation.collapse_genes_and_phenotypes(
        variants_table, gene_list, phenotype_table)

    if correction_outside_gene_set is not None:
        # If correction outside gene set is provided, we will add the column to the annotated_collapsed_genes
        annotated_collapsed_genes = data_preparation.add_outside_gene_set_correction(
            annotated_collapsed_genes, variants_table, gene_list, correction_outside_gene_set)

    # Compute the burden model
    model_result, log_model = Burden_fcts.prepare_and_run_model(
        annotated_collapsed_genes, options, keep_all_regression_results, list_covariantes, interaction_column, column_for_conditional)
    print(log_model)
    return model_result


def FunBurd_one_gene_set_wih_name(phenotype_table: pd.DataFrame,
                                  phenotype_column_name: str,
                                  variants_table: pd.DataFrame,
                                  genes_item: list,
                                  keep_all_regression_results: bool = False,
                                  list_covariantes: list = None,
                                  interaction_column: list = None,
                                  column_for_conditional: str = None,
                                  correction_outside_gene_set: dict = None):

    result = FunBurd_one_gene_set(
        phenotype_table=phenotype_table,
        phenotype_column_name=phenotype_column_name,
        variants_table=variants_table,
        gene_list=genes_item[1],
        keep_all_regression_results=keep_all_regression_results,
        list_covariantes=list_covariantes,
        interaction_column=interaction_column,
        column_for_conditional=column_for_conditional,
        correction_outside_gene_set=correction_outside_gene_set
    )
    result['Gene_set'] = genes_item[0]
    result['Variable'] = result.index
    return result


def FunBurd_multiple_gene_sets(phenotype_table: pd.DataFrame,
                               phenotype_column_name: str,
                               variants_table: pd.DataFrame,
                               gene_sets_dict: dict,
                               keep_all_regression_results: bool = False,
                               list_covariantes: list = None,
                               interaction_column: list = None,
                               column_for_conditional: str = None,
                               correction_outside_gene_set: dict = None,
                               n_cpus: int = 1):

    # Load the gene sets from the file

    # Check if the number of CPUs is greater than 1 for parallel processing
    if n_cpus < 1:
        raise ValueError("Number of CPUs must be at least 1.")
    if n_cpus == 1:
        print("Running in single-threaded mode. For parallel processing, set n_cpus > 1.")

    if n_cpus > 1:
        gene_sets_items = list(gene_sets_dict.items())

        results = Parallel(n_jobs=n_cpus, batch_size=4)(delayed(FunBurd_one_gene_set_wih_name)(
            phenotype_table=phenotype_table,
            phenotype_column_name=phenotype_column_name,
            variants_table=variants_table,
            genes_item=gene_set_item,
            keep_all_regression_results=keep_all_regression_results,
            list_covariantes=list_covariantes,
            interaction_column=interaction_column,
            column_for_conditional=column_for_conditional,
            correction_outside_gene_set=correction_outside_gene_set) for gene_set_item in tqdm(gene_sets_items, desc="Processing gene sets"))
    else:
        results = []
        for gene_set_item in gene_sets_dict.items():
            print(f"Processing gene set: {gene_set_item[0]}")

            result = FunBurd_one_gene_set_wih_name(
                phenotype_table=phenotype_table,
                phenotype_column_name=phenotype_column_name,
                variants_table=variants_table,
                genes_item=gene_set_item,
                keep_all_regression_results=keep_all_regression_results,
                list_covariantes=list_covariantes,
                interaction_column=interaction_column,
                column_for_conditional=column_for_conditional,
                correction_outside_gene_set=correction_outside_gene_set
            )
            results.append(result)

    results_df = pd.concat(results, ignore_index=True)
    return results_df


def bootstraped_FunBurd(phenotype_table: pd.DataFrame,
                        phenotype_column_name: str,
                        variants_table: pd.DataFrame,
                        gene_list: list,
                        n_bootstrap: int = 1000,
                        n_cpus: int = 1,
                        seed: int = 42,
                        bootstrap_size: int = None,
                        keep_all_regression_results: bool = False,
                        list_covariantes: list = None,
                        interaction_column: list = None,
                        column_for_conditional: str = None,
                        correction_outside_gene_set: dict = None):

    if bootstrap_size is None:
        bootstrap_size = len(gene_list) / 2
    else:
        if bootstrap_size > len(gene_list):
            raise ValueError(
                "Bootstrap size cannot be larger than the number of genes in the gene list.")

    # Create the bootstrap samples
    bootstrap_samples = data_preparation.create_bootstrap_samples(
        gene_list=gene_list,
        n_bootstrap=n_bootstrap,
        bootstrap_size=bootstrap_size,
        seed=seed)

    result = FunBurd_one_gene_set(
        phenotype_table=phenotype_table,
        phenotype_column_name=phenotype_column_name,
        variants_table=variants_table,
        gene_list=gene_list,
        keep_all_regression_results=keep_all_regression_results,
        list_covariantes=list_covariantes,
        interaction_column=interaction_column,
        column_for_conditional=column_for_conditional
    )

    return 'TODO'

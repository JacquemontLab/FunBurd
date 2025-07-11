import pandas as pd

import scripts.Burden_fcts as Burden_fcts
import scripts.data_preparation as data_preparation
import scripts.bootstap_functions as bootstrap_functions

from joblib import Parallel, delayed
from tqdm import tqdm
import statsmodels.stats.multitest as ssm


def FunBurd_one_gene_set(phenotype_table: pd.DataFrame,
                         phenotype_column_name: str,
                         variants_table: pd.DataFrame,
                         gene_list: list,
                         keep_all_regression_results: bool = False,
                         list_covariantes: list = None,
                         interaction_column: list = None,
                         column_for_conditional: str = None,
                         correction_outside_gene_set: dict = None) -> pd.DataFrame:
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
    if annotated_collapsed_genes.empty:
        print(
            "No variants found for the provided gene list. Please check the gene names and the variants table.")
        return pd.DataFrame()

    if correction_outside_gene_set is not None:
        # If correction outside gene set is provided, we will add the column to the annotated_collapsed_genes
        annotated_collapsed_genes = data_preparation.add_outside_gene_set_correction(
            annotated_collapsed_genes, variants_table, gene_list, correction_outside_gene_set)

    # Compute the burden model
    model_result, log_model = Burden_fcts.prepare_and_run_model(
        annotated_collapsed_genes, options, keep_all_regression_results, list_covariantes, interaction_column, column_for_conditional)
    print(log_model)
    print('-' * 50)
    print('')
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
    """
    Perform a burden analysis for a single gene set with its name.
    Args:
        phenotype_table (pd.DataFrame): DataFrame containing the phenotype data.
        phenotype_column_name (str): Name of the column in the phenotype table that contains the phenotype scores (Can be either continuous or binary, the model will adapt accordingly).
        variants_table (pd.DataFrame): DataFrame containing the variants altering each gene for each individual of the analysis.
        genes_item (list): List containing the name of the gene set and the list of genes to collapse as a unique functional group for the burden analysis. The first element is the name of the gene set and the second element is the list of genes. The names of the genes have to match the 'Gene' column in the variants table.
        keep_all_regression_results (bool): If False, will only keep the results of the regression model for the gene burden. If True, will keep all the results of the regression model (including covariates and interactions). Default is False.
        list_covariantes (list): List of covariates to include in the regression model in addition of the gene burden (Can be ancestry, sex, age, etc.). If not None the column names have to be in the phenotype table. Default is None.
        interaction_column (list): List of columns for interaction terms in the regression model. If not None the column name has to be in the phenotype table. Default is None.
        column_for_conditional (str): Column name for conditional analysis. If not None, the column name has to be in the phenotype table. Default is None.
        correction_outside_gene_set (dict): Dictionary containing the categories of the genes in the genome. Could be LOEUF categories for example (LOEUF below or Above 1), if only one category included in the dictionnary, will corrrect for the average of the burden analysis outside the gene set. If None, no correction will be applied. This correction is mostly useful for multigenic variants such as CNVs. Default is None.
    Returns:
        pd.DataFrame: DataFrame containing the results of the burden analysis with each line representing a covariate association with the phenotype for the gene set.
    """
    
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
    if result is not None and not result.empty:
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
                               n_cpus: int = 1) -> pd.DataFrame:
    """
    Perform a burden analysis for multiple gene sets in parallel.
    Args:
        phenotype_table (pd.DataFrame): DataFrame containing the phenotype data.
        phenotype_column_name (str): Name of the column in the phenotype table that contains the phenotype scores (Can be either continuous or binary, the model will adapt accordingly).
        variants_table (pd.DataFrame): DataFrame containing the variants altering each gene for each individual of the analysis.
        gene_sets_dict (dict): Dictionary where keys are gene set names and values are lists of genes to collapse as unique functional groups for the burden analysis. The names of the genes have to match the 'Gene' column in the variants table.
        keep_all_regression_results (bool): If False, will only keep the results of the regression model for the gene burden. If True, will keep all the results of the regression model (including covariates and interactions). Default is False.
        list_covariantes (list): List of covariates to include in the regression model in addition of the gene burden (Can be ancestry, sex, age, etc.). If not None the column names have to be in the phenotype table. Default is None.
        interaction_column (list): List of columns for interaction terms in the regression model. If not None the column name has to be in the phenotype table. Default is None.
        column_for_conditional (str): Column name for conditional analysis. If not None, the column name has to be in the phenotype table. Default is None.
        correction_outside_gene_set (dict): Dictionary containing the categories of the genes in the genome. Could be LOEUF categories for example (LOEUF below or Above 1), if only one category included in the dictionnary, will corrrect for the average of the burden analysis outside the gene set. If None, no correction will be applied. This correction is mostly useful for multigenic variants such as CNVs. Default is None.
        n_cpus (int): Number of CPUs to use for parallel processing. Default is 1 (single-threaded).
    Returns:
        pd.DataFrame: DataFrame containing the results of the burden analysis with each line representing a covariate association with the phenotype for each gene set.
    """

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
                        bootstrap_size: int = None,
                        pvalue_threshold: float = 0.05,
                        return_effect_size_df: bool = False,
                        n_cpus: int = 1,
                        seed: int = None,
                        list_covariantes: list = None,
                        interaction_column: list = None,
                        column_for_conditional: str = None,
                        correction_outside_gene_set: dict = None) -> dict:
    """
    Perform a bootstrap burden analysis for a list of genes.
    Args:
        phenotype_table (pd.DataFrame): DataFrame containing the phenotype data.
        phenotype_column_name (str): Name of the column in the phenotype table that contains the phenotype scores (Can be either continuous or binary, the model will adapt accordingly).
        variants_table (pd.DataFrame): DataFrame containing the variants altering each gene for each individual of the analysis.
        gene_list (list): List of genes to collapse as a unique functional group and used for the burden analysis. The names of the genes have to match the 'Gene' column in the variants table.
        n_bootstrap (int): Number of bootstrap samples to create. Default is 1000.
        bootstrap_size (int): Size of each bootstrap sample. If None, will be set to half the size of the gene list. If larger than the gene list, will raise an error.
        pvalue_threshold (float): P-value threshold for significance. Default is 0.05.
        return_effect_size_df (bool): If True, will return the DataFrame with the effect sizes of the gene sets. Default is False.
        n_cpus (int): Number of CPUs to use for parallel processing. Default is 1 (single-threaded).
        seed (int): Random seed for reproducibility. Default is None.
        list_covariantes (list): List of covariates to include in the regression model in addition of the gene burden (Can be ancestry, sex, age, etc.). If not None the column names have to be in the phenotype table. Default is None.
        interaction_column (list): List of columns for interaction terms in the regression model. If not None the column name has to be in the phenotype table. Default is None.
        column_for_conditional (str): Column name for conditional analysis. If not None, the column name has to be in the phenotype table. Default is None.
        correction_outside_gene_set (dict): Dictionary containing the categories of the genes in the genome. Could be LOEUF categories for example (LOEUF below or Above 1), if only one category included in the dictionnary, will corrrect for the average of the burden analysis outside the gene set. If None, no correction will be applied. This correction is mostly useful for multigenic variants such as CNVs. Default is None.
    Returns:
        dict: A dictionary containing the results of the bootstrap burden analysis. If return_effect_size_df is True, the dictionary will also contain a DataFrame with the effect sizes of the gene sets.
    """

    if bootstrap_size is None:
        bootstrap_size = len(gene_list) / 2
    else:
        if bootstrap_size > len(gene_list):
            raise ValueError(
                "Bootstrap size cannot be larger than the number of genes in the gene list.")

    # Create the bootstrap samples
    bootstrap_samples = bootstrap_functions.create_bootstrap_samples(
        gene_list=gene_list,
        n_bootstrap=n_bootstrap,
        bootstrap_size=bootstrap_size,
        seed=seed)

    # Convert the bootstrap samples to a DataFrame
    # Each row corresponds to a bootstrap sample and contains the genes in that sample
    # The index of the DataFrame corresponds to the bootstrap sample identifier
    bootstrap_samples_df = pd.DataFrame({'bootstrap_idx': bootstrap_samples.keys(),
                                         'genes': bootstrap_samples.values()})

    results = FunBurd_multiple_gene_sets(
        phenotype_table=phenotype_table,
        phenotype_column_name=phenotype_column_name,
        variants_table=variants_table,
        gene_sets_dict=bootstrap_samples,
        keep_all_regression_results=False,
        list_covariantes=list_covariantes,
        interaction_column=interaction_column,
        column_for_conditional=column_for_conditional,
        correction_outside_gene_set=correction_outside_gene_set,
        n_cpus=n_cpus)

    results = results[results['Variable'] == 'GeneSet']
    results['pvalue_adj'] = ssm.fdrcorrection(results['pvalue'])[1]

    # keep only sig results
    sig_results = results[results['pvalue_adj'] < pvalue_threshold]
    print(f"Number of significant subsets: {len(sig_results)}")
    count_in_sig = bootstrap_functions.get_count_genes_in_subset(
        sig_results, bootstrap_samples_df)
    
    non_sig_results = results[~results['Gene_set'].isin(sig_results['Gene_set'])]

    count_in_non_sig = bootstrap_functions.get_count_genes_in_subset(
        non_sig_results, bootstrap_samples_df)
    res_enrich = bootstrap_functions.compute_enrichment_in_sig_subset(count_in_sig, count_in_non_sig)

    dict_res = {'Gene enrichment': res_enrich}

    if return_effect_size_df:
        dict_res['Subset effect size'] = results
    
    return dict_res
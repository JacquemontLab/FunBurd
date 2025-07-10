import random
import numpy as np
import pandas as pd
import itertools
from scipy import stats


def create_bootstrap_samples(gene_list: list, n_bootstrap: int, bootstrap_size: int, seed: int) -> dict:
    """
    Create bootstrap samples from the gene list.

    :param gene_list: List of genes to sample from.
    :param n_bootstrap: Number of bootstrap samples to create.
    :param bootstrap_size: Size of each bootstrap sample.
    :param seed: Random seed for reproducibility.
    :return: A dictionary where keys are bootstrap sample identifiers and values are lists of genes in each sample.
    """
    if seed is not None:
        random.seed(seed)

    bootstrap_samples = []
    for _ in range(n_bootstrap):
        sample = np.random.choice(gene_list, size=bootstrap_size, replace=False)
        bootstrap_samples.append(sample)

    bootstrap_sample_dict = {
        f"Bootstrap_{i}": sample for i, sample in enumerate(bootstrap_samples)}
    return bootstrap_sample_dict


def compute_Fisher(count_sig, count_non_sig, rest_sig, rest_non_sig, alternative):
    # Create a 2x2 contingency table
    contingency_table = np.array([[count_sig, count_non_sig],
                                  [rest_sig, rest_non_sig]])
    # Perform Fisher's exact test
    odds_ratio, p_value = stats.fisher_exact(
        contingency_table, alternative=alternative)
    return odds_ratio, p_value


def get_count_genes_in_subset(subset_es: pd.DataFrame, bootstrap_samples_df: pd.DataFrame) -> pd.DataFrame:

    subset_es = pd.merge(
        subset_es, bootstrap_samples_df, left_on='Gene_set', right_on='bootstrap_idx', how='left')
    full_list_genes_in_subset = list(
        itertools.chain.from_iterable(subset_es.genes.values.tolist()))
    full_list_genes_in_subset = pd.Series(full_list_genes_in_subset)

    # Count the number of occurrences of each gene in the subset
    count_significant_genes = full_list_genes_in_subset.value_counts().reset_index()
    count_significant_genes.columns = ['Gene', 'Count']
    count_significant_genes['rest_count'] = subset_es.shape[0] - count_significant_genes['Count']
    return count_significant_genes


# non_significant_es = results[(results['pvalue_adj'] >= 0.05)]
# non_significant_es = pd.merge(
#     non_significant_es, bootstrap_samples_df, left_on='Gene_set', right_index=True, how='left')

# sig_genes = sig_gene_sets.genes.str.split(';', expand=True).stack(
# ).reset_index(level=1, drop=True).rename('gene_id')

# non_sig_genes = non_sig_gene_sets.genes.str.split(
#     ';', expand=True).stack().reset_index(level=1, drop=True).rename('gene_id')

# # Count the number of occurrences of each gene in the significant and not significant gene-sets
# sig_genes = sig_genes.value_counts().reset_index()
# sig_genes.columns = ['gene_id', 'count']
# sig_genes['rest_count'] = sig_gene_sets.shape[0] - sig_genes['count']
# non_sig_genes = non_sig_genes.value_counts().reset_index()
# non_sig_genes.columns = ['gene_id', 'count']
# non_sig_genes['rest_count'] = non_sig_gene_sets.shape[0] - \
#     non_sig_genes['count']

def compute_enrichment_in_sig_subset(count_genes_in_significant: pd.DataFrame,
                                     count_genes_in_non_significant: pd.DataFrame) -> pd.DataFrame:
    # Merge the counts of the significant and non-significant gene-sets
    merged_counts = pd.merge(
        count_genes_in_significant, count_genes_in_non_significant, on='Gene', how='outer', suffixes=('_sig', '_non_sig'))
    merged_counts.fillna(0, inplace=True)

    merged_counts['N_subsets'] = merged_counts['Count_sig'] + merged_counts['Count_non_sig']

    # Compute the Fisher's exact test for each gene
    merged_counts['Enrichment_OR'], merged_counts['pvalue_enrichment'] = zip(*merged_counts.apply(lambda x: compute_Fisher(
        x['Count_sig'], x['Count_non_sig'], x['rest_count_sig'], x['rest_count_non_sig'], 'greater'), axis=1))
    
    return merged_counts[['Gene', 'N_subsets', 'Count_sig', 'Count_non_sig', 'Enrichment_OR', 'pvalue_enrichment']]

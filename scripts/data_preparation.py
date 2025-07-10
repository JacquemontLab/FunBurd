import pandas as pd


def get_phenotype_type(phenotype_df, pheno_column_name):
    # Get the phenotype type
    number_of_unique_values = phenotype_df[pheno_column_name].nunique()
    if number_of_unique_values == 2:
        pheno_type = 'binary'
    elif number_of_unique_values > 2:
        pheno_type = 'continuous'
    else:
        raise ValueError('Phenotype score not binary nor continuous..')
    print(f'The phenotype column {pheno_column_name} has {number_of_unique_values} unique values, '
          f'therefore it is considered as a {pheno_type} phenotype.')
    return pheno_type


def check_mandatory_columns(phenotype_table, variants_table, phenotype_column_name, list_covariantes, interaction_term, family_id_column):
    # Check the mandatory columns in the different tables
    if 'SampleID' not in phenotype_table.columns:
        raise ValueError("Phenotype table must contain 'SampleID' column.")
    if 'SampleID' not in variants_table.columns:
        raise ValueError("Variants table must contain 'SampleID' column.")
    if 'Gene' not in variants_table.columns:
        raise ValueError("Variants table must contain 'Gene' column.")
    if phenotype_column_name not in phenotype_table.columns:
        raise ValueError(
            f"Phenotype table must contain '{phenotype_column_name}' column.")
    if list_covariantes is not None:
        for covariate in list_covariantes:
            if covariate not in phenotype_table.columns:
                raise ValueError(
                    f"Phenotype table must contain the '{covariate}' column.")
    if interaction_term is not None and interaction_term not in phenotype_table.columns:
        raise ValueError(
            f"Phenotype table must contain the '{interaction_term}' column for interaction term.")
    if family_id_column is not None and family_id_column not in phenotype_table.columns:
        raise ValueError(
            f"Phenotype table must contain the '{family_id_column}' column for the conditional model to be run.")


def collapse_genes_and_phenotypes(variants_table: pd.DataFrame, gene_list: list, phenotype_table: pd.DataFrame):
    # Collapse the number of variants per gene per individual
    variants_in_gene_list = variants_table[variants_table['Gene'].isin(
        gene_list)]

    variants_in_gene_list = variants_in_gene_list[['SampleID', 'Gene']]
    collpased_genes = variants_in_gene_list.groupby(
        ['SampleID']).size().reset_index(name='GeneSet')

    annotated_collapsed_genes = pd.merge(
        phenotype_table, collpased_genes, on='SampleID', how='left')
    annotated_collapsed_genes = annotated_collapsed_genes.fillna(
        0)  # Fill NaN values with 0 for individuals without variants
    return annotated_collapsed_genes


def add_outside_gene_set_correction(annotated_collapsed_genes, variants_table, gene_list, correction_dictionary):
    # Add the correction for the outside gene set
    variants_outside_gene_set = variants_table[~variants_table['Gene'].isin(
        gene_list)]

    list_cat_outisde = list(set(correction_dictionary.values()))
    for cat in list_cat_outisde:
        list_genes_in_cat = [
            gene for gene, category in correction_dictionary.items() if category == cat]

        variants_in_cat = variants_outside_gene_set[variants_outside_gene_set['Gene'].isin(
            list_genes_in_cat)]
        if not variants_in_cat.empty:
            counts_per_individual = variants_in_cat.groupby(
                'SampleID').size().reset_index(name=f'OutsideGeneSet_{cat}')
            annotated_collapsed_genes = pd.merge(
                annotated_collapsed_genes, counts_per_individual, on='SampleID', how='left')
            annotated_collapsed_genes[f'OutsideGeneSet_{cat}'] = annotated_collapsed_genes[f'OutsideGeneSet_{cat}'].fillna(
                0)
        else:
            print(
                f"No variants found outside of the gene set for category '{cat}'.")
    return annotated_collapsed_genes

import numpy as np
import pandas as pd
import time
import statsmodels.formula.api as sm
from statsmodels.discrete.conditional_models import ConditionalLogit

import warnings
from statsmodels.tools.sm_exceptions import PerfectSeparationWarning, ConvergenceWarning
warnings.simplefilter('ignore', PerfectSeparationWarning)
warnings.simplefilter('ignore', RuntimeWarning)


def print_model_output_logs(print_output_model, variant_type, pheno_score, gene_set_name, model_type):
    ''' Print the model output logs
    Args:
        print_output_model (str): The output of the model.
    '''
    print('--------------')
    print(
        f'Model run for {variant_type} in {gene_set_name} for {pheno_score} with {model_type} model')
    print(print_output_model)
    print('--------------', flush=True)
    return None


def add_interaction_term_and_covariates(formula, other_covariates, interaction_term):
    ''' Add interaction terms to the first element of the formula
    Args:
        formula (str): The formula for the regression.
        interaction_term (str): The interaction term to add to the formula.
    Returns:
        str: The formula with the interaction term added.
    '''
    pheno = formula.split('~')[0].strip()
    main_covariates = formula.split('~')[1].strip()
    main_covariates = main_covariates.split(' + ')

    if interaction_term is not None:
        main_covariates.append(f'{main_covariates[0]} * {interaction_term}')

    if other_covariates is not None:
        formula_new = f'{pheno} ~ ' + \
            ' + '.join(main_covariates) + ' + ' + ' + '.join(other_covariates)
        print(f'Formula: {formula_new}')
    else:
        formula_new = formula
    return formula_new


def output_model_summary(reg, covariates, keep_all_regression_results):
    indexes = reg.params.index.tolist()

    if not keep_all_regression_results:
        # Keep only the main covariates and the interaction term
        index_val_of_cov = [indexes.index(idx)
                            for idx in indexes if idx not in covariates]
        reg.params = reg.params.iloc[index_val_of_cov]
        reg.bse = reg.bse.iloc[index_val_of_cov]
        reg.pvalues = reg.pvalues.iloc[index_val_of_cov]
        indexes = [indexes[i] for i in index_val_of_cov]

    # Create a DataFrame with the results
    final_df = pd.DataFrame({'Estimate': reg.params,
                             'SE': reg.bse,
                             'pvalue': reg.pvalues,
                             }, index=indexes
                            )
    return final_df


def prepare_and_run_model(
        data: pd.DataFrame,
        options: dict,
        keep_all_regression_results: bool,
        covariates: list,
        interaction_term: str,
        conditional_column: str):

    pheno_score = options['pheno_score']
    pheno_type = options['pheno_type']
    if options['correction_outside_gene_set'] is not None:
        correction_outside_gene_set = options['correction_outside_gene_set']

    list_main_covariates = [
        col for col in data.columns if 'GeneSet' in col]
    if len(list_main_covariates) == 0:
        raise ValueError(
            "No 'GeneSet' column found in the data. Please ensure that the data contains a column with the number of genes disrupted by category (gene-set and/or LOEUF cat) per individual.")

    formula = "{} ~ {}".format(pheno_score,
                               ' + '.join(list_main_covariates))

    if interaction_term is not None or covariates is not None:
        formula = add_interaction_term_and_covariates(
            formula, covariates, interaction_term)

    # Run the linear regression following formula and based on data. Capture error and print exception
    reg, print_output_model = run_regression_model(
        data, formula, pheno_type, conditional_column)

    # Return the results of the regression
    if reg is not None:
        final_df = output_model_summary(
            reg, covariates, keep_all_regression_results)

        return final_df, print_output_model
    else:
        return None, None


def run_logistic_regression_model(data, formula):
    try:
        with warnings.catch_warnings(record=True) as w:
            reg = sm.logit(formula, data=data).fit(maxiter=100, disp=0)
            # ignore any non-custom warnings that may be in the list
            w = filter(lambda i: issubclass(i.category, ConvergenceWarning), w)

            if any(w):
                # do something with the first warning
                print(
                    "Convergence warning encountered. This usually means that the model did not converge.")
                print("Please be careful with the following results.")

            return reg

    except Exception as e:
        if e.args[0] == 'Singular matrix':
            print("Error: Singular matrix encountered. This usually means that there is a perfect multicollinearity in the data.")
            print("Please check your data and remove any redundant variables.")
            return None
        else:
            print("Error: ", e)
            return None


def run_conditional_logit_model(data, formula):
    try:
        covariates = [var.strip() for var in formula.split("~")[1].split(
            "+")] + [formula.split("~")[0].strip()] + ['FID']
        interactions = [var for var in covariates if '*' in var]
        for i in interactions:
            data[i] = data.eval(i)
        X = data[covariates]
        X = X.dropna()
        y = X[formula.split("~")[0].strip()]
        groups = X['FID']

        X = X.drop(columns=[formula.split("~")[0].strip()] + ['FID'])

        X = pd.get_dummies(X, drop_first=True)

        model = ConditionalLogit(y, X, groups=groups)
        list_families = [list(x) for x in model._endog_grp]
        asd_status = pd.Series(np.concatenate(list_families).tolist())
        asd_counts = asd_status.value_counts()
        reg = model.fit(disp=0)
        print_output = f'Number of individuals {asd_status.shape[0]}: with ASD: {asd_counts[1]}, without ASD: {asd_counts[0]}\n'
        return reg, print_output

    except Exception as e:
        print('----------')
        print("Error: ", e)
        print(formula)
        print(
            f'Pheno distri: {data[formula.split(" ~ ")[0]].value_counts()}')
        print(
            f'Variable counts\n{data[formula.split(" ~ ")[1].split(" + ")].sum(axis=0)}')
        print('--------------\n')
        return None, ''


def run_regression_model(data, formula, pheno_type, conditional_column):
    ''' Run the regression model
    Args:
        data (pd.DataFrame): The dataframe containing the number of genes disrupted by category (gene-set and/or LOEUF cat) per individual + their phenotypes.
        formula (str): The formula for the regression.
        pheno_type (str): The type of phenotype.
    Returns:
        pd.DataFrame: The results of the regression.
    '''
    start = time.time()
    # Check if the formula is correct (remove the double spaces)
    formula = formula.replace('  ', ' ')

    if pheno_type == 'binary':
        if conditional_column is None:
            reg = run_logistic_regression_model(data, formula)
            print_output = ''
        else:
            reg, print_output = run_conditional_logit_model(data, formula)

    elif pheno_type == 'continuous':
        try:
            reg = sm.gls(formula, data=data).fit()
        except Exception as e:
            print("Error: ", e)
            return pd.DataFrame(), None
    else:
        raise ValueError(
            'The type of phenotype is unknown, it has to be binary or continuous')
    end = time.time()
    print_output = print_output + \
        f'Time taken for the model: {end - start:.2f} seconds'
    return reg, print_output

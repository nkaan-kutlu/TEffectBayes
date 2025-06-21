process BNM_CALCULATION {
    tag "bnm_calculation"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'nkaankutlu/pytools:latest'

    publishDir "${params.outdir}/BNM/calculations/", mode: 'copy'
    
    input:
    path samplesheet
    path samplesheet_chip
    path bnm_gene_input_pickles
    
    output:
    path "BNM/calculations/*.pkl", emit: gene_model_pickles
    
    script:
    """
    mkdir -p BNM/calculations

    python <<EOF    
    import os
    import pickle
    from itertools import permutations
    import pandas as pd
    import bnlearn as bn
    import numpy as np
    from sklearn.feature_selection import mutual_info_regression
    from sklearn.utils import resample
    from sklearn.preprocessing import KBinsDiscretizer
    import logging

    # Set up logging
    logging.basicConfig(level=logging.INFO)

    input_dir = ("${bnm_gene_input_pickles}")
    output_dir = "BNM/calculations"
    os.makedirs(output_dir, exist_ok=True)

    # Get biological condition mapping from samplesheet
    try:
        samplesheet_df = pd.read_csv("${samplesheet}")
        condition_order = samplesheet_df['condition'].dropna().drop_duplicates().tolist()
        condition_mapping = {cond: idx for idx, cond in enumerate(condition_order)}
        logging.info(f"Condition mapping: {condition_mapping}")
    except Exception as e:
        logging.error(f"Failed to read samplesheet or generate condition mapping: {e}")
        condition_mapping = {}
    
    # Read chip samplesheet
    chip_df = pd.read_csv("${samplesheet_chip}")

    # Get unique antibody names (case-insensitive matching for safety)
    antibody_names = chip_df['antibody'].dropna().unique().tolist()

    # Get all .pkl files in the input directory
    input_files = [f for f in os.listdir(input_dir) if f.endswith(".pkl")]

    for file_name in input_files:
        try:
            file_path = os.path.join(input_dir, file_name)
            logging.info(f"Processing file: {file_name}")
            
            # Load the pickle file
            with open(file_path, "rb") as f:
                df = pickle.load(f)
            
            # Apply condition mapping
            df['biological_condition'] = df['biological_condition'].astype(str)
            df['biological_condition'] = df['biological_condition'].map(condition_mapping)
            if df['biological_condition'].isnull().any():
                raise ValueError(f"Unmapped biological_condition values found in {file_name}")
            
            # Detect columns dynamically
            repeat_columns = [col for col in df.columns if '_TPM' in col]
            histone_columns = [col for col in df.columns if '_rpm' in col]
            gene_column = [col for col in df.columns if '_FPKM' in col]
            histone_repeat_columns = [col for col in df.columns if any(h.lower() in col.lower() for h in antibody_names)]
            specific_gene_column = gene_column[0]
            
            columns_of_interest = gene_column + histone_columns + repeat_columns + histone_repeat_columns + ['biological_condition']
            df = df[columns_of_interest]
            
            # Stratified bootstrapping
            classes = df['biological_condition'].unique()
            bootstrapped_df = pd.concat([
                resample(df[df['biological_condition'] == c], replace=True, n_samples=int(len(df) * 3 / len(classes)), random_state=42)
                for c in classes
            ])
            
            # Mutual information calculation
            mi = mutual_info_regression(bootstrapped_df.drop(columns=[specific_gene_column]), bootstrapped_df[specific_gene_column])
            
            # Define whitelist for edges
            white_list = list(set([("biological_condition", col) for col in gene_column + histone_columns + repeat_columns]))
            white_list += [(histone, specific_gene_column) for histone in histone_columns]
            white_list += [(repeat, specific_gene_column) for repeat in repeat_columns]
            
            # Generate blacklist
            nodes = bootstrapped_df.columns.tolist()
            all_possible_edges = list(permutations(nodes, 2))
            black_list = [edge for edge in all_possible_edges if edge not in white_list]
            fixed = white_list.copy()
            
            def add_conditional_edges(df, fixed, histone_columns, repeat_columns):
                for histone_repeat in histone_repeat_columns:
                    parts = histone_repeat.split('_', 1)
                    if len(parts) < 2:
                        continue
                    histone_name, repeat_name = parts
                    histone_col = next((col for col in histone_columns if histone_name.lower() in col.lower()), None)
                    repeat_col = f'{repeat_name}_TPM'                
                    if histone_col and repeat_col in df.columns:
                        fixed.append((histone_col, repeat_col))
                return fixed
            
            fixed = add_conditional_edges(bootstrapped_df, fixed, histone_columns, repeat_columns)
            
            # Learn Bayesian network structure
            model = bn.structure_learning.fit(bootstrapped_df, n_jobs=12, methodtype='hc', scoretype='bic',
                                            bw_list_method='List', white_list=white_list, black_list=black_list, fixed_edges=fixed)        
                       
            # Discretize data
            discretizer = KBinsDiscretizer(n_bins=2, encode='ordinal', strategy='uniform')
            discretized_data = discretizer.fit_transform(bootstrapped_df.drop(columns=['biological_condition']))
            discretized_df = pd.DataFrame(discretized_data, columns=bootstrapped_df.columns[:-1])
            discretized_df['biological_condition'] = bootstrapped_df['biological_condition'].values
            
            # Perform parameter learning using Maximum Likelihood Estimation
            model2 = bn.parameter_learning.fit(model,
                                            discretized_df,
                                            n_jobs=12,
                                            methodtype='maximumlikelihood',
                                            scoretype='k2')
                        
            # Learn Bayesian network parameters
            model_up = bn.parameter_learning.fit(model=model2, df=discretized_df, n_jobs=12, methodtype='bayes',
                                                scoretype='dirichlet', smooth=1, verbose=5)
                    
            # Save the model
            model_up_path = os.path.join(output_dir, f"{file_name.replace('bnm_input.pkl', '_calculation.pkl')}")
            bn.save(model_up, filepath=model_up_path, overwrite=True)
                                                            
            logging.info(f"Completed processing for {file_name}")
        
        except Exception as e:
            logging.error(f"Error processing {file_name}: {e}")


    EOF
    """
}
process BNM_INFERENCE {
    tag "bnm_inference"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'nkaankutlu/pytools:latest'

    publishDir "${params.outdir}/BNM/inferences/", mode: 'copy'
    
    input:
    path gene_model_pickles

    output:
    path "BNM/inferences/*.csv", emit: gene_model_inferences

    script:
    """
    mkdir -p BNM/inferences
    python <<EOF
    import pickle
    import pandas as pd
    import bnlearn as bn
    from itertools import product

    model_dir = ("${gene_model_pickles}")
    output_dir = "BNM/inferences"
    os.makedirs(output_dir, exist_ok=True)

    def run_inference_for_model(model_path, target_var_suffix='_FPKM'):
        with open(model_path, 'rb') as f:
            model = pickle.load(f)

        gene_name = os.path.basename(model_path).replace('_calculation.pkl', '')
        target_variable = f"{gene_name}{target_var_suffix}"

        model_variables = model['model'].states.keys() if isinstance(model, dict) else model.states.keys()
        model_variables = list(model_variables)

        evidence_variables = [v for v in model_variables if v != target_variable]

        cpd_lookup = model['model'].cpds if isinstance(model, dict) else model.cpds
        value_dict = {}

        for cpd in cpd_lookup:
            if cpd.variable in evidence_variables:
                value_dict[cpd.variable] = list(map(int, cpd.state_names[cpd.variable]))

        combinations = list(product(*value_dict.values()))
        evidence_keys = list(value_dict.keys())

        results = []
        for comb in combinations:
            evidence = dict(zip(evidence_keys, comb))
            inference_result = bn.inference.fit(model, variables=[target_variable], evidence=evidence)
            for _, row in inference_result.df.iterrows():
                result_row = {**evidence, target_variable: row[target_variable], 'Probability': row['p']}
                results.append(result_row)

        return pd.DataFrame(results)
    
    # Loop over all models
    for model_path in glob(os.path.join(model_dir, "*_calculation.pkl")):
        gene_name = os.path.basename(model_path).replace('_calculation.pkl', '')
        print(f"Running inference for {gene_name}")
        
        try:
            df = run_inference_for_model(model_path)
            output_path = os.path.join(output_dir, f"{gene_name}_inference_results.csv")
            df.to_csv(output_path, index=False)
        except Exception as e:
            print(f" Error processing {gene_name}: {e}")
    
    EOF
    """
}

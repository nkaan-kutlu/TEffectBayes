process BNM_VISUALIZATION {
    tag "bnm_visualization"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'nkaankutlu/pytools:latest'

    publishDir "${params.outdir}/BNM/visualizations/", mode: 'copy'
    
    input:
    path gene_model_pickles
    
    output:
    path "BNM/visualizations/*.pkl", emit: gene_model_graphs
    
    script:
    """
    mkdir -p BNM/visualizations

    python <<EOF
    import os
    import glob
    import pickle
    import bnlearn as bn
    import networkx as nx
    import matplotlib.pyplot as plt
    from networkx.drawing.nx_agraph import graphviz_layout

    model_dir = ("${gene_model_pickles}")
    output_dir = "BNM/visualizations"
    os.makedirs(output_dir, exist_ok=True)

    # Get all *_final.pkl files in the directory
    model_files = glob.glob(os.path.join(model_dir, "*_calculation.pkl"))

    # Function to rename nodes in the model
    def rename_nodes(model):
        new_adjmat = model['adjmat'].copy()
        new_names = {node: node.replace("_TPM", "").replace("_FPKM", "").replace("_rpm", "") for node in new_adjmat.columns}
        
        # Rename the adjacency matrix columns and index
        new_adjmat.rename(columns=new_names, index=new_names, inplace=True)
        
        # Update the model dictionary
        model['adjmat'] = new_adjmat
        return model

    # Process each model file
    for model_file in model_files:
        # Extract the base filename (without extension) for saving outputs
        base_name = os.path.basename(model_file).replace("_calculation.pkl", "")

        # Load the Bayesian model
        with open(model_file, 'rb') as f:
            model = pickle.load(f)
            
        node_properties = bn.get_node_properties(model)
        
        # Extract edges from the Bayesian model
        edges = model["model_edges"]  # Ensure your model contains the correct edges
        
        # Create a mapping to rename edges
        rename_mapping = {node: node.replace("_TPM", "").replace("_FPKM", "").replace("_rpm", "") for node in model['adjmat'].columns}
        
        # Rename edges
        edges = [(rename_mapping[src], rename_mapping[dst]) for src, dst in edges]
        
        # Rename nodes
        model = rename_nodes(model)
        
        # Create a NetworkX directed graph
        G = nx.DiGraph()
        G.add_edges_from(edges)
            
        # Use Graphviz layout for better edge positioning
        pos = graphviz_layout(G, prog="dot")  # 'dot' for hierarchical layout

        # Draw the graph
        plt.figure(figsize=(10, 8))
        nx.draw(G, pos, with_labels=True, node_color="lightblue", edge_color="gray", node_size=2000, font_size=10,
                connectionstyle="arc3,rad=0.1")


        # Save the figure as a PDF
        output_file = os.path.join(output_dir, f"{base_name}_model.pdf")
        plt.savefig(output_file, format="pdf", dpi=300)
        plt.close()  # Close the figure to free memory

        print(f"Saved {output_file}")

    print("All models processed successfully!")
    EOF
    """
}
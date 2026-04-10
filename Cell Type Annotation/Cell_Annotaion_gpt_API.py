from openai import OpenAI
import os
import json
import anndata
import scanpy as sc

import warnings
warnings.filterwarnings("ignore")

slide = "0050586"

#0050585
# NL = [1,2,13,14,25,26,37,38,49,50,61,62,73,74,85,86,97,98,109,110] 
# TC = [3,4,15,16,27,28,39,40,51,52,63,64,75,76,87,88,99,100,111,112]
# IM = [5,6,17,18,29,30,41,42,53,54,65,66,77,78,89,90,101,102,113,114]


#0050586
NL = [7,8,19,20,31,32,43,44,55,56,67,68,79,80,91,92,103,104,115,116]
TC = [9,10,21,22,33,34,45,46,57,58,69,70,81,82,93,94,105,106,117,118]   
IM = [11,12,23,24,35,36,47,48,59,60,71,72,83,84,95,96,107,108,119,120]



for core in TC:


    file_path = f"/home/longinus723/home/Anndata_Cluster/{slide}/clustered_core_{core}.h5ad"

    if not os.path.exists(file_path):
        print(f"core {core} 파일 없음⚠️")
        continue


    print("=========================================")
    print(f"🔄 core {core} Gene 추출 시작")
    adata = anndata.read_h5ad(file_path)

    sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon")
    ranked_genes = adata.uns["rank_genes_groups"]
    clusters = ranked_genes["names"].dtype.names  # 클러스터 이름 목록
    
    save_dir = f"/home/longinus723/home/Cell_Annotation/json/{slide}/core_{core}"
    os.makedirs(save_dir, exist_ok=True)
    
    
    
    for cluster in clusters:
        marker_genes = ranked_genes["names"][cluster][:15]
        print(f"  Cluster {cluster}: {marker_genes}")

        client = OpenAI(api_key="")

        dataset_type = "single-cell RNA-seq"
        sample_type = "CRC patient colon Tissue"
        # marker_genes = ["STMN1", "TYMS", "PCLAF", "CDT1", "TK1", "RRM2", "DHFR", "CENPM", "MKI67", "MCM4", "GINS2", "MCM2", "ZWINT", "CLSPN", "ASF1B", "CENPU", "E2F1", "UHRF1", "MYBL2", "NUSAP1", "HIST1H1B", "CDCA7", "FANCI", "TCF19", "CENPF"]

        # Define JSON Schema for structured response
        json_schema = {
            "type": "object",
            "properties": {
                "cell_type": {"type": "string", "description": "Predicted cell type based on marker genes"},
                "confidence_score": {"type": "number", "description": "Confidence score of the prediction"},
                "marker_gene_details": {
                    "type": "array",
                    "items": {
                        "type": "object",
                        "properties": {
                            "gene": {"type": "string", "description": "Marker gene name"},
                            "expression_pattern": {"type": "string", "description": "Expression pattern in the identified cell type"}
                        },
                        "required": ["gene", "expression_pattern"]
                    }
                }
            },
            "required": ["cell_type", "confidence_score", "marker_gene_details"]
        }

        # Define messages for GPT API request
        """messages = [
            {"role": "system", "content": "You are a helpful assistant for identifying cell types based on marker genes. "
                                        "Your response must be formatted as a JSON object."},
            {"role": "user", "content": f"I performed a DEGs analysis and clustering on my {dataset_type} data "
                                        f"from {sample_type} samples. "
                                        f"I found the following marker genes that are highly expressed in a specific cluster: "
                                        f"{marker_genes}. "
                                        f"Based on these marker genes, how do these marker genes specifically indicate the cell type of this cluster?"}
        ]"""
        ### Modify
        messages = [
            {"role": "system", "content": "You are a helpful assistant that identifies cell types based on marker genes. "
                                        "Your response MUST be a JSON object following this structure:\n\n"
                                        "{\n"
                                        '  "cell_type": "Predicted cell type",\n'
                                        '  "confidence_score": 0.95,\n'
                                        '  "marker_gene_details": [\n'
                                        '    {"gene": "PMEL", "expression_pattern": "Highly expressed in melanocytes"},\n'
                                        '    {"gene": "MLANA", "expression_pattern": "Specific to melanocytes"}\n'
                                        "  ]\n"
                                        "}\n\n"
                                        "Strictly follow this format. Do NOT include extra explanations, only return JSON."},
            {"role": "user", "content": f"I performed a DEGs analysis and clustering on my {dataset_type} data "
                                    f"from {sample_type} samples. "
                                    f"I found the following marker genes that are highly expressed in a specific cluster: "
                                    f"{marker_genes}. "
                                    "These genes are well-known markers for specific cell types. "
                                    "Predict the most likely cell type based on these markers. "
                                    "Provide a structured JSON response with a confidence score between 0.0 and 1.0, "
                                    "and include details on how each gene contributes to the prediction. "
                                    "Format the response as a JSON object with the following keys:\n\n"
                                    "- `cell_type`: A string indicating the predicted cell type.\n"
                                    "- `confidence_score`: A number between 0.0 and 1.0 representing confidence in the prediction.\n"
                                    "- `marker_gene_details`: A list of objects, each containing:\n"
                                    "    - `gene`: The marker gene name.\n"
                                    "    - `expression_pattern`: A description of its expression in the identified cell type.\n\n"
                                    "Return only valid JSON. No additional explanations or text."}
        ]

        try:
            completion = client.beta.chat.completions.parse(
                model="gpt-4o-2024-08-06",
                messages=messages,
                response_format={"type": "json_object"},  # Enforce structured output
                temperature=0  # ✅ Forces a deterministic response
            )

            # Extract and validate response
            #response_data = completion.model_dump() # Convert response opject to dictionary

            ### Modify 
            # Extract response content (string)
            response_content = completion.choices[0].message.content  # This is a string
            response_data = json.loads(response_content)  # ✅ Fix: Convert string to JSON object

            """if "choices" in response_data and response_data["choices"]:
                result = response_data["choices"][0]["message"]["content"]
                print("Response: ", result)

                # Save response as a JSON file
                with open("cell_type_response.json", "w") as json_file:
                    json.dump({"response": result}, json_file, indent=4)
                print("Response saved to 'cell_type_response.json'")

            else:
                print("Unexpected response format:", response_data)"""
            
            # Key normalization: Standardize keys for consistency
            """normalized_response = {
                "cell_type": response_data.get("choices", [{}])[0].get("message", {}).get("content", {}).get("cell_type", "Unknown"),
                "confidence_score": response_data.get("choices", [{}])[0].get("message", {}).get("content", {}).get("confidence_score", 0.0),
                "marker_gene_details": response_data.get("choices", [{}])[0].get("message", {}).get("content", {}).get("marker_gene_details", [])
            }"""
            ### Modify
            normalized_response = {
                "cell_type": response_data.get("cell_type", "Unknown"),
                "confidence_score": response_data.get("confidence_score", 0.0),
                "marker_gene_details": response_data.get("marker_gene_details", [])
            }

            # Print the structured response
            # print("Structured Response:", json.dumps(normalized_response, indent=4))

            # Save response as a JSON file
            


            json_file_name = f"core_{core}_cluster_{cluster}_cell.json"
            json_file_path = os.path.join(save_dir, json_file_name)

            with open(json_file_path, "w") as json_file:
                json.dump(normalized_response, json_file, indent=4)
            
            print(f"Response saved to '{json_file_name}.json'✅")
            
        except Exception as e:
            print(f"Error: {e}")

    print(f"core {core} Gene 추출 & Annotation 완료✅")
print("전체 core 진행 완료✅")
from openai import OpenAI
import os
import sys
import json
import anndata
import scanpy as sc

import warnings
warnings.filterwarnings("ignore")

# ✅ h5ad 파일 경로 직접 지정
file_path = sys.argv[1] if len(sys.argv) > 1 else input("h5ad 파일 경로 입력: ").strip()

if not os.path.exists(file_path):
    print(f"파일 없음: {file_path}")
    sys.exit(1)

# 저장 디렉토리: 입력 파일과 같은 위치에 annotation 폴더 생성
base_name = os.path.splitext(os.path.basename(file_path))[0]
save_dir = os.path.join(os.path.dirname(file_path), f"{base_name}_annotation")
os.makedirs(save_dir, exist_ok=True)

print(f"🔄 {base_name} Gene 추출 시작")
adata = anndata.read_h5ad(file_path)

sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon")
ranked_genes = adata.uns["rank_genes_groups"]
clusters = ranked_genes["names"].dtype.names

for cluster in clusters:
    marker_genes = ranked_genes["names"][cluster][:15]
    print(f"  Cluster {cluster}: {marker_genes}")

    client = OpenAI(api_key="")

    dataset_type = "single-cell RNA-seq"
    sample_type = "CRC patient colon Tissue"

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
            response_format={"type": "json_object"},
            temperature=0
        )

        response_content = completion.choices[0].message.content
        response_data = json.loads(response_content)

        normalized_response = {
            "cell_type": response_data.get("cell_type", "Unknown"),
            "confidence_score": response_data.get("confidence_score", 0.0),
            "marker_gene_details": response_data.get("marker_gene_details", [])
        }

        json_file_name = f"cluster_{cluster}_cell.json"
        json_file_path = os.path.join(save_dir, json_file_name)

        with open(json_file_path, "w") as json_file:
            json.dump(normalized_response, json_file, indent=4)

        print(f"  저장 완료: {json_file_name} ✅")

    except Exception as e:
        print(f"  Cluster {cluster} Error: {e}")

print(f"\n{base_name} Annotation 완료 ✅")
print(f"결과 저장 위치: {save_dir}")

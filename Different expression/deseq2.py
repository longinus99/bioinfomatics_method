# deseq2.py

import pandas as pd
import os
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
from rpy2.robjects.conversion import localconverter
from rpy2.robjects import default_converter

# rpy2 <-> pandas 자동 변환 활성화
pandas2ri.activate()
DESeq2 = importr("DESeq2")

def run_deseq2_from_tsv(
    tsv_path,
    group1="1",
    group2="2",
    logfc_thresh=1.0,
    padj_thresh=0.05,
    return_filtered=False,
    output_csv_path=None
):
    """
    TSV 파일에서 DESeq2 분석 실행
    :param tsv_path: TSV 파일 경로 (행: geneID, 열: 샘플명. e.g. '1-1', '2-3' 등)
    :param group1: 기준 그룹 (예: "1")
    :param group2: 비교 그룹 (예: "2")
    :param logfc_thresh: |log2FC| 필터링 기준
    :param padj_thresh: padj 필터링 기준
    :param return_filtered: 필터링된 DEG만 반환할지 여부
    :return: pandas DataFrame (DEG 결과)
    """
    # 1. TSV 파일 로딩
    df = pd.read_csv(tsv_path, sep="\t", index_col=0)
    df = df.round(0).astype(int)
    # 2. 그룹 추출 (e.g., "1-1" → "1")
    sample_names = df.columns.astype(str)
    groups = sample_names.str.split("-", n=1).str[0]
    groups = pd.Categorical(groups)
    coldata = pd.DataFrame({"group": groups}, index=sample_names)
    ro.r('rm(list=ls())')  # R 환경 초기화
    # 3. R로 전달 (수동 변환: Jupyter 안정성 ↑)
    with localconverter(default_converter + pandas2ri.converter):
        r_counts = pandas2ri.py2rpy(df)
        r_coldata = pandas2ri.py2rpy(coldata)

    ro.globalenv['countData'] = r_counts
    ro.globalenv['colData'] = r_coldata

    # 4. DESeq2 분석

    ro.r(f'''
        library(DESeq2)
        dds <- DESeqDataSetFromMatrix(countData = countData,
                                      colData = colData,
                                      design = ~ group)
        dds <- DESeq(dds)
        res <- results(dds, contrast = c("group", "{group2}", "{group1}"))
    ''')

    # 5. 결과 반환
    with localconverter(default_converter + pandas2ri.converter):
        res_df = (ro.r('as.data.frame(res)'))

    res_df.index = df.index

    if return_filtered:
        res_df = res_df[
            (res_df['padj'] < padj_thresh) &
            (res_df['log2FoldChange'].abs() > logfc_thresh)
        ]

    if output_csv_path is not None:
        os.makedirs(os.path.dirname(output_csv_path), exist_ok=True)
        res_df.reset_index().to_csv(output_csv_path, index=False)

    return res_df

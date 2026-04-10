import pandas as pd
import os
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects import default_converter

# rpy2 <-> pandas 변환 활성화
pandas2ri.activate()


def run_limma_from_csv(
    csv_path,
    case_label="Case",
    control_label="Control",
    group_col="SampleGroup",
    sample_col="SampleId",
    logfc_thresh=1.0,
    padj_thresh=0.05,
    return_filtered=False,
    output_csv_path=None
):
    """
    CSV 입력에서 limma 분석 실행 (Case vs Control)
    :param csv_path: CSV 파일 경로
    :param case_label: 환자군 라벨
    :param control_label: 대조군 라벨
    :param group_col: 그룹 지정 컬럼명
    :param sample_col: 샘플 ID 컬럼명
    :param logfc_thresh: |logFC| 필터링 기준
    :param padj_thresh: adj.P.Val 필터링 기준
    :param return_filtered: 필터링된 결과만 반환 여부
    :param output_csv_path: CSV 저장 경로
    """

    ro.r(f"""
        library(limma)

        # CSV 불러오기
        df <- read.csv("{csv_path}", header=TRUE, row.names=NULL, check.names=FALSE)
        df <- df[, colnames(df) != ""]
        # 메타 정보
        meta <- df[, c("{sample_col}", "{group_col}")]
        exprs <- as.matrix(df[, !(names(df) %in% c("{sample_col}", "{group_col}"))])
        rownames(exprs) <- df${sample_col}
        colnames(exprs) <- colnames(df)[!(colnames(df) %in% c("{sample_col}", "{group_col}"))]

        # 그룹 factor
        groups <- factor(meta[[ "{group_col}" ]], levels=c("{control_label}", "{case_label}"))

        # design matrix
        design <- model.matrix(~0 + groups)
        colnames(design) <- levels(groups)

        # limma 분석
        fit <- lmFit(t(exprs), design)
        contrast <- makeContrasts({case_label}-{control_label}, levels=design)
        fit2 <- contrasts.fit(fit, contrast)
        fit2 <- eBayes(fit2)
        res <- topTable(fit2, number=Inf, sort.by="none")

        # gene 심볼 추가
        res$Gene <- rownames(res)
        res <- res[, c("Gene", colnames(res)[colnames(res)!="Gene"])]
    """)

    with localconverter(default_converter + pandas2ri.converter):
        res_df = ro.r("as.data.frame(res)")

    if return_filtered:
        res_df = res_df[
            (res_df["adj.P.Val"] < padj_thresh) &
            (res_df["logFC"].abs() > logfc_thresh)
        ]

    if output_csv_path is not None:
        os.makedirs(os.path.dirname(output_csv_path), exist_ok=True)
        res_df.to_csv(output_csv_path, index=False)

    return res_df

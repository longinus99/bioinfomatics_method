import os
import pandas as pd
import scipy.io


# patient_id = "TA7047"
patient_id = "TA7048"
patient_num = "0041134" if patient_id == "TA7047" else "0053073"


base_path = "/home/longinus723/home"

region_folders = os.path.join(base_path, f"250210_CRC/{patient_id}")
region_folders = os.listdir(region_folders)


patient_output_path = os.path.join(base_path, "Demulti", patient_num)
os.makedirs(patient_output_path, exist_ok=True)


for region_folder in region_folders:
    region = region_folder.split('__')[2]
    region_num = region[-1]
    region_path = os.path.join(base_path, "core", f"{patient_id}_region")
    print(f"✅ {region} 처리 중...")

    if not os.path.exists(region_path):
        print(f"⚠️ {region_path} 없음. 스킵합니다.")
        continue

    # 10x 데이터 관련 경로
    cell_base_path = os.path.join(base_path, f"250323_CRC_BJM_2/{patient_id}")
    cell_matrix_path = os.path.join(
        cell_base_path,
        f"output-XETG00274__{patient_num}__{region}__20250319__044828",
        "cell_feature_matrix"
    )
    matrix_path = os.path.join(cell_matrix_path, "matrix.mtx.gz")
    barcodes_path = os.path.join(cell_matrix_path, "barcodes.tsv.gz")
    features_path = os.path.join(cell_matrix_path, "features.tsv.gz")
    cells_path = os.path.join(
        cell_base_path,
        f"output-XETG00274__{patient_num}__{region}__20250319__044828",
        "cells.csv.gz"
    )

    if not os.path.exists(matrix_path):
        print(f"⚠️ {matrix_path} 없음. 스킵합니다.")
        continue

    # (1) 전체 세포 정보 불러오기
    cells_df = pd.read_csv(cells_path, compression="gzip")
    total_cells_in_region = len(cells_df)  # region 내 전체 세포 개수

    # (2) 바코드 정보
    barcodes = pd.read_csv(barcodes_path, sep="\t", header=None)[0].values
    barcodes_set = set(barcodes)

    # (3) 이 region에서 분리된 코어의 총 세포 수를 세기 위한 변수
    region_core_cells_sum = 0

    coord_path = os.path.join(region_path, f"coordinates/{patient_id}_region{region_num}_coordinates.csv")

    if not os.path.exists(coord_path):
        print(f"⚠️ {coord_path} 없음. 스킵합니다.")
        continue 

    coord_df  = pd.read_csv(coord_path, skiprows=2)
    cores =coord_df ["Selection"].unique()

    # region 폴더 내 코어 좌표 파일 탐색
    for core in cores:
        print(f" {core} 처리 중...")
        core_coord  = coord_df [coord_df ["Selection"] == core]
        core_number = core.split("_")[0]  # 예: "1_coordinates.csv" → "1"
        core_id = f"core_{core_number}"
        core_output_path = os.path.join(patient_output_path, core_id)
        os.makedirs(core_output_path, exist_ok=True)

            
        # x, y 좌표의 최소/최댓값 추출
        x_min = core_coord ["X"].min()
        x_max = core_coord ["X"].max()
        y_min = core_coord ["Y"].min()
        y_max = core_coord ["Y"].max()

        # (5) cells_df 중 사각형 범위 내부인 세포 선별
        inside = (
            (cells_df["x_centroid"] >= x_min) &
            (cells_df["x_centroid"] <= x_max) &
            (cells_df["y_centroid"] >= y_min) &
            (cells_df["y_centroid"] <= y_max)
        )
        core_cells_df = cells_df[inside]
        core_cells_count = len(core_cells_df)  # 이 코어에 속하는 세포 수

        if core_cells_count == 0:
            print(f"⚠️ {region_folder} - Core {core_number}: 사각형 범위 내 cell 없음!")
            continue

        # region_core_cells_sum에 누적
        region_core_cells_sum += core_cells_count


        # (6) 10x 바코드 필터링
        matched_cell_ids = core_cells_df["cell_id"].values
        matched_barcodes = list(set(matched_cell_ids).intersection(barcodes_set))
        core_indices = [i for i, bc in enumerate(barcodes) if bc in matched_barcodes]

        if not core_indices:
            print(f"⚠️ {region_folder} - Core {core_number}: 매칭된 바코드 없음!")
            continue

        # (7) matrix.mtx.gz에서 코어 해당 세포만 추출
        matrix = scipy.io.mmread(matrix_path).tocsc()
        features = pd.read_csv(features_path, sep="\t", header=None)
        core_matrix = matrix[:, core_indices]

        # (8) 코어별 결과 저장
        scipy.io.mmwrite(os.path.join(core_output_path, "matrix.mtx"), core_matrix)
        pd.DataFrame(matched_barcodes).to_csv(
            os.path.join(core_output_path, "barcodes.tsv"),
            sep="\t", index=False, header=False
        )
        features.to_csv(
            os.path.join(core_output_path, "features.tsv"),
            sep="\t", index=False, header=False
        )

        # (9) 코어 처리 결과 출력 (코어 내 세포 수 표시)
        print(
            f"✅ {region_folder} - Core {core_number} 데이터 저장 완료! "
            f"({core_output_path}), cell count: {core_cells_count}"
        )

    # (10) region 전체 처리 결과: region_core_cells_sum vs total_cells_in_region
    print(
        f"🎉 {region_folder} 처리 완료! "
        f"코어 분리된 세포 합계: {region_core_cells_sum} / "
        f"전체 세포 수: {total_cells_in_region}"
    )

print(f"🎉 환자 `{patient_id}`의 모든 Region에 대한 Demultiplexing 완료!")

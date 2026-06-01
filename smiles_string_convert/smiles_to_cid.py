"""
SMILES → CID 변환 스크립트
모든 assay folds.json에서 SMILES 수집 후 PubChem API로 CID 조회
결과: smiles_cid_map.json
"""

import json
import time
import requests
from pathlib import Path
from urllib.parse import quote

# ---------------------- Config ---------------------- #

data_dir = Path('/home/longinus723/home/TOX-AI/comparing_gnns/01_assay_data_processing')
output_file = data_dir / 'smiles_cid_map.json'

assays = [
    'ATG_PXRE_CIS', 'LTEA_HepaRG_CYP2C19', 'LTEA_HepaRG_UGT1A1',
    'CCTE_Simmons_CellTiterGLO_HEK293T', 'CEETOX_H295R_OHPROG',
    'CLD_CYP1A1_48hr', 'NVS_ENZ_hBACE'
]

# ---------------------- SMILES 수집 ---------------------- #

all_smiles = set()
for assay in assays:
    fpath = data_dir / f'{assay}_folds.json'
    if not fpath.exists():
        print(f'⚠ {assay}_folds.json 없음, 스킵')
        continue
    folds = json.load(open(fpath))
    for fold in folds['folds'].values():
        all_smiles.update(fold['train'].keys())
        all_smiles.update(fold['test'].keys())

all_smiles = sorted(all_smiles)
print(f'총 SMILES 수: {len(all_smiles)}')

# ---------------------- 기존 결과 로드 (재시작 대비) ---------------------- #

if output_file.exists():
    smiles_cid_map = json.load(open(output_file))
    print(f'기존 결과 {len(smiles_cid_map)}개 로드')
else:
    smiles_cid_map = {}

# ---------------------- PubChem API 조회 ---------------------- #

def get_cid(smiles):
    try:
        url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/cids/JSON"
        r = requests.post(url, data={'smiles': smiles}, timeout=10)
        if r.status_code == 200:
            return r.json()['IdentifierList']['CID'][0]
    except Exception:
        pass
    return None

todo = [s for s in all_smiles if s not in smiles_cid_map]
print(f'조회 필요: {len(todo)}개')

for i, smiles in enumerate(todo):
    cid = get_cid(smiles)
    smiles_cid_map[smiles] = cid

    if (i + 1) % 50 == 0:
        json.dump(smiles_cid_map, open(output_file, 'w'))
        found = sum(1 for v in smiles_cid_map.values() if v is not None)
        print(f'[{i+1}/{len(todo)}] 저장 완료 | CID 찾음: {found}개')

    time.sleep(0.2)  # 초당 5요청 제한

# 최종 저장
json.dump(smiles_cid_map, open(output_file, 'w'))
found = sum(1 for v in smiles_cid_map.values() if v is not None)
print(f'\n완료! 총 {len(smiles_cid_map)}개 중 CID 찾음: {found}개 ({found/len(smiles_cid_map)*100:.1f}%)')
print(f'저장: {output_file}')

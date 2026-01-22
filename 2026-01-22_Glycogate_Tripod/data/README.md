# 3-Arm Tripod SMILES 구조

## 파일 설명

- **3arm_peg24_glc.smi** - 3-arm Tripod 구조의 SMILES (수소 제외)
- **3arm_peg24_glc_with_h.smi** - 수소가 추가된 SMILES (CGenFF 업로드용)
- **add_hydrogens.py** - RDKit을 사용하여 수소를 추가하는 스크립트

## 구조 설명

### TRIS 중심 (Tripod Core)
```
        Arm 1 (PEG24-Glucose)
         |
    N---C---Arm 2 (PEG24-Glucose)
         |
        Arm 3 (PEG24-Glucose)
```

### 각 팔 구성
- **PEG24:** 24개의 에틸렌 옥사이드 반복 단위
- **Glucose:** L-Glucose (입체화학 유지)
- **연결:** 에테르 결합 (O-glycosidic)

## 원자 수

- **수소 제외:** 410개 원자
- **수소 포함:** 889개 원자

## Bipod와의 비교

| 특성 | Bipod (2-arm) | Tripod (3-arm) |
|------|---------------|----------------|
| 팔 개수 | 2개 | 3개 |
| 원자 수 (H 제외) | ~273개 | 410개 |
| 원자 수 (H 포함) | ~593개 | 889개 |
| 대칭성 | C2 | C3 |
| 크기 | 중간 | 큰 |

## 입체화학

각 Glucose 단위는 L-형태를 유지:
- `[C@@H]` - S 입체중심
- `[C@H]` - R 입체중심

## 다음 단계

1. **CGenFF 업로드**
   - https://cgenff.umaryland.edu/
   - `3arm_peg24_glc_with_h.smi` 업로드
   - `lig.rtf`, `lig.prm` 다운로드

2. **Ligand Reader**
   - CHARMM-GUI Ligand Reader
   - 3D 구조 생성

3. **Solution Builder**
   - 수용액 시스템 구축

4. **MD 시뮬레이션**
   - 20ns Production MD

## 주의사항

- Tripod는 Bipod보다 크고 복잡하여 계산 비용이 더 높습니다
- 3개 팔의 대칭성을 유지하는 것이 중요합니다
- CUDA 가속 필수 (CPU로는 매우 느림)

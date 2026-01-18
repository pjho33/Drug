# SMILES 구조 분석

## 원본: Tris-(PEG2-β‑L‑Glc)3

```
NC(COCCOCCOO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H](O)[C@H]1O)(COCCOCCOO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H](O)[C@H]1O)COCCOCCOO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H](O)[C@H]1O
```

## 구조 분석

### TRIS 중심
- `NC(...)(...)(...)` - 중심 탄소에 3개 팔 + NH2

### 각 팔 구조 (3개 동일)
```
C - O - CC - O - CC - O - O - [Glucose]
│   │    │    │    │    │
팔  PEG  PEG  PEG  PEG  glycosidic bond
```

**PEG2 = 2개 EO 유닛:**
- `OCCOC` = -O-CH2-CH2-O-CH2-CH2-O-

**Glucose (β-L-Glc):**
```
O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H](O)[C@H]1O
```
- Pyranose ring (6원환)
- 5개 OH + 1개 CH2OH
- β-L 입체화학

---

## 1-Arm 버전: TRIS-PEG24-β‑L‑Glc

### 변경사항
1. **3개 팔 → 1개 팔**
   - 나머지 2개는 OH로 캡핑
   
2. **PEG2 → PEG24**
   - `OCCOC` (2 EO) → `OCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCC` (24 EO)
   - 24개 `-CH2-CH2-O-` 유닛

### 구조
```
        OH
        |
NH2 - C - CH2-O-(CH2-CH2-O)24-O-[Glucose]
        |
        OH
```

---

## PEG24 SMILES 생성

**PEG24 = 24개 EO 유닛:**
```
OCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCC
```

**확인:**
- `OCC` = -O-CH2-CH2- (1개 EO)
- 24번 반복 + 마지막 O

---

## 최종 1-Arm SMILES

```
NC(CO)(CO)COCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H](O)[C@H]1O
```

### 구조 설명
- `NC(CO)(CO)` - TRIS 중심 (NH2 + 2개 OH 팔)
- `COCCOC...` - PEG24 (24개 EO)
- `O[C@@H]1...` - β-L-Glucose

---

## 검증

### 원자 수 계산

**TRIS 부분:**
- C: 1 (중심)
- N: 1
- O: 2 (2개 OH)
- H: 계산 필요

**PEG24 부분:**
- C: 48 (24 × 2)
- O: 25 (24 + 1 glycosidic)
- H: 98 (24 × 4 + 2)

**Glucose 부분:**
- C: 6
- O: 6
- H: 12

**총합:**
- C: 1 + 48 + 6 = 55
- N: 1
- O: 2 + 25 + 6 = 33
- H: 계산 필요

---

## 다음 단계

1. SMILES 검증 (RDKit/OpenBabel)
2. 3D 구조 생성
3. CGenFF 파라미터 생성
4. Penalty score 확인

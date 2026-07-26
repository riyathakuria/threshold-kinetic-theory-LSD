# HuMicA v2.0.0 - Structure Report

Source file: `HuMicA_v2.0.0.h5ad`

## Dimensions

- **Cells (obs)**: 90,716
- **Genes in `.X` (var)**: 3,000
- **Genes in `.raw`**: 3,000
- **`.X` dtype**: float64

## Interpretation

- `.raw` present with 3,000 genes (not larger than `.X`); `.X` is not a strict HVG subset of `.raw`.

## `.var` (X gene annotation)

- index name: `None`; columns: ['features', 'SCT_features']
- sample gene IDs: ['ENSG00000152208', 'ENSG00000106366', 'ENSG00000118785', 'ENSG00000184226', 'ENSG00000124491', 'ENSG00000059804', 'ENSG00000183230', 'ENSG00000080824', 'ENSG00000148948', 'ENSG00000144724']

## `.raw.var` (raw gene annotation)

- index name: `None`; columns: ['_index']
- sample gene IDs: ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']

## `.obs` (cell metadata)

- 21 columns: ['TAG', 'orig.ident', 'nCount_RNA', 'nFeature_RNA', 'percent.mt', 'nCount_SCT', 'nFeature_SCT', 'Sample_ID', 'Source', 'Study', 'Group', 'Sex', 'Age', 'Tissue_type', 'Tissue_condition', 'Diagnosis', 'Death', 'Methodology', 'PMI', 'Tissue_group', 'V1_clusters']

### Annotation column value counts (<=60 categories)

- **TAG**: 90716 unique values (dtype object)
- **orig.ident**: HuMicA=90716
- **nCount_RNA**: 8216 unique values (dtype float64)
- **nFeature_RNA**: 3823 unique values (dtype int32)
- **percent.mt**: 49817 unique values (dtype float64)
- **nCount_SCT**: 4798 unique values (dtype float64)
- **nFeature_SCT**: 2682 unique values (dtype int32)
- **Sample_ID**: 241 unique values (dtype object)
- **Source**: 208 unique values (dtype object)
- **Study**: Mancuso 2019=18419, Franjic 2022=14951, Olah 2020=14263, Feleke 2021=7472, Lau 2020=5746, Leng 2021=4893, Thrupp 2020=3740, Morabito 2021=3597, Smajic 2022=3545, Tran 2021=3231, Zhou 2020=3159, Velmeshev 2019=2695 ...
- **Group**: No Neuropathology=33021, Epilepsy=24638, AD=20036, LBD=8877, COVID-19=1795, MS=1212, ASD=1137
- **Sex**: Male=45993, Female=30460, F=8461, M=5802
- **Age**: 98 unique values (dtype object)
- **Tissue_type**: Prefrontal Cortex=28660, Temporal Neocortex=24638, Anterior Cingulate Cortex=8462, Dentate Gyrus=6713, Hippocampus=5600, Entorhinal Cortex=4853, Midbrain=3545, Subiculum=2226, Superior Frontal Gyrus=1673, Medial Frontal Cortex=1666, Amygdala=968, White Matter=615 ...
- **Tissue_condition**: Postmortem=66078, Ex vivo=20285, Surgical Resection=4353
- **Diagnosis**: NA=23990, AD=9893, Epilepsy caused by brain tumor=9261, Therapy resistant epilepsy=6525, Epilepsy=6161, Mild Cognitive Impairment=5010, PD=4712, AD Braak 6=4136, DLB=2120, PDD=2045, Temporal Lobe Epilepsy=1866, AD TREM R62H=1311 ...
- **Death**: NA=70259, Myocardial infarction=8742, Atherosclerotic cardiovascular disease=4194, Cardiovascular disease, endstage renal failure=1448, Suicide, Hanging=778, Respiratory insufficiency and/or myocardial infarction=567, Interstitial pneumonia, multiorgan failure; arterial hypertonia, diabetes type II without insulin supplementation=333, Drowning=310, End organ failure due to COVID-19 pneumonia=298, Interstitial pneumonia, hypovolemic shock; arterial hypertonia, atrial fibrillation, pulmonary emphysema with barrel chest=262, Dilated cardiomegaly=253, Suicide, hanging=250 ...
- **Methodology**: nucleus=58034, cell=32682
- **PMI**: 89 unique values (dtype object)
- **Tissue_group**: Frontal_Parietal=32301, Limbic=29239, Temporal=25016, Midbrain=3545, White_matter=615
- **V1_clusters**: 0.0=25563, 1.0=21974, 2.0=16186, 3.0=8635, 4.0=8559, 5.0=3195, 6.0=2616, 7.0=2425, 8.0=1563

## Embeddings (`.obsm`)

- `X_pca`: shape [90716, 50]
- `X_umap`: shape [90716, 2]

## Layers

- ['SCT']

## `.uns` keys

- ['neighbors']

## Raw HDF5 top-level structure

```
/X: {'type': 'dataset', 'shape': (90716, 3000), 'dtype': 'float64'}
/layers: {'type': 'group', 'keys': ['SCT']}
/layers/SCT: {'type': 'dataset', 'shape': (90716, 3000), 'dtype': 'float64'}
/obs: {'type': 'group', 'keys': ['Age', 'Death', 'Diagnosis', 'Group', 'Methodology', 'PMI', 'Sample_ID', 'Sex', 'Source', 'Study', 'TAG', 'Tissue_condition', 'Tissue_group', 'Tissue_type', 'V1_clusters', '_index', 'nCount_RNA', 'nCount_SCT', 'nFeature_RNA', 'nFeature_SCT', 'orig.ident', 'percent.mt']}
/obs/Age: {'type': 'dataset', 'shape': (90716,), 'dtype': 'object'}
/obs/Death: {'type': 'dataset', 'shape': (90716,), 'dtype': 'object'}
/obs/Diagnosis: {'type': 'dataset', 'shape': (90716,), 'dtype': 'object'}
/obs/Group: {'type': 'dataset', 'shape': (90716,), 'dtype': 'object'}
/obs/Methodology: {'type': 'dataset', 'shape': (90716,), 'dtype': 'object'}
/obs/PMI: {'type': 'dataset', 'shape': (90716,), 'dtype': 'object'}
/obs/Sample_ID: {'type': 'dataset', 'shape': (90716,), 'dtype': 'object'}
/obs/Sex: {'type': 'dataset', 'shape': (90716,), 'dtype': 'object'}
/obs/Source: {'type': 'dataset', 'shape': (90716,), 'dtype': 'object'}
/obs/Study: {'type': 'dataset', 'shape': (90716,), 'dtype': 'object'}
/obs/TAG: {'type': 'dataset', 'shape': (90716,), 'dtype': 'object'}
/obs/Tissue_condition: {'type': 'dataset', 'shape': (90716,), 'dtype': 'object'}
/obs/Tissue_group: {'type': 'dataset', 'shape': (90716,), 'dtype': 'object'}
/obs/Tissue_type: {'type': 'dataset', 'shape': (90716,), 'dtype': 'object'}
/obs/V1_clusters: {'type': 'dataset', 'shape': (90716,), 'dtype': 'float64'}
/obs/_index: {'type': 'dataset', 'shape': (90716,), 'dtype': 'object'}
/obs/nCount_RNA: {'type': 'dataset', 'shape': (90716,), 'dtype': 'float64'}
/obs/nCount_SCT: {'type': 'dataset', 'shape': (90716,), 'dtype': 'float64'}
/obs/nFeature_RNA: {'type': 'dataset', 'shape': (90716,), 'dtype': 'int32'}
/obs/nFeature_SCT: {'type': 'dataset', 'shape': (90716,), 'dtype': 'int32'}
/obs/orig.ident: {'type': 'dataset', 'shape': (90716,), 'dtype': 'object'}
/obs/percent.mt: {'type': 'dataset', 'shape': (90716,), 'dtype': 'float64'}
/obsm: {'type': 'group', 'keys': ['X_pca', 'X_umap']}
/obsm/X_pca: {'type': 'dataset', 'shape': (90716, 50), 'dtype': 'float64'}
/obsm/X_umap: {'type': 'dataset', 'shape': (90716, 2), 'dtype': 'float64'}
/raw: {'type': 'group', 'keys': ['X', 'var']}
/raw/X: {'type': 'group', 'keys': ['data', 'indices', 'indptr']}
/raw/var: {'type': 'group', 'keys': ['_index']}
/uns: {'type': 'group', 'keys': ['neighbors']}
/uns/neighbors: {'type': 'group', 'keys': ['distances', 'params']}
/var: {'type': 'group', 'keys': ['SCT_features', '_index', 'features']}
/var/SCT_features: {'type': 'dataset', 'shape': (3000,), 'dtype': 'object'}
/var/_index: {'type': 'dataset', 'shape': (3000,), 'dtype': 'object'}
/var/features: {'type': 'dataset', 'shape': (3000,), 'dtype': 'object'}
/varm: {'type': 'group', 'keys': ['PCs']}
/varm/PCs: {'type': 'dataset', 'shape': (3000, 50), 'dtype': 'float64'}
```
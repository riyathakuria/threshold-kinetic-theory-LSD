# Threshold-kinetic-theory-LSD
Computational modeling and kinetic analysis of lysosomal storage disorder (Niemann-Pick Disease), focusing on sphingolipid metabolism thresholds. 

**Proteomic Refinement Funnel**
We executed a multi-stage filtration pipeline with biological relavance. _Proteomic_filtration_funnel.py_

**1. The Foundational Atlas**
The project began with the Brain Lysosome Atlas-microglia (BLA), a comprehensive dataset of microglial-specific lysosomal proteomics. This ensured that all downstream analysis was grounded in the unique metabolic profile of microglia rather than generic neuronal or systemic lysosomal data. (n=5353), _Count_unique_brain_lysosome_proteins.py_

**2. Systematic Literature Curation**
To ensure pathological relevance, the BLA was cross-referenced against a manually curated protein library (n = 1207) of 126 peer-reviewed studies strictly focused on Lysosomal Storage Disorders (LSDs). This manual verification step identified a high-confidence list of proteins that are  directly implicated in LSD pathogenesis.(n=647), _Brain lysosome atlas mapped to LSD proteins.py_


**3.Deep-Learning Palmitoylation Prediction Screening**
To identify post-translational modification sites, we screened the LSD candidates for S-palmitoylation:

Prediction: Using MusiteDeep, we performed deep-learning-based predictions on canonical UniProt sequences.

Validation: Candidates were verified against the SwissPalm database to differentiate between SwissProt Verified and Novel Prediction status.
_MusiteDeep_Palmitoylation_Only.py_

**4. SwissPalm Verified**
SwissPalm Validation: Predictions were annotated against the SwissPalm (Release 5) database. This categorized each protein as either:

SwissProt Verified: Experimentally confirmed in high-confidence databases.

Novel Prediction: Identified by MusiteDeep but not yet recorded in SwissPalm. _MusiteDeep+SwissProt_masterfile.py_

**5. Manual Protein-Protein Interaction using Stringdb v12.0**
   Network Construction: A Protein-Protein Interaction (PPI) network was constructed from palmitoylation-positive candidates (n=626).

**6.Spatiotemporal Spanning (Allen Brain Atlas)**
Spatiotemporal Mapping: Cluster expression was mapped across 21 developmental timepoints from 8 pcw to 11 years.

Regional Resolution: Focused on four key brain regions: Cerebellum, Hippocampus, Striatum, and Visual Cortex. 

Visualization: Row-wise z-score scaling heatmap.
_Spatiotemporal_mapping_of_PPI_candidates.py_, _ABA_heatmap.py_

For downstream analysis, the developmental timeline was grouped into 4 developmental groups. _Splitting_to_developmentalstages_postABA.py_

**7.Kinetic Stability**
Variance-Based Classification: Genes were categorized as Stable or Dynamic based on their expression variance across the developmental timeline.

Regional Ranking: Identified the Stable Genes (constitutively expressing) and Dynamic Genes (developmentally fluctuating) for each region and visualised by stacked bar plot. 
_Mapping_kinetic_stability_to_brainregions.py_, Stable_and_dynamic_genes_computingand_stackedbarplots

**8. Region-Based Prioritisation & MRI Mapping**
Mean expression and regional enrichment-based prioritisation of Stable Anchors and Dynamic Drivers._Computing_stableanchors_and_dynamicdrivers.py_

MRI mapping- Top1 stable anchor and top1 dynamic drivers were visualised using MNI coordinates. _MRI_viewer_stableanddynamic.py_

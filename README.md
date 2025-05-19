<div align="center">

# 🚀 `DisProtBench`:  Disorder-Aware, Task-Rich Benchmark for Evaluating Protein Structure Prediction in Realistic Biological Contexts

**Xinyue Zeng¹**, **Tuo Wang¹**, **Adithya Kulkarni¹**, **Alexander Lu¹**, **Alexandra Ni¹**, **Phoebe Xing¹**, **Junhan Zhao²³**, **Siwei Chen⁴⁵**, **Dawei Zhou¹**

¹ Virginia Tech, ² Harvard Medical School, ³ Harvard T.H. Chan School of Public Health, ⁴ Broad Institute of MIT and Harvard, ⁵ Massachusetts General Hospital

<!-- Stylish Buttons -->
<p>
  <img src="Figures/Benchmark.png" alt="DisProtBench" width="100%">
</p>

</div>

---

## 📌 Abstract
Recent advances in protein structure prediction have achieved near-atomic accuracy for well-folded proteins. However, current benchmarks inadequately assess model performance in biologically challenging contexts, especially those involving intrinsically disordered regions (IDRs), limiting their utility in applications like drug discovery, disease variant interpretation, and protein interface design. We introduce DisProtBench, a comprehensive benchmark for evaluating protein structure prediction models (PSPMs) under structural disorder and complex biological conditions. DisProtBench spans three key axes: 

(1) **Data complexity**—covering disordered regions, G protein-coupled receptors (GPCR)–ligand pairs, and multimeric complexes; 

(2) **Task diversity**—benchmarking twelve leading PSPMs across structure-based tasks with unified classification, regression, and interface metrics; 

(3) **Interpretability**—via the [DisProtBench Portal](http://zhoulab-1.cs.vt.edu:8501/), offering precomputed 3D structures and visual error analyses. 

Our results reveal significant variability in model robustness under disorder, with low-confidence regions linked to functional prediction failures. Notably, global accuracy metrics often fail to predict task performance in disordered settings, emphasizing the need for function-aware evaluation. DisProtBench establishes a reproducible, extensible, and biologically grounded framework for assessing next-generation PSPMs in realistic biomedical scenarios.

### 🔹 `DisProtBench` - A Unified Benchmark for IDR Investigation
We introduce DisProtBench with the following key contributions:

(1) **Database Development:** We curate a large benchmark dataset spanning biologically complex IDR scenarios, including thousands of disease-associated human proteins, GPCR–ligand interactions, and multimeric complexes with disorder-mediated interfaces. It captures structural heterogeneity essential for assessing model robustness in realistic contexts.

(2) **Task and Toolbox Development:** We introduce a unified evaluation toolbox to benchmark eleven PSPMs on disorder-sensitive tasks, using consistent metrics across PPI prediction, ligand binding, and contact mapping. Incorporating pLDDT-based stratification, DisProtBench uniquely isolates model behavior in ambiguous regions across tasks and model families.

(3) **Visual-Interactive Interface Development:** The DisProtBench Portal provides 3D visualizations, model comparison heatmaps, and interactive results to explore structure–function links, assess disorder-specific performance, and support hypothesis generation—without local setup.

## 📂 Datasets
We open-sourced our benchmark on [Kaggle](https://doi.org/10.34740/kaggle/ds/7400098), consisting of the following subsets:

| **Dataset**                  | **Description**                   | **# Number of Protein Only** | **Source**                                                                 |
|-----------------------------|-----------------------------------|---------------|---------------------------------------------------------------------------|
| **DisProt-Based Dataset**   | Disorder in human disease         | 3,060         | First proposed in our work                                                |
| **Protein Interaction Dataset** | Disorder-mediated interfaces      | 1,200         | [GitHub](https://github.com/ohuelab/SpatialPPI/tree/main)                |
| **Individual Protein Dataset** | Disorder and ligand binding       | 20            | [GitHub](https://github.com/ChengF-Lab/LISA-CPI?tab=readme-ov-file)       |
              
## 🏗️ Toolbox
### 📂 Models Toolbox
We benchmark state-of-the-art PSPMs spanning diverse architectures, inputs, and structural representations across protein-related tasks, as summarized below:

| **PSPM**       | **Task**       | **Architecture**     | **Input**       | **Source**         | **Structural Representation**                |
|----------------|----------------|-----------------------|------------------|---------------------|------------------------|
| AF2            | PPI, Drug      | Evoformer             | MSA              | [Paper](https://www.nature.com/articles/s41586-021-03819-2)            | Atomic                 |
| AF3            | PPI, Drug      | Evoformer+LLM         | MSA + Seq        | [Paper](https://www.nature.com/articles/s41586-024-07487-w)            | Atomic + ligand        |
| OpenFold       | PPI, Drug      | Evoformer             | MSA              | [Paper](https://www.nature.com/articles/s41592-024-02272-z)               | Atomic                 |
| UniFold        | PPI, Drug      | Evoformer             | MSA              | [Paper](https://www.biorxiv.org/content/10.1101/2022.08.04.502811v3)        | Atomic                 |
| Boltz          | PPI, Drug      | Transformer           | Seq-only         | [Paper](https://www.biorxiv.org/content/10.1101/2024.11.19.624167v4)      | Coarse-grained         |
| Chai           | PPI, Drug      | Transformer           | Seq-only         | [Paper](https://www.biorxiv.org/content/10.1101/2024.10.10.615955v1)      | Coarse-grained         |
| Protenix       | PPI, Drug      | Transformer+          | Seq-only         | [Paper](https://www.biorxiv.org/content/10.1101/2025.01.08.631967v1)           | Atomic + ligand        |
| ESMFold        | Drug           | Transformer           | Seq-only         | [Paper](https://www.biorxiv.org/content/10.1101/2022.07.20.500902v3)             | Coarse-grained         |
| OmegaFold      | Drug           | Transformer           | Seq-only         | [Paper](https://www.biorxiv.org/content/10.1101/2022.07.21.500999v1)             | Coarse-grained         |
| RoseTTAFold    | Drug           | Hybrid (CNN+Attn)     | MSA              | [Paper](https://pubmed.ncbi.nlm.nih.gov/34282049/)           | Atomic                 |
| DeepFold       | Drug           | Custom DL             | Seq-only         | [Paper](https://pubmed.ncbi.nlm.nih.gov/37995286/)       | Atomic                 |

### 📊 Evaluation ToolBox
We evaluate model performance using a comprehensive set of classification, regression, and structural interface metrics, defined as follows:

| **Metric**                        | **Definition / Formula**                                                                                                           |
|----------------------------------|-------------------------------------------------------------------------------------------------------------------------------------|
| **Classification Metrics**       |                                                                                                                                     |
| Precision (Positive Predictive Value) | $\displaystyle \frac{\text{TP}}{\text{TP} + \text{FP}}$                                                                      |
| Recall (Sensitivity)             | $\displaystyle \frac{\text{TP}}{\text{TP} + \text{FN}}$                                                                      |
| F1 Score                         | $\displaystyle \frac{2 \cdot \text{TP}}{2 \cdot \text{TP} + \text{FP} + \text{FN}}$                                           |
| Accuracy                         | $\displaystyle \frac{\text{TP} + \text{TN}}{\text{TP} + \text{TN} + \text{FP} + \text{FN}}$                                   |
| **Regression Metrics**           |                                                                                                                                     |
| Mean Absolute Error (MAE)        | $\displaystyle \frac{1}{N} \sum_{i=1}^{N} \left| y_i - \hat{y}_i \right|$                                                     |
| Mean Squared Error (MSE)         | $\displaystyle \frac{1}{N} \sum_{i=1}^{N} \left( y_i - \hat{y}_i \right)^2$                                                    |
| Pearson Correlation Coefficient ($R$) | $\displaystyle \frac{ \sum_{i=1}^{N} (y_i - \bar{y})(\hat{y}_i - \bar{\hat{y}}) }{ \sqrt{ \sum_{i=1}^{N} (y_i - \bar{y})^2 } \cdot \sqrt{ \sum_{i=1}^{N} (\hat{y}_i - \bar{\hat{y}})^2 } }$ |
| **Structural Interface Metrics** |                                                                                                                                     |
| Receptor Precision (RP)          | $\displaystyle \frac{|\text{True Receptor Interface} \cap \text{Predicted Receptor Interface}|}{|\text{Predicted Receptor Interface}|}$ |
| Receptor Recall (RR)             | $\displaystyle \frac{|\text{True Receptor Interface} \cap \text{Predicted Receptor Interface}|}{|\text{True Receptor Interface}|}$    |
| Ligand Precision (LP)            | $\displaystyle \frac{|\text{Predicted Ligand Interface} \cap \text{True Receptor Interface}|}{|\text{Predicted Ligand Interface}|}$   |
| Ligand Recall (LR)               | $\displaystyle \frac{|\text{True Ligand Interface} \cap \text{Predicted Receptor Interface}|}{|\text{True Ligand Interface}|}$        |

For the definitions of Receptor and Ligand, we follow the work of [*Multi-level analysis of intrinsically disordered protein docking methods*](https://pubmed.ncbi.nlm.nih.gov/35609776/):

<p align="center"> <img src="Figures/ligand_receptor.jpg" alt="Definition" width="100%"> </p>

## 🎨 Visualize Portal
For more visualizations please link to 

🔗 **Visual-Interactive Interface**: [DisProtBench](http://zhoulab-1.cs.vt.edu:8501/) 

<p align="center"> <img src="Figures/UI.png" alt="Visualization Examples" width="100%"> </p>

## 🏗️ Results
### PPI Prediction
We use PPI prediction as a downstream task for evaluating PSPMs with the following pipeline:

<p align="center"> <img src="Figures/PPI.png" alt="PPI Pipeline" width="100%"> </p>

We examine the robustness of PSPMs in predicting PPI, a setting where disordered regions frequently mediate transient or flexible binding. The results are as follows:

|          | **Original** |          |         |        | **pLDDT ≥ 30** |          |         |        | **pLDDT ≥ 50** |          |         |        |
| -------- | ------------ | -------- | ------- | ------ | -------------- | -------- | ------- | ------ | -------------- | -------- | ------- | ------ |
| **PSPM** | **Acc**      | **Prec** | **Rec** | **F1** | **Acc**        | **Prec** | **Rec** | **F1** | **Acc**        | **Prec** | **Rec** | **F1** |
| AF2      | 0.793        | 0.783    | 0.799   | 0.791  | 0.802          | 0.791    | 0.812   | 0.801  | 0.818          | 0.809    | 0.825   | 0.817  |
| AF3      | 0.900        | 0.888    | 0.915   | 0.901  | 0.905          | 0.893    | 0.905   | 0.906  | 0.913          | 0.899    | 0.930   | 0.914  |
| Boltz    | 0.850        | 0.848    | 0.853   | 0.850  | 0.858          | 0.853    | 0.863   | 0.858  | 0.869          | 0.870    | 0.868   | 0.869  |
| Chai     | 0.850        | 0.841    | 0.863   | 0.852  | 0.858          | 0.847    | 0.873   | 0.860  | 0.869          | 0.857    | 0.887   | 0.871  |
| OpenFold | 0.624        | 0.605    | 0.605   | 0.605  | 0.643          | 0.622    | 0.638   | 0.630  | 0.671          | 0.656    | 0.651   | 0.653  |
| Proteinx | 0.810        | 0.809    | 0.812   | 0.810  | 0.819          | 0.820    | 0.818   | 0.819  | 0.834          | 0.834    | 0.835   | 0.834  |
| UniFold  | 0.552        | 0.378    | 0.667   | 0.483  | 0.567          | 0.389    | 0.667   | 0.491  | 0.597          | 0.417    | 0.714   | 0.526  |

<p align="center"> <img src="Figures/overall_pvalue_heatmap_PPI.png" alt="PPI Heatmap" width="100%"> </p>

### Drug Discovery
We use drug discovery as a downstream task for evaluating PSPMs with the following pipeline:

<p align="center"> <img src="Figures/Drug.png" alt="Drug Discovery Pipeline" width="100%"> </p>

We examine the robustness of PSPMs in discovering drug, a setting where disordered regions frequently mediate transient or flexible binding. The results are as follows:

|             | **Original** |       | **pLDDT ≥ 30** |       | **pLDDT ≥ 50** |       |
| ----------- | ------------ | ----- | -------------- | ----- | -------------- | ----- |
| **Model**   | **MAE**      | **R** | **MAE**        | **R** | **MAE**        | **R** |
| AlphaFold3  | 0.048        | 0.999 | 0.048          | 0.999 | 0.049          | 0.999 |
| Proteinx    | 0.072        | 0.997 | 0.072          | 0.997 | 0.072          | 0.997 |
| Boltz       | 0.079        | 0.996 | 0.079          | 0.996 | 0.080          | 0.996 |
| Chai        | 0.096        | 0.995 | 0.096          | 0.995 | 0.097          | 0.995 |
| DeepFold    | 0.144        | 0.988 | 0.144          | 0.988 | 0.144          | 0.988 |
| ESMFold     | 0.150        | 0.987 | 0.150          | 0.987 | 0.151          | 0.987 |
| OmegaFold   | 0.151        | 0.987 | 0.151          | 0.987 | 0.152          | 0.987 |
| OpenFold    | 0.151        | 0.987 | 0.151          | 0.987 | 0.151          | 0.987 |
| AlphaFold2  | 0.160        | 0.985 | 0.160          | 0.985 | 0.160          | 0.985 |
| UniFold     | 0.183        | 0.981 | 0.183          | 0.981 | 0.184          | 0.981 |
| RoseTTAFold | 0.190        | 0.979 | 0.190          | 0.979 | 0.190          | 0.979 |

<p align="center"> <img src="Figures/overall_pvalue_heatmap_Drug.png" alt="Drug Heatmap" width="100%"> </p>

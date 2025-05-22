# DisProtBench: Benchmarking Protein-Protein Interactions and Drug Binding Predictions

This repository contains the data and code used for our manuscript on benchmarking protein-protein interactions (PPI) and drug binding predictions. The repository is organized into two main directories: `PPI/` and `Drug/`, each containing analysis scripts, data, and generated figures.

## Repository Structure

### PPI Directory
The `PPI/` directory contains analysis and visualization code for protein-protein interaction predictions:

- **Scripts:**
  - `barplot_model.py`: Generates bar plots comparing different model performances
  - `barplot_boltz_vertical.py`: Creates vertical bar plots for Boltzmann-based metrics
  - `heatmap.py`: Generates heatmaps for p-value analysis

- **Generated Figures:**
  - `boltz_metrics_vertical_barplot.png`: Vertical bar plot of Boltzmann metrics
  - `overall_pvalue_heatmap.png`: Heatmap visualization of p-values
  - `ligand_metrics_overall_barplot.png`: Bar plot of ligand binding metrics

- **Data Directories:**
  - `Full/`: Contains complete dataset analysis
  - `Disorder/`: Contains disorder-specific analysis

### Drug Directory
The `Drug/` directory contains analysis and visualization code for drug binding predictions:

- **Scripts:**
  - `pvalue.py`: Calculates and visualizes p-values for drug binding predictions
  - `MAE.py`: Computes Mean Absolute Error metrics
  - `MSE.py`: Computes Mean Squared Error metrics

- **Generated Figures:**
  - `overall_pvalue_heatmap.png`: Heatmap of p-values across different conditions
  - `MAE_R_overall_horizontal.png`: Horizontal bar plot of MAE metrics
  - `MAE_R_overall.png`: Overall MAE visualization
  - `MSE_barplot.png`: Bar plot of MSE metrics

- **Data Files:**
  - `MSE_train.csv`: Training MSE data for full dataset
  - `MSE_train_30.csv`: Training MSE data for 30% disorder threshold
  - `MSE_train_50.csv`: Training MSE data for 50% disorder threshold

- **Data Directories:**
  - `Prediction_LISA/`: Contains LISA-based prediction results
  - `full/`: Complete dataset analysis
  - `disorder_30/`: Analysis for 30% disorder threshold
  - `disorder_50/`: Analysis for 50% disorder threshold

## Usage

To reproduce the analyses and figures:

1. Ensure you have Python installed with required dependencies
2. Navigate to the respective directory (`PPI/` or `Drug/`)
3. Run the Python scripts to generate figures and analyses

## Citation

If you use this code or data in your research, please cite our manuscript:

[Citation information to be added]

## Contact

For questions or issues, please contact [Contact information to be added] 
# DisProtBench Downstream Evaluation Toolbox

This toolbox provides a unified interface for downstream evaluation tasks, including Protein-Protein Interaction (PPI) and Compound-Protein Interaction (CPI) prediction.

## Directory Structure

- `toolbox.py` : Unified command-line interface for all tasks
- `ppi/` : PPI prediction codebase (TensorFlow, 3D CNNs)
- `cpi/` : CPI/drug discovery codebase (PyTorch, GPCR, AlphaFold, etc.)
- `data/` : Datasets
- `scripts/` : Utility scripts

---

## Unified CLI Usage (`toolbox.py`)

All major PPI and CPI tasks can be run from a single entry point:

### General format
```bash
python toolbox.py <ppi|cpi> <subcommand> [options]
```

### PPI Tasks

#### Train a PPI model
```bash
python toolbox.py ppi train \
  --model DenseNet3D \
  --datapath ./data/example_dataset/distance \
  --train_set ./data/example_dataset/part_0_train.csv \
  --test_set ./data/example_dataset/part_0_val.csv \
  --savingPath ./models/example_model
```

#### Test a PPI model
```bash
python toolbox.py ppi test \
  --model DenseNet3D \
  --datapath ./data/example_dataset/distance \
  --weights ./models/example_model/<timestamp>/best_model.weights.h5 \
  --output ./models/example_model/<timestamp>/preds.npy \
  --test_set ./data/example_dataset/part_0_test.csv
```

### CPI/Drug Discovery Tasks

#### Prepare CPI dataset
```bash
python toolbox.py cpi prepare_data \
  --dataset data/original/top20_raw.csv \
  --gpcr-col uniprot_id \
  --smiles-col smiles \
  --label-col pKi \
  --file-name-col inchi_key \
  --rep-path data/representations/top20/{}.npy \
  --save-path data/ligands/top20/imgs \
  --anno-path data/ligands/top20/anno \
  -j 12 \
  --test-size 0.3 \
  --task regression
```

#### Train a CPI model
```bash
python toolbox.py cpi train --cfg configs/train/top20.yml
```

#### Predict with a CPI model
```bash
python toolbox.py cpi predict \
  --cfg configs/prediction/pain.yml \
  --data-dir data/pred/fda \
  --rep-path data/representations/pain/P08908.npy \
  --out-dir output/prediction/fda/P08908
```

---

## References
- [LISA-CPI GitHub](https://github.com/ChengF-Lab/LISA-CPI?tab=readme-ov-file)
- [SpatialPPI GitHub](https://github.com/ohuelab/SpatialPPI/tree/main)
- [AlphaFold GitHub](https://github.com/deepmind/alphafold)

---

**Tip:**
You can further automate workflows by writing additional wrapper scripts in `toolbox/scripts/` for custom tasks if needed. 
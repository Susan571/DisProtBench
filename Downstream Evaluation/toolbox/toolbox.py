import argparse
import subprocess
import sys
import os

TOOLBOX_ROOT = os.path.dirname(os.path.abspath(__file__))

def run_command(cmd, cwd=None):
    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, cwd=cwd)
    if result.returncode != 0:
        sys.exit(result.returncode)

def ppi_train(args):
    cmd = [
        "python", "train.py",
        "--model", args.model,
        "--datapath", args.datapath,
        "--train_set", args.train_set,
        "--test_set", args.test_set,
        "--savingPath", args.savingPath
    ]
    run_command(cmd, cwd=os.path.join(TOOLBOX_ROOT, "ppi"))

def ppi_test(args):
    cmd = [
        "python", "test.py",
        "--model", args.model,
        "--datapath", args.datapath,
        "--weights", args.weights,
        "--output", args.output,
        "--test_set", args.test_set
    ]
    run_command(cmd, cwd=os.path.join(TOOLBOX_ROOT, "ppi"))

def cpi_prepare_data(args):
    cmd = [
        "python", "scripts/prepare_dataset.py",
        "--dataset", args.dataset,
        "--gpcr-col", args.gpcr_col,
        "--smiles-col", args.smiles_col,
        "--label-col", args.label_col,
        "--file-name-col", args.file_name_col,
        "--rep-path", args.rep_path,
        "--save-path", args.save_path,
        "--anno-path", args.anno_path,
        "-j", str(args.workers),
        "--test-size", str(args.test_size),
        "--task", args.task
    ]
    run_command(cmd, cwd=os.path.join(TOOLBOX_ROOT, "cpi"))

def cpi_train(args):
    cmd = [
        "python", "scripts/train.py",
        "--cfg", args.cfg
    ]
    run_command(cmd, cwd=os.path.join(TOOLBOX_ROOT, "cpi"))

def cpi_predict(args):
    cmd = [
        "python", "scripts/prediction.py",
        "--cfg", args.cfg,
        "--data-dir", args.data_dir,
        "--rep-path", args.rep_path,
        "--out-dir", args.out_dir
    ]
    run_command(cmd, cwd=os.path.join(TOOLBOX_ROOT, "cpi"))

def main():
    parser = argparse.ArgumentParser(description="Unified Toolbox for PPI and CPI tasks")
    subparsers = parser.add_subparsers(dest="command")

    # PPI subcommands
    ppi_parser = subparsers.add_parser("ppi", help="PPI tasks")
    ppi_subparsers = ppi_parser.add_subparsers(dest="subcommand")

    ppi_train_parser = ppi_subparsers.add_parser("train", help="Train PPI model")
    ppi_train_parser.add_argument("--model", required=True)
    ppi_train_parser.add_argument("--datapath", required=True)
    ppi_train_parser.add_argument("--train_set", required=True)
    ppi_train_parser.add_argument("--test_set", required=True)
    ppi_train_parser.add_argument("--savingPath", required=True)
    ppi_train_parser.set_defaults(func=ppi_train)

    ppi_test_parser = ppi_subparsers.add_parser("test", help="Test PPI model")
    ppi_test_parser.add_argument("--model", required=True)
    ppi_test_parser.add_argument("--datapath", required=True)
    ppi_test_parser.add_argument("--weights", required=True)
    ppi_test_parser.add_argument("--output", required=True)
    ppi_test_parser.add_argument("--test_set", required=True)
    ppi_test_parser.set_defaults(func=ppi_test)

    # CPI subcommands
    cpi_parser = subparsers.add_parser("cpi", help="CPI tasks")
    cpi_subparsers = cpi_parser.add_subparsers(dest="subcommand")

    cpi_prepare_parser = cpi_subparsers.add_parser("prepare_data", help="Prepare CPI dataset")
    cpi_prepare_parser.add_argument("--dataset", required=True)
    cpi_prepare_parser.add_argument("--gpcr-col", required=True)
    cpi_prepare_parser.add_argument("--smiles-col", required=True)
    cpi_prepare_parser.add_argument("--label-col", required=True)
    cpi_prepare_parser.add_argument("--file-name-col", required=True)
    cpi_prepare_parser.add_argument("--rep-path", required=True)
    cpi_prepare_parser.add_argument("--save-path", required=True)
    cpi_prepare_parser.add_argument("--anno-path", required=True)
    cpi_prepare_parser.add_argument("-j", "--workers", type=int, default=12)
    cpi_prepare_parser.add_argument("--test-size", type=float, default=0.3)
    cpi_prepare_parser.add_argument("--task", default="regression")
    cpi_prepare_parser.set_defaults(func=cpi_prepare_data)

    cpi_train_parser = cpi_subparsers.add_parser("train", help="Train CPI model")
    cpi_train_parser.add_argument("--cfg", required=True)
    cpi_train_parser.set_defaults(func=cpi_train)

    cpi_predict_parser = cpi_subparsers.add_parser("predict", help="Predict with CPI model")
    cpi_predict_parser.add_argument("--cfg", required=True)
    cpi_predict_parser.add_argument("--data-dir", required=True)
    cpi_predict_parser.add_argument("--rep-path", required=True)
    cpi_predict_parser.add_argument("--out-dir", required=True)
    cpi_predict_parser.set_defaults(func=cpi_predict)

    args = parser.parse_args()
    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()

if __name__ == "__main__":
    main() 
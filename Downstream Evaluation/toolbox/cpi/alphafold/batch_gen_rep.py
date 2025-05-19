import os

def parse_args():

    import argparse

    parser = argparse.ArgumentParser(
        description='batch generate intermediate representations'
    )
    parser.add_argument('--data-dir', default='/home/user/GPCR/alphafold/data', type=str, help='absolute path to alphafold data')
    parser.add_argument('--fasta-path', default='/home/user/GPCR/alphafold/fastas', type=str, help='absolute path to folder include fasta files')
    parser.add_argument('--out-dir', default='/home/user/GPCR/alphafold/output', type=str, help='absolute output path')

    args = parser.parse_args()

    return args

def main(args):

    print(args)

    fasta_list = [os.path.join(args.fasta_path, fasta) for fasta in os.listdir(args.fasta_path)]
    fasta_paths = ','.join(fasta_list)
    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir)
    
    cmd = 'python docker/run_docker.py \
    --fasta_paths={} \
    --max_template_date=2021-11-01 \
    --data_dir={} \
    --output_dir={} \
    --docker_image_name=alphafold_test'.format(fasta_paths, args.data_dir, args.out_dir)

    os.system(cmd)

if __name__ == '__main__':

    args = parse_args()
    main(args)
import difflib

def compare_files(file1_path, file2_path):
    with open(file1_path, 'r') as file1, open(file2_path, 'r') as file2:
        file1_lines = file1.readlines()
        file2_lines = file2.readlines()

    differ = difflib.Differ()
    diff = list(differ.compare(file1_lines, file2_lines))

    differences = []
    for line in diff:
        if line.startswith('- ') or line.startswith('+ ') or line.startswith('? '):
            differences.append(line)

    return differences

def main():
    file1_path = '/home/xyzeng/Data/Uniprot/final_selected_data.txt'
    file2_path = '/home/xyzeng/Data/Uniprot/final_selected_data_2.txt'

    differences = compare_files(file1_path, file2_path)

    if not differences:
        print("The files are identical.")
    else:
        print("The files are different. Here are the differences:")
        for diff in differences:
            print(diff)

if __name__ == "__main__":
    main()
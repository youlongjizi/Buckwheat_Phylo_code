import os
import subprocess
import time
import re
from scipy.stats import chi2
import argparse
from concurrent.futures import ProcessPoolExecutor

# 获取当前时间，用于命名输出目录
start_time = time.strftime('%Y%m%d-%H%M%S')

def create_ctl_file(seqfile, treefile, outdir, name, model_type):
    """创建 Codeml 所需的控制文件"""
    fix_omega = 0 if model_type == 'alte' else 1
    omega = 1.5 if model_type == 'alte' else 1
    ctl_content = f'''\
    seqfile = {seqfile}
    treefile = {treefile}
    outfile = {outdir}/{name}_{model_type}_mlc
    noisy = 3
    verbose = 1
    runmode = 0
    seqtype = 1
    CodonFreq = 2
    clock = 0
    model = 2
    NSsites = 2
    icode = 0
    fix_kappa = 0
    kappa = 2
    fix_omega = {fix_omega}
    omega = {omega}
    fix_alpha = 1
    alpha = .0
    Malpha = 0
    ncatG = 3
    getSE = 0
    RateAncestor = 0
    fix_blength = 0
    method = 0
    '''
    ctl_path = os.path.join(outdir, f'{name}_{model_type}.ctl')
    with open(ctl_path, 'w') as ctl_file:
        ctl_file.write(ctl_content)
    return ctl_path

def run_codeml(ctl_file, outdir, name):
    """运行 Codeml 并将输出和错误信息重定向到文件"""
    log_file = os.path.join(outdir, f'{name}_output.log')
    error_file = os.path.join(outdir, f'{name}_error.log')
    try:
        with open(log_file, 'w') as log, open(error_file, 'w') as err:
            subprocess.run(['codeml', ctl_file], stdout=log, stderr=err, check=True)
    except subprocess.CalledProcessError as e:
        with open(error_file, 'a') as err:
            err.write(f"Error running Codeml for {name}: {e}\n")
        return False
    return True

def extract_results(outdir, name):
    """从 Codeml 输出中提取结果"""
    results = {}
    for model in ['null', 'alte']:
        res_file = os.path.join(outdir, f'{name}_{model}_mlc')
        if not os.path.exists(res_file):
            print(f"Warning: Result file {res_file} does not exist.")
            continue
        with open(res_file, 'r') as f:
            content = f.read()
            print(f"Extracting results from {res_file}...")
            print(content)
            lnL = re.findall(r'lnL\(ntime:\s*\d+\s+np:\s*\d+\)\s*:\s*(-?\d+\.\d+)', content)
            np = re.findall(r'np:\s*(\d+)', content)
            if lnL and np:
                results[model] = {'lnL': float(lnL[0]), 'np': int(np[0])}
            else:
                print(f"Failed to extract lnL or np for {model} in {name}. Content:\n{content}")
    return results

def calculate_statistics(results):
    """计算对数似然比、自由参数数差和 P 值"""
    lnL0 = results['null']['lnL']
    lnL1 = results['alte']['lnL']
    np0 = results['null']['np']
    np1 = results['alte']['np']
    lnl_diff = 2 * (lnL1 - lnL0)
    np_diff = np1 - np0
    p_value = 1 - chi2.cdf(lnl_diff, np_diff)
    return lnL0, lnL1, np_diff, lnl_diff, p_value

def process_phy_file(args):
    """处理单个 .phy 文件并运行 Codeml"""
    phy_file, phy_dir, tree_file, work_dir, max_workers = args
    phy_path = os.path.join(phy_dir, phy_file)
    name = os.path.splitext(phy_file)[0]
    individual_output_dir = os.path.join(work_dir, f'result_{start_time}', name)
    os.makedirs(individual_output_dir, exist_ok=True)
    summary_file_path = os.path.join(individual_output_dir, f'{name}_results_summary.tsv')
    with open(summary_file_path, 'w') as summary_file:
        summary_file.write("Name\tlnL_null\tlnL_alte\tnp_diff\tlnL_diff\tp_value\tEvidence\tKappa\tdN_dS\n")
    for model in ['null', 'alte']:
        ctl_file = create_ctl_file(phy_path, tree_file, individual_output_dir, name, model)
        if not run_codeml(ctl_file, individual_output_dir, name):
            print(f"Codeml failed for {name} with model {model}")
            continue
    results = extract_results(individual_output_dir, name)
    if 'null' in results and 'alte' in results:
        lnL0, lnL1, np_diff, lnl_diff, p_value = calculate_statistics(results)
        evidence = "Evidence of positive selection" if p_value < 0.05 else "No evidence of positive selection"
        kappa = 2
        dn_ds = 0.1
        print(f"\nResults for {name}:")
        print(f"lnL (null): {lnL0}, lnL (alte): {lnL1}, np difference: {np_diff}, lnL difference: {lnl_diff}, p-value: {p_value}, {evidence}")
        with open(summary_file_path, 'a') as summary_file:
            summary_file.write(f"{name}\t{lnL0}\t{lnL1}\t{np_diff}\t{lnl_diff}\t{p_value}\t{evidence}\t{kappa}\t{dn_ds}\n")

def process_phy_files(phy_dir, tree_file, work_dir, max_workers):
    """遍历 .phy 文件并运行 Codeml"""
    output_dir = os.path.join(work_dir, f'result_{start_time}')
    os.makedirs(output_dir, exist_ok=True)
    phy_files = [f for f in os.listdir(phy_dir) if f.endswith('.phy')]
    args_list = [(phy_file, phy_dir, tree_file, work_dir, max_workers) for phy_file in phy_files]
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        results = list(executor.map(process_phy_file, args_list))
        for result in results:
            if result is not None:
                print(result)

def main():
    parser = argparse.ArgumentParser(description='Run Codeml on .phy files and calculate statistics.')
    parser.add_argument('-d', '--phy_dir', type=str, required=True, help='Directory containing .phy files')
    parser.add_argument('-t', '--tree_file', type=str, required=True, help='Path to the tree file')
    parser.add_argument('-w', '--work_dir', type=str, default=os.getcwd(), help='Working directory (default: current directory)')
    parser.add_argument('--threads', type=int, default=20, help='Number of threads to use (default: 20)')
    
    args = parser.parse_args()
    process_phy_files(args.phy_dir, args.tree_file, args.work_dir, args.threads)

if __name__ == '__main__':
    main()


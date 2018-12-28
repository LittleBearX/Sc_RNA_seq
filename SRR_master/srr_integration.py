__doc__ = """
            INPUT：SRR_INFO、SRR_DATA
            PROCESS：按GSE整合SRR
            OUTPUT：GSE_DATA
            """

import os
import sys
import getopt
import pandas as pd
from functools import partial
from multiprocessing import Pool
from Sc_RNA_seq.helper import log
from Sc_RNA_seq.helper import make_dir
from Sc_RNA_seq.helper import add_color
from Sc_RNA_seq.SRR_master.srr_helper import distribute_srr
from Sc_RNA_seq.SRR_master.srr_helper import integration_worker

# 读取配置文件
from Sc_RNA_seq.SRR_master.configure import SEP
from Sc_RNA_seq.SRR_master.configure import VALUE
from Sc_RNA_seq.SRR_master.configure import THRESHOLD
from Sc_RNA_seq.SRR_master.configure import PRINT_DETAIL
from Sc_RNA_seq.SRR_master.configure import N_INTEGRATION


@log
def get_gse_dict(srr_file, sep=SEP):
    """读取gse：srr的字典"""

    print('读取SRR_INFO：%s' % srr_file)
    with open(srr_file, 'r', encoding='utf8') as f:
        srr_info = pd.read_csv(f, sep=sep)
    if srr_info.shape[0] == 0:
        print('%s：文件里什么信息都没有！' % srr_file)
        exit()
    gse_dict = {}
    srr_count = 0
    for group in srr_info.groupby('GSE'):
        gse = group[0]
        data = group[1]
        if len(set(data['organism'])) > 1:
            text = f'{gse}中的SRR属于多个物种，将不予处理！'
            print(add_color(text, 'yellow'))
        else:
            # 对每个GSE的SRR排下序
            srr_list = data['SRR'].tolist()
            srr_list.sort()
            gse_dict[gse] = srr_list
            srr_count += len(data['SRR'])
    print('读取到%d个GSE的%d个SRR信息' % (len(gse_dict), srr_count))
    if len(gse_dict) == 0:
        exit()

    return gse_dict


@log
def srr_pool(srr_data):
    """将所有result文件下以SRR开头的文件移动到all_result_path中"""

    result_path = os.path.join(srr_data, 'result')
    all_result_path = os.path.join(result_path, 'all_result')
    make_dir(all_result_path)

    print('正在移动文件')
    for dir_path, _, _ in os.walk(result_path):
        abs_all_result_path = os.path.abspath(all_result_path)
        abs_dir_path = os.path.abspath(dir_path)
        if abs_dir_path != abs_all_result_path:
            file_list = [srr for srr in os.listdir(dir_path) if srr.startswith('SRR')]
            if len(file_list) > 0:
                print(f'{abs_dir_path}/*\t>\t{abs_all_result_path}/')
                command = f'mv {abs_dir_path}/* {abs_all_result_path}/'
                os.system(command)

    if len(os.listdir(all_result_path)) == 0:
        text = f'{result_path}中没有数据！'
        print(add_color(text, 'red'))
        exit()


@log
def integrate_srr(gse_dict, srr_data, gse_data, print_detail=PRINT_DETAIL):
    """按GSE整合SRR"""

    result_path = os.path.join(srr_data, 'result')
    all_result_path = os.path.join(result_path, 'all_result')
    all_result_srr = set(os.listdir(all_result_path))

    for gse_count, (gse, srr_list) in enumerate(gse_dict.items(), 1):
        print(f'===================={gse_count}/{len(gse_dict)}====================')
        finished_srr = [srr for srr in srr_list if srr in all_result_srr]

        # ============输出处理报告============
        completion = f'({len(finished_srr)}/{len(srr_list)})'
        if len(finished_srr) == 0:
            text = f'{gse}中的SRR完全没有被处理{completion}！'
            print(add_color(text, 'red'))
        else:
            # 输出每个GSE中SRR处理信息
            error_srr = set(srr_list) - set(finished_srr)
            if len(error_srr) == 0:
                text = f'{gse}中的SRR处理完全{completion}！'
                print(add_color(text, 'green'))
            else:
                text = f'{gse}中的SRR缺失{completion}！'
                print(add_color(text, 'yellow'))
                if print_detail:
                    print(add_color(error_srr, 'yellow'))

            # ============对完成度大于阈值的GSE进行整合============
            # 一个GSE中的SRR处理完成度大于这个阈值才会整合成表达矩阵
            if len(finished_srr) / len(srr_list) < THRESHOLD:
                text = f'{gse}中SRR缺失过多，不予整合！'
                print(add_color(text, 'yellow'))
            else:
                gse_dir = os.path.join(gse_data, gse)
                make_dir(gse_dir)
                matrix_file = os.path.join(gse_dir, 'matrix.csv')

                # 如果当前GSE已经存在一个matrix.csv文件
                # 会检查处理好的srr是否已经写好了，从而避免多次处理
                if os.path.exists(matrix_file):
                    with open(matrix_file, 'r', encoding='utf8') as f:
                        matrix_srr = [srr for srr in f.readline().strip().split(',') if srr]
                        if set(matrix_srr) == set(finished_srr):
                            print('已找到整合好的matrix.csv文件！')
                            continue

                # 首先抽一个文件读取gene的列表
                # 经过验证，每个featureCounts出来的文件基因的次序和个数是相同的
                with open(os.path.join(all_result_path, finished_srr[0]), 'r', encoding='utf8') as f:
                    # 读掉开始两行
                    f.readline()
                    f.readline()
                    genes_length = []
                    genes_list = []
                    for line in f:
                        genes_length.append(int(line.strip().split('\t')[-2]))
                        genes_list.append(line.strip().split('\t')[0])

                # 多进程整合SRR成GSE的表达矩阵
                print('开启%d个进程整合SRR数据！' % N_INTEGRATION)
                n_worker, srr_per_worker = distribute_srr(finished_srr, n_worker=N_INTEGRATION)
                pool = Pool(processes=n_worker)
                new_integration_worker = partial(integration_worker, all_result_path=all_result_path)
                result = pool.map(new_integration_worker, srr_per_worker)
                gse_data_dict = {key: every_dict[key] for every_dict in result for key in every_dict}
                gse_matrix = pd.DataFrame(gse_data_dict, index=genes_list)  # 创建一个DataFrame写入表达矩阵
                if VALUE == 'RPKM':
                    cells_numi = gse_matrix.sum(axis=0)
                    gse_matrix = gse_matrix.div(cells_numi, axis=1).div(genes_length, axis=0) * 10**9
                elif VALUE == 'TPM':
                    foo = gse_matrix.div(genes_length, axis=0) * 1000
                    foo_numi = foo.sum(axis=0)
                    gse_matrix = foo.div(foo_numi) * 10**6
                print('整合完毕，保存数据中！')

                with open(matrix_file, 'w', encoding='utf8') as f:
                    # TPM、RPKM的数量级最小大概就是三位小数，所以文件保存保留三位小数
                    gse_matrix.to_csv(f, sep=',', header=True, index=True, float_format='%.3f')
                    text = '保存成功：%s' % os.path.abspath(matrix_file)
                    print(add_color(text, 'green'))


def main(srr_file, srr_data, gse_data):
    """主函数"""

    gse_dict = get_gse_dict(srr_file)  # INPUT：从SRR_INFO读取信息
    srr_pool(srr_data)  # 移动所有SRR到一个all_result目录下
    # <PROCESS：按GSE整合SRR><OUTPUT：GSE_DATA>
    integrate_srr(gse_dict, srr_data, gse_data)


if __name__ == '__main__':
    options, _ = getopt.getopt(sys.argv[1:], '', ['file=', 'data='])
    options_dict = dict(options)
    miss_key = {'--file', '--data'} - options_dict.keys()
    if miss_key:
        print('缺失了参数：%s' % str(miss_key))
        exit()

    SRR_INFO = options_dict['--file']  # INPUT: SRR_INFO
    SRR_DATA = options_dict['--data']  # INPUT: SRR_DATA
    OUTPUT = './GSE_data'  # OUTPUT: GSE_DATA
    main(SRR_INFO, SRR_DATA, OUTPUT)

__doc__ = """
            INPUT：GSE_INFO，GSE_DATA
            PROCESS：统计GSE处理结果
            OUTPUT：GSE_SUMMARY
            """

import os
import sys
import getopt
import pandas as pd

from Sc_RNA_seq.helper import log
from Sc_RNA_seq.helper import add_color

# 读取配置文件
from Sc_RNA_seq.GSE_master.configure import SEP


@log
def get_gse_organism_dict(gse_info, sep=SEP):
    """从GSE_info读取数据"""

    print('读取GSE_info:', gse_info)
    with open(gse_info, encoding='utf8') as f:
        gse_info_data = pd.read_csv(f, sep=sep)

    columns = gse_info_data.columns
    assert 'Organism' in columns
    assert 'Series' in columns
    gse_organism_dict = dict(zip(gse_info_data['Series'], gse_info_data['Organism']))
    print('得到%d个GSE的物种信息' % len(gse_organism_dict))
    if len(gse_organism_dict) == 0:
        exit()

    return gse_organism_dict


@log
def get_gse_summary(gse_organism_dict, gse_data):
    """统计GSE处理结果，行是基因，列是样本"""

    print('开始统计GSE处理结果')
    gse_summary_dict = {}
    for gse in os.listdir(gse_data):
        gse_dir = os.path.join(gse_data, gse)
        gse_file = os.path.join(gse_dir, 'matrix.csv')
        gse_pca_file = os.path.join(gse_dir, 'pca.csv')
        if os.path.isdir(gse_dir):
            # ==========判断GSE的物种==========
            if gse in gse_organism_dict:
                organism = gse_organism_dict[gse]
            else:
                organism = 'None'
                print(add_color(f'不存在{gse}的物种信息', 'yellow'))
            # ==========计算GSE表达矩阵的基因数目和细胞数目==========
            if os.path.isfile(gse_file):
                with open(gse_file, 'r', encoding='utf8') as f:
                    samples = len(f.readline().strip().split(',')) - 1
                    genes = 0
                    for _ in f:
                        genes += 1
            else:
                samples = genes = 'None'
                print(add_color(f'不存在{gse_file}', 'yellow'))
            # ==========计算聚类数目==========
            if os.path.isfile(gse_pca_file):
                with open(gse_pca_file, 'r', encoding='utf8') as f:
                    gse_pca_data = pd.read_csv(f, index_col=0)
                    clusters = len(set(gse_pca_data['cluster'][1:]))
            else:
                clusters = 'None'
                print(add_color(f'不存在{gse_pca_file}', 'yellow'))
            gse_summary_dict[gse] = {'organism': organism, 'samples': samples, 'genes': genes, 'clusters': clusters}

    print('统计到了%d个GSE处理的结果' % len(gse_summary_dict))
    gse_summary = pd.DataFrame(gse_summary_dict).reindex(['organism', 'samples', 'genes', 'clusters']).transpose()
    gse_summary.index.name = 'GSE'

    return gse_summary


@log
def save_gse_summary(gse_summary, gse_summary_file, sep=SEP):
    """保存GSE_summary文件"""

    print('正在保存：%s' % gse_summary_file)
    with open(gse_summary_file, 'w', encoding='utf8') as f:
        gse_summary.to_csv(f, sep=sep)
    print('保存成功')


def main(gse_info, gse_data, gse_summary_file):
    """主函数"""

    gse_organism_dict = get_gse_organism_dict(gse_info)  # INPUT：从GSE_INFO读取数据
    gse_summary = get_gse_summary(gse_organism_dict, gse_data)  # PROCESS：统计GSE处理结果
    save_gse_summary(gse_summary, gse_summary_file)  # OUTPUT：GSE_SUMMARY


if __name__ == '__main__':
    options, _ = getopt.getopt(sys.argv[1:], '', ['gseInfo=', 'gseData='])
    options_dict = dict(options)
    miss_key = {'--gseInfo', '--gseData'} - options_dict.keys()
    if miss_key:
        print('缺失了参数：%s' % str(miss_key))
        exit()

    GSE_INFO = options_dict['--gseInfo']  # INPUT: GSE_INFO
    GSE_DATA = options_dict['--gseData']  # INPUT: GSE_DATA
    OUTPUT = './GSE_summary.txt'  # OUTPUT: GSE_SUMMARY
    main(GSE_INFO, GSE_DATA, OUTPUT)

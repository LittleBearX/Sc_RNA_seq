__doc__ = """
            INPUT：GENE_INFO，GSE_INFO，GSE_DATA
            PROCESS：根据GENE_INFO分割数据为CODING_DATA和NCODING_data两个文件夹
            OUTPUT：CODING_DATA，NCODING_DATA
            """

import os
import sys
import getopt
import pandas as pd

from Sc_RNA_seq.helper import log
from Sc_RNA_seq.helper import make_dir
from Sc_RNA_seq.helper import add_color

# 读取配置文件
from Sc_RNA_seq.GSE_master.configure import SEP


@log
def get_organism_genes_dict(gene_info, sep=SEP):
    """从GENE_info读取数据"""

    print('读取GENE_info:', gene_info)
    with open(gene_info, encoding='utf8') as f:
        gene_data = pd.read_csv(f, sep=sep)

    columns = gene_data.columns
    assert 'Gene_stable_ID' in columns
    assert 'Organism' in columns
    assert 'Coding' in columns
    organism_genes_dict = {}
    for organism, data in gene_data.groupby('Organism'):
        foo = {'coding': set(data['Gene_stable_ID'][data['Coding'] == 'coding']),
               'ncoding': set(data['Gene_stable_ID'][data['Coding'] == 'ncoding'])}
        organism_genes_dict[organism] = foo
    print('读取到下面物种基因信息: \n%s' % set(organism_genes_dict.keys()))
    if len(organism_genes_dict) == 0:
        exit()

    return organism_genes_dict


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
    # 将空格替换为下划线
    gse_organism_dict = {series: organism.replace(' ', '_') for series, organism in gse_organism_dict.items()}
    print('得到%d个GSE的物种信息' % len(gse_organism_dict))
    if len(gse_organism_dict) == 0:
        exit()

    return gse_organism_dict


@log
def split_gse(gse_data, gse_organism_dict, organism_genes_dict, coding_data, ncoding_data):
    """切分表达矩阵获得编码和非编码两个文件"""

    print('切分表达矩阵获得编码和非编码两个文件')
    make_dir(coding_data)
    make_dir(ncoding_data)
    all_gse_data = os.listdir(gse_data)
    for count, gse in enumerate(all_gse_data, 1):
        print('========================================')
        gse_dir = os.path.join(gse_data, gse)
        gse_file = os.path.join(gse_dir, 'matrix.csv')
        if os.path.isdir(gse_dir) and not os.path.isfile(gse_file):
            # 存在文件夹，但文件夹里面没有matrix.csv会报错
            text = f'不存在{gse_file}'
            print(add_color(text, 'red'))
        else:
            if gse not in gse_organism_dict:
                text = f'GSE_info中没有{gse}的物种信息！'
                print(add_color(text, 'red'))
                continue

            organism = gse_organism_dict[gse].replace(' ', '_')
            if organism not in organism_genes_dict:
                text = f'{gse}: GENE_info中没有{organism}的基因信息！'
                print(add_color(text, 'red'))
                continue

            file_size = '%.3fM' % (os.path.getsize(gse_file)/(10**6))
            text = f'正在处理: {gse} {organism} {file_size} ({count}/{len(all_gse_data)})'
            print(add_color(text, 'yellow'))
            coding = organism_genes_dict[organism]['coding']
            ncoding = organism_genes_dict[organism]['ncoding']
            with open(gse_file) as f:
                matrix_data = pd.read_csv(f, index_col=0)

            # 判断表达矩阵行名中哪些是编码基因，哪些是非编码基因
            coding_genes = [gene for gene in matrix_data.index if gene in coding]
            ncoding_genes = [gene for gene in matrix_data.index if gene in ncoding]
            # 保存编码基因矩阵
            if coding_genes:
                print('找到%d个Coding genes' % len(coding_genes))
                coding_dir = os.path.join(coding_data, gse)
                coding_file = os.path.join(coding_dir, 'matrix.csv')
                make_dir(coding_dir)
                with open(coding_file, 'w') as f:
                    foo = matrix_data.loc[coding_genes, :]
                    foo.to_csv(f, sep=',')
            else:
                text = f'{gse_file}: 未发现Coding genes'
                print(add_color(text, 'yellow'))

            # 保存非编码基因矩阵
            if ncoding_genes:
                print('找到%d个Non coding genes' % len(ncoding_genes))
                ncoding_dir = os.path.join(ncoding_data, gse)
                ncoding_file = os.path.join(ncoding_dir, 'matrix.csv')
                make_dir(ncoding_dir)
                with open(ncoding_file, 'w') as f:
                    foo = matrix_data.loc[ncoding_genes, :]
                    foo.to_csv(f, sep=',')
            else:
                text = f'{gse_file}: 未发现Non coding genes'
                print(add_color(text, 'yellow'))

            text = f'处理完毕: {gse}'
            print(add_color(text, 'green'))


def main(gene_info, gse_info, gse_data, coding_data, ncoding_data):
    """主函数"""

    organism_genes_dict = get_organism_genes_dict(gene_info)  # INPUT：从GENE_INFO读取数据
    gse_organism_dict = get_gse_organism_dict(gse_info)  # INPUT：从GSE_INFO读取数据
    # <INPUT：从GSE_DATA> <OUTPUT：CODING_DATA> <OUTPUT：NCODING_DATA>
    split_gse(gse_data, gse_organism_dict, organism_genes_dict, coding_data, ncoding_data)


if __name__ == '__main__':
    options, _ = getopt.getopt(sys.argv[1:], '', ['geneInfo=', 'gseInfo=', 'gseData='])
    options_dict = dict(options)
    miss_key = {'--geneInfo', '--gseInfo', '--gseData'} - options_dict.keys()
    if miss_key:
        print('缺失了参数：%s' % str(miss_key))
        exit()

    GENE_INFO = options_dict['--geneInfo']  # INPUT: GENE_INFO
    GSE_INFO = options_dict['--gseInfo']  # INPUT: GSE_INFO
    GSE_DATA = options_dict['--gseData']  # INPUT: GSE_DATA
    CODING_DATA = './CODING_data'  # OUTPUT: CODING_DATA
    NCODING_DATA = './NCODING_data'  # OUTPUT: NCODING_DATA
    main(GENE_INFO, GSE_INFO, GSE_DATA, CODING_DATA, NCODING_DATA)

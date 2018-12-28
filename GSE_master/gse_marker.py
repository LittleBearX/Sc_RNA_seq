__doc__ = """
            INPUT：CELL_MARKER，GSE_INFO，GSE_DATA
            PROCESS：利用Fisher exact text检验为聚类结果表明类别
            """

import os
import sys
import getopt
import pandas as pd
from scipy import stats

from Sc_RNA_seq.helper import log
from Sc_RNA_seq.helper import add_color

# 读取配置文件
from Sc_RNA_seq.GSE_master.configure import SEP


@log
def get_cell_marker_dict(cell_marker, sep):
    """获取含有cell_marker的数据，返回物种每种细胞类型的marker基因"""

    print('读取cell_marker:', cell_marker)
    with open(cell_marker, encoding='utf8') as f:
        cell_marker_data = pd.read_csv(f, sep=sep)

    columns = cell_marker_data.columns
    for col_name in ['speciesType', 'tissueType', 'cellName', 'ensemblID']:
        try:
            assert col_name in columns
        except AssertionError:
            text = 'cell_marker must have column: ' + col_name
            print(add_color(text, 'red'))
            exit()

    cell_marker_dict = {}
    for organism, data in cell_marker_data.groupby('speciesType'):
        organism_cell_dict = {}
        for _, cell in data.iterrows():
            cell_type = (cell['tissueType'], cell['cellName'])
            gene_set = set(cell['ensemblID'].strip().split(','))
            organism_cell_dict[cell_type] = gene_set
        cell_marker_dict[organism] = organism_cell_dict

    print('读取到下面物种的marker基因信息: \n%s' % set(cell_marker_dict.keys()))
    if len(cell_marker_dict) == 0:
        exit()

    # 计算每个物种所有marker gene
    sum_dict = {}
    for organism in cell_marker_dict:
        organism_all_gene_set = set()
        for cell_type, gene_set in cell_marker_dict[organism].items():
            organism_all_gene_set.update(gene_set)
        sum_dict[organism] = organism_all_gene_set
    cell_marker_dict['all'] = sum_dict

    return cell_marker_dict


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
def gse_marker_handle(gse_data, gse_organism_dict, cell_marker_dict, odds_ratio_threshold=2,
                      p_value_threshold=0.01, method='greater'):
    """利用Fisher exact text检验为聚类结果标明类别"""

    assert method in {'two-sided', 'less', 'greater'}
    all_gse_data = os.listdir(gse_data)
    for count, gse in enumerate(all_gse_data, 1):
        print('========================================')
        gse_dir = os.path.join(gse_data, gse)
        marker_genes_file = os.path.join(gse_dir, 'marker_genes.csv')
        if os.path.isdir(gse_dir) and not os.path.isfile(marker_genes_file):
            # 存在文件夹，但文件夹里面没有matrix.csv会报错
            text = f'不存在{marker_genes_file}'
            print(add_color(text, 'red'))
        else:
            if gse not in gse_organism_dict:
                text = f'GSE_info中没有{gse}的物种信息！'
                print(add_color(text, 'red'))
                continue

            organism = gse_organism_dict[gse].replace(' ', '_')
            if organism not in cell_marker_dict:
                text = f'{gse}: cell_marker中没有{organism}的marker基因信息！'
                print(add_color(text, 'red'))
                continue

            text = f'正在处理: {gse} {organism} ({count}/{len(all_gse_data)})'
            print(add_color(text, 'yellow'))
            with open(marker_genes_file, 'r', encoding='utf8') as f:
                marker_genes_data = pd.read_csv(f, sep=',')

            item_list = []
            all_marker = cell_marker_dict['all'][organism]  # 总marker
            for cluster, data in marker_genes_data.groupby('cluster'):
                cluster_marker = set(data['gene']) & all_marker  # 某个cluster的marker基因
                n_all_marker = len(all_marker)
                n_cluster_marker = len(cluster_marker)
                if n_cluster_marker == 0:
                    continue
                cluster_marker_prop = n_cluster_marker / n_all_marker  # cluster的marker占总marker的百分比
                for cell_type, cell_type_marker in cell_marker_dict[organism].items():
                    n_cell_type_marker = len(cell_type_marker)  # 某一cell_type的marker基因数目
                    # 随机情况下期望cluster的marker基因hit到cell_type的marker基因的数目
                    n_expected_hit = cluster_marker_prop * n_cell_type_marker
                    n_hit = len(cluster_marker & cell_type_marker)  # 实际hit到的数目
                    odds_ratio = n_hit / n_expected_hit  # 实际hit到的高于随机hit到的比例
                    if odds_ratio > odds_ratio_threshold:
                        # 构建列联表并进行Fisher exact test
                        n_non_hit_cell_type_marker = n_cell_type_marker - n_hit
                        n_non_hit_cluster_marker = n_cell_type_marker - n_hit
                        n_other_marker = n_all_marker - n_hit - n_non_hit_cell_type_marker - n_non_hit_cluster_marker
                        table = [[n_other_marker, n_non_hit_cell_type_marker], [n_non_hit_cluster_marker, n_hit]]
                        p_value = stats.fisher_exact(table, method)[1]
                        if p_value < p_value_threshold:
                            item = [cluster, n_all_marker, n_cluster_marker, n_cell_type_marker, n_hit,
                                    n_expected_hit, odds_ratio, p_value, organism, cell_type[0], cell_type[1]]
                            item_list.append(item)
            if item_list:
                item_data = pd.DataFrame(item_list)
                item_data.columns = ['cluster', 'n_all_marker', 'n_cluster_marker', 'n_cell_type_marker', 'n_hit',
                                     'n_expected_hit', 'odds_ratio', 'p_value', 'organism', 'tissueType', 'cellName']
                item_data.sort_values(by=['cluster', 'p_value'], inplace=True)
                cells_type_file = os.path.join(gse_dir, 'cells_type.csv')
                with open(cells_type_file, 'w', encoding='utf8') as f:
                    item_data.to_csv(f, index=False)
                text = f'处理完毕: {gse}'
                print(add_color(text, 'green'))
            else:
                text = f'没有cluster可以标记cell_type: {gse}'
                print(add_color(text, 'yellow'))

    text = '所有GSE都处理完毕！'
    print(add_color(text, 'green'))


def main(cell_marker, gse_info, gse_data):
    """主函数"""

    cell_marker_dict = get_cell_marker_dict(cell_marker, sep=SEP)  # INPUT：从cell_marker中读取数据
    gse_organism_dict = get_gse_organism_dict(gse_info)  # INPUT：从GSE_INFO读取数据
    gse_marker_handle(gse_data, gse_organism_dict, cell_marker_dict)  # PROCESS：利用Fisher exact text检验为聚类结果表明类别


if __name__ == '__main__':
    options, _ = getopt.getopt(sys.argv[1:], '', ['cellMarker=', 'gseInfo=', 'gseData='])
    options_dict = dict(options)
    miss_key = {'--cellMarker', '--gseInfo', '--gseData'} - options_dict.keys()
    if miss_key:
        print('缺失了参数：%s' % str(miss_key))
        exit()

    CELL_MARKER = options_dict['--cellMarker']  # INPUT: CELL_MARKER
    GSE_INFO = options_dict['--gseInfo']  # INPUT: GSE_INFO
    GSE_DATA = options_dict['--gseData']  # INPUT: GSE_DATA
    main(CELL_MARKER, GSE_INFO, GSE_DATA)

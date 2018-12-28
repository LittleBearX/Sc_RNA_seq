__doc__ = """
            INPUT：GSE_INFO
            PROCESS：利用API来获取pubmed的信息
            OUTPUT：PUBMED_INFO
            """

import re
import sys
import pandas as pd
from Sc_RNA_seq.helper import log
from Sc_RNA_seq.GSE_GSM_SRR.crawler import get_urls_dict
from Sc_RNA_seq.GSE_GSM_SRR.get_pubmed_info_helper import get_pubmed_info

# 读取配置文件
from Sc_RNA_seq.GSE_GSM_SRR.configure import SEP
from Sc_RNA_seq.GSE_GSM_SRR.configure import SLEEP
from Sc_RNA_seq.GSE_GSM_SRR.configure import METHOD
from Sc_RNA_seq.GSE_GSM_SRR.configure import N_REQUEST
from Sc_RNA_seq.GSE_GSM_SRR.configure import N_CRAWLER
from Sc_RNA_seq.GSE_GSM_SRR.configure import PRINT_DETAIL
from Sc_RNA_seq.GSE_GSM_SRR.configure import PUBMED_PREFIX
from Sc_RNA_seq.GSE_GSM_SRR.configure import COLLECT_PUBMED_INFO


@log
def get_pubmed_list(file_name, sep=SEP):
    """读取GSE文件获取pubmed ID信息"""

    print('读取GSE_info：%s' % file_name)
    with open(file_name, 'r', encoding='utf8') as f:
        gse_data = pd.read_csv(f, sep=sep, index_col=0)

    re_obj = re.compile('^[0-9]{8}$')  # PMID的合法格式
    pubmed_list = []
    for i in gse_data['Citation']:
        sub_pubmed_list = eval(i)
        for pubmed in sub_pubmed_list:
            if re_obj.match(pubmed):
                pubmed_list.append(pubmed)
    pubmed_list = list(set(pubmed_list))
    print('读取到%d个pubmed ID' % len(pubmed_list))
    if len(pubmed_list) == 0:
        exit()

    return pubmed_list


@log
def get_pubmed_data(pubmed_list):
    """利用API来获取pubmed的信息"""

    print('开始获取pubmed信息')
    if N_REQUEST > 100:
        print('API一次请求的网页过多了！')
        print('可修改：（configure/N_REQUEST)')
    split = range(0, len(pubmed_list), N_REQUEST)
    urls = []
    for start in split:
        end = start + N_REQUEST
        suffix = ','.join(pubmed_list[start: end])
        url = PUBMED_PREFIX + suffix
        urls.append(url)

    urls_dict = get_urls_dict(urls, get_pubmed_info, method=METHOD, n_crawler=N_CRAWLER,
                              print_detail=PRINT_DETAIL, sleep=SLEEP)

    # 整合结果并排好顺序
    pubmed_dict = {key: url_dict[key] for url_dict in urls_dict.values() for key in url_dict}
    pubmed_data = pd.DataFrame(pubmed_dict)
    pubmed_data = pubmed_data.transpose().reindex(columns=COLLECT_PUBMED_INFO)
    print('搜集到了%d个pubmed的信息' % len(pubmed_data['PMID']))

    # 输出错误报告
    error = set(pubmed_list) - set(pubmed_data['PMID'])
    if error:
        print('%d个pubmed信息没有被正确处理' % len(error))
        print(error)
    else:
        print('所有pubmed的信息都正确处理')

    return pubmed_data


@log
def save_pubmed_data(pubmed_data, file_name, sep=SEP):
    """将pubmed_data保存"""

    print('保存数据：%s' % file_name)
    with open(file_name, 'w', encoding='utf8') as f:
        pubmed_data.to_csv(f, sep=sep, index=False)
    print('已保存')


def main(input_file, output_file):
    """主函数"""

    pubmed_list = get_pubmed_list(input_file)  # INPUT：GSE_INFO
    pubmed_data = get_pubmed_data(pubmed_list)  # PROCESS：PROCESS：利用API来获取pubmed的信息
    save_pubmed_data(pubmed_data, output_file)  # OUTPUT：PUBMED_INFO


if __name__ == '__main__':
    # INPUT and OUTPUT
    INPUT = sys.argv[1]
    OUTPUT = './PUBMED_info.txt'
    main(INPUT, OUTPUT)

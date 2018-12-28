__doc__ = """
            INPUT：GSE
            PROCESS：处理GSE的html文件得到GSE的主要信息
            OUTPUT：GSE_INFO
            """

import re
import sys
from Sc_RNA_seq.helper import log
from Sc_RNA_seq.GSE_GSM_SRR.crawler import get_urls_dict
from Sc_RNA_seq.GSE_GSM_SRR.get_gse_info_helper import get_gse_info

# 读取配置文件
from Sc_RNA_seq.GSE_GSM_SRR.configure import SEP
from Sc_RNA_seq.GSE_GSM_SRR.configure import SLEEP
from Sc_RNA_seq.GSE_GSM_SRR.configure import METHOD
from Sc_RNA_seq.GSE_GSM_SRR.configure import N_CRAWLER
from Sc_RNA_seq.GSE_GSM_SRR.configure import GSE_PREFIX
from Sc_RNA_seq.GSE_GSM_SRR.configure import PRINT_DETAIL
from Sc_RNA_seq.GSE_GSM_SRR.configure import COLLECT_GSE_INFO


@log
def get_gse_list(file_name):
    """从文件中读取GSE"""

    print('读取GSE：%s' % file_name)
    re_obj = re.compile('GSE[0-9]+')  # 合法的GSE格式
    gse_list = set()
    with open(file_name, 'r', encoding='utf8') as f:
        for line in f:
            foo = re_obj.search(line.strip())
            if foo:
                gse_list.add(foo.group())

    gse_list = list(gse_list)
    print('读取到合法GSE格式%s个' % len(gse_list))
    if len(gse_list) == 0:
        exit()
    return gse_list


@log
def get_gse_dict(gse_list):
    """爬取gse_list的html信息"""

    print('开始获取GSE信息')
    urls = [GSE_PREFIX + i for i in gse_list]
    gse_dict = get_urls_dict(urls, func=get_gse_info, method=METHOD,
                             n_crawler=N_CRAWLER, print_detail=PRINT_DETAIL, sleep=SLEEP)
    print('%d个gse的html被爬取' % len(gse_dict))

    # 输出错误报告
    re_obj = re.compile('GSE[0-9]+')
    gse_list2 = [re_obj.search(i).group() for i in gse_dict.keys()]
    error = set(gse_list) - set(gse_list2)
    if error:
        print('%d个GSE没有被正确爬取' % len(error))
        print(error)
    else:
        print('所有GSE都被正确爬取')
    return gse_dict


@log
def save_gse_dict(gse_dict, file_name, sep=SEP):
    """将gse_dict保存"""

    print('保存数据：%s' % file_name)
    with open(file_name, 'w', encoding='utf8') as f:
        all_data = [sep.join(COLLECT_GSE_INFO)]
        for gse in gse_dict:
            data = sep.join([str(gse_dict[gse][item]) for item in COLLECT_GSE_INFO])
            all_data.append(data)
        f.write('\n'.join(all_data))
    print('已保存')


def main(input_file, output_file):
    """主函数"""

    gse_list = get_gse_list(input_file)  # INPUT：GSE
    gse_dict = get_gse_dict(gse_list)  # PROCESS：爬取并处理gse的html得到gse_dict
    save_gse_dict(gse_dict, output_file)  # OUTPUT：GSE_INFO


if __name__ == '__main__':
    # INPUT and OUTPUT
    INPUT = sys.argv[1]
    OUTPUT = './GSE_info.txt'
    main(INPUT, OUTPUT)

__doc__ = """
            INPUT：GSM_INFO
            PROCESS：利用API来获取SRR的信息
            OUTPUT：SRR_INFO
            """

import re
import sys
import pandas as pd
from Sc_RNA_seq.helper import log
from Sc_RNA_seq.GSE_GSM_SRR.crawler import get_urls_dict
from Sc_RNA_seq.GSE_GSM_SRR.get_srr_info_helper import get_srr_info

# 读取配置文件
from Sc_RNA_seq.GSE_GSM_SRR.configure import SEP
from Sc_RNA_seq.GSE_GSM_SRR.configure import SLEEP
from Sc_RNA_seq.GSE_GSM_SRR.configure import METHOD
from Sc_RNA_seq.GSE_GSM_SRR.configure import N_REQUEST
from Sc_RNA_seq.GSE_GSM_SRR.configure import N_CRAWLER
from Sc_RNA_seq.GSE_GSM_SRR.configure import PRINT_DETAIL
from Sc_RNA_seq.GSE_GSM_SRR.configure import REQUEST_PREFIX
from Sc_RNA_seq.GSE_GSM_SRR.configure import COLLECT_SRR_INFO


@log
def get_srx_list(file_name):
    """读取GSM文件获取SRX信息"""

    print('读取GSM_info：%s' % file_name)

    # 因为有些GSM（即SRX）属于多个GSE，所以为了避免出错，
    # 只这里提取SRX信息，其他信息后面从API里面解析出来
    re_obj = re.compile('SRX[0-9]+')  # 合法的SRX格式
    srx_list = set()
    with open(file_name, 'r', encoding='utf8') as f:
        for line in f:
            foo = re_obj.search(line.strip())
            if foo:
                srx_list.add(foo.group())

    srx_list = list(srx_list)
    print('读取到合法SRX格式%d个' % len(srx_list))
    if len(srx_list) == 0:
        exit()

    return srx_list


@log
def get_srr_data(srx_list):
    """利用API来获取SRR的信息"""

    print('开始获取SRR的信息')
    if N_REQUEST > 100:
        print('API一次请求的网页过多了！')
        print('可修改：（configure/N_REQUEST)')
    split = range(0, len(srx_list), N_REQUEST)
    urls = []
    for start in split:
        end = start + N_REQUEST
        suffix = ','.join(srx_list[start: end])
        url = REQUEST_PREFIX + suffix
        urls.append(url)
    urls_dict = get_urls_dict(urls, get_srr_info, method=METHOD, n_crawler=N_CRAWLER,
                              print_detail=PRINT_DETAIL, sleep=SLEEP)

    # 整合结果并排好顺序
    srr_dict = {key: url_dict[key] for url_dict in urls_dict.values() for key in url_dict}
    srr_data = pd.DataFrame(srr_dict)
    srr_data = srr_data.transpose().reindex(columns=COLLECT_SRR_INFO)
    srr_data.sort_values('GSM', inplace=True)
    srr_data.sort_values('GSE', inplace=True)
    print('搜集到了%d个SRR的信息' % len(srr_data['SRR']))

    # 输出错误报告
    error = set(srx_list) - set(srr_data['SRX'])
    if error:
        print('%d个SRX信息没有被正确处理' % len(error))
        print(error)
    else:
        print('所有SRX的信息都正确处理')

    return srr_data


@log
def save_srr_data(srr_data, file_name, sep=SEP):
    """将srr_data保存"""

    print('保存数据：%s' % file_name)
    with open(file_name, 'w', encoding='utf8') as f:
        srr_data.to_csv(f, sep=sep, index=False)
    print('已保存')


def main(input_file, output_file):
    """主函数"""

    srx_list = get_srx_list(input_file)  # INPUT：GSM_INFO
    srr_data = get_srr_data(srx_list)  # PROCESS：利用API来获取SRR的信息
    save_srr_data(srr_data, output_file)  # OUTPUT：SRR_INFO


if __name__ == '__main__':
    # INPUT and OUTPUT
    INPUT = sys.argv[1]
    OUTPUT = './SRR_info.txt'
    main(INPUT, OUTPUT)

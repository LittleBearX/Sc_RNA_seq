__doc__ = """
            # INPUT：GSE_INFO
            # PROCESS：解析Series matrix文件得到GSM的信息
            # OUTPUT：GSM_INFO
            """

import sys
import pandas as pd
from Sc_RNA_seq.helper import log
from Sc_RNA_seq.GSE_GSM_SRR.get_gsm_info_helper import get_all_gsm_data

# 读取配置文件
from Sc_RNA_seq.GSE_GSM_SRR.configure import SEP
from Sc_RNA_seq.GSE_GSM_SRR.configure import PRINT_DETAIL


@log
def get_gse_data(file_name, sep=SEP):
    """读取GSE的信息"""

    print('读取GSE_info：%s' % file_name)
    with open(file_name, 'r', encoding='utf8') as f:
        gse_data = pd.read_csv(f, sep=sep, index_col=0)
    print('读取到%d个GSE' % gse_data.shape[0])
    if gse_data.shape[0] == 0:
        exit()
    return gse_data


@log
def get_gsm_data(gse_data, print_detail=PRINT_DETAIL):
    """从GSE的Series matrix文件和SRX的html得到GSM的信息"""

    print('开始处理%d个GSE的GSM信息' % gse_data.shape[0])

    # 从Series matrix得到GSM的主要信息
    print('从Series matrix得到GSM的主要信息')
    gsm_data = get_all_gsm_data(gse_data, print_detail)
    print('处理了%d个GSE的GSM信息' % len(set(gsm_data['Series'])))

    # 输出错误报告
    error = set(gse_data.index) - set(gsm_data['Series'])
    if error:
        print('%d个GSE的GSM信息没有被正确处理' % len(error))
        print(error)
    else:
        print('所有GSE的GSM信息都正确处理')

    return gsm_data


@log
def save_gsm_data(gsm_data, file_name, sep=SEP):
    """将gsm_data保存"""

    print('保存数据：%s' % file_name)
    with open(file_name, 'w', encoding='utf8') as f:
        gsm_data.to_csv(f, sep=sep)
    print('已保存')


def main(input_file, output_file):
    """主函数"""

    gse_data = get_gse_data(input_file)  # INPUT：GSE_INFO
    gsm_data = get_gsm_data(gse_data)  # PROCESS：解析Series matrix文件得到GSM的信息
    save_gsm_data(gsm_data, output_file)  # OUTPUT：GSM_INFO


if __name__ == '__main__':
    # INPUT and OUTPUT
    INPUT = sys.argv[1]
    OUTPUT = './GSM_info.txt'
    main(INPUT, OUTPUT)

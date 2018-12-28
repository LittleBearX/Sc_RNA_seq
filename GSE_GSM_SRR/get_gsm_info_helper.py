__doc__ = """此模块辅助获取GSM的信息"""

import os
import re
import pandas as pd
from ftplib import FTP
from Sc_RNA_seq.helper import unpack_gz

# 读取配置文件
from Sc_RNA_seq.GSE_GSM_SRR.configure import COLLECT_GSM_INFO


# ====================处理Series matrix的主要函数====================
def download_series_matrix(ftp):
    """通过ftp下载Series matrix数据"""

    ftp_server = 'ftp.ncbi.nlm.nih.gov'
    ftp_file_path = ftp[26:]
    ftp_client = FTP(ftp_server)
    ftp_client.login()
    ftp_client.cwd(ftp_file_path)
    ftp_files = ftp_client.nlst()
    for ftp_file in ftp_files:
        ftp_client.retrbinary('RETR ' + ftp_file, open(ftp_file, 'wb').write, blocksize=8192)
    return ftp_files


def unpack_series_matrix(ftp_files):
    """解压Series matrix数据，并删除原始文件"""

    series_matrix_data = []
    for ftp_file in ftp_files:
        data = unpack_gz(ftp_file)
        series_matrix_data.append(data)
        os.remove(ftp_file)
    return series_matrix_data


def get_one_series_matrix_data(data):
    """解析Series matrix文件的主要函数"""

    data_dict = {}
    relation = []
    lines = data.split('\n')
    collect_gsm_info = set(COLLECT_GSM_INFO)  # 集合的处理速度比列表快
    re_obj1 = re.compile('[\n\r]')
    re_obj2 = re.compile('(!Sample_)|(_ch1)|(_[0-9]+)')
    for line in lines:
        line = re_obj1.sub(' ', line.replace('"', ''))  # 替换掉奇怪的字符
        split_line = line.split('\t')
        line_key = split_line[0]
        line_value = split_line[1:]
        if '!Sample' in line_key:
            new_line_key = re_obj2.sub('', line_key)  # 将命名替换为标准格式
            # relation这条要特殊处理，分为BioSample和SRA处理
            # 有些奇葩文件BioSample和SRA交错写，所以下面代码会比较奇怪

            if new_line_key != 'relation':
                if new_line_key in collect_gsm_info:
                    # 有些信息会跨过多行，这么写可以把不同行信息拼起来
                    if new_line_key not in data_dict:
                        data_dict[new_line_key] = line_value
                    else:
                        sep = ';'
                        data_dict[new_line_key] = [i + sep + j for i, j in zip(data_dict[new_line_key], line_value)]
            else:
                relation.append(line_value)
    if relation:
        bio_sample = ['None'] * len(relation[0])
        sra = ['None'] * len(relation[0])
        for line_value in relation:
            for i, j in enumerate(line_value):
                if j.startswith('BioSample'):
                    bio_sample[i] = j.split(': ')[1]
                elif j.startswith('SRA'):
                    sra[i] = j.split(': ')[1]
        data_dict['BioSample'] = bio_sample
        data_dict['SRA'] = sra

    index = data_dict.pop('geo_accession')
    columns = COLLECT_GSM_INFO[:]
    columns.remove('geo_accession')
    data_frame = pd.DataFrame(data_dict, index=index, columns=columns)
    return data_frame


def get_series_matrix_data(gse_name, data_list):
    """整合一个GSE中所有Series matrix的文件的结果"""

    data_frame_list = [get_one_series_matrix_data(data) for data in data_list]
    data_frame = pd.concat(data_frame_list, axis=0)
    data_frame.fillna('None', inplace=True)
    new_columns = ['Series'] + data_frame.columns.tolist()
    data_frame = data_frame.reindex(columns=new_columns, fill_value=gse_name)
    return data_frame


def get_one_gsm_data(gse_name, gse_info):
    """下载，解压，解析一个GSE的所有Series matrix文件"""

    ftp = gse_info['Series matrix']
    ftp_files = download_series_matrix(ftp)
    data_list = unpack_series_matrix(ftp_files)
    data_frame = get_series_matrix_data(gse_name, data_list)
    return data_frame


def get_all_gsm_data(gse_data, print_detail):
    """下载，解压，解析所有GSE的所有Series matrix文件"""

    data_frame_list = []
    for gse_name, gse_info in gse_data.iterrows():
        # 因为网络中断，压缩包错误，ftp是None，ip被封等问题导致错误由try来处理
        try:
            data_frame = get_one_gsm_data(gse_name, gse_info)
            data_frame_list.append(data_frame)
            if print_detail:
                print('√ %s' % gse_name)
        except Exception:
            if print_detail:
                print('× %s' % gse_name)
    gsm_data = pd.concat(data_frame_list, axis=0)
    gsm_data.index.name = 'geo_accession'
    return gsm_data

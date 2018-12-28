__doc__ = """此模块辅助获取SRR的信息"""

import re
from bs4 import BeautifulSoup

# 读取配置文件
from Sc_RNA_seq.GSE_GSM_SRR.configure import COLLECT_SRR_INFO


def get_srr_info(html):
    """处理SRR的html得到srr_info"""

    if isinstance(html, bytes):
        html = html.decode('utf8')
    srr_info = get_srr_page_list(html)
    return srr_info


def get_key_value(srr_item, name, attrs, key, value):
    """辅助解析xml标签里面的信息，并以key:value形式返回"""

    foo_list = []
    for read in srr_item.findAll(name, attrs):
        _key = read.attrs[key]
        _value = read.attrs[value]
        foo_list.append('%s:%s' % (_key, _value))
    if foo_list:
        return ';'.join(foo_list)
    else:
        return 'None'


def get_attr(info, key, start, bad_choice):
    """尝试用key寻找info的信息，如果符合一定格式则返回"""

    if key in info and info[key].startswith(start):
        return info[key]
    else:
        # 如果在固定位置找不到符合一定格式的目标，下策是使用正则全文本匹配寻找
        re_obj = re.search(start + '[0-9]+', bad_choice)
        if re_obj:
            return re_obj.group()
        else:
            return 'None'


def get_srr_page_list(html):
    """处理html得到里面的SRR信息"""

    srr_page_dict = {}

    bs_obj = BeautifulSoup(html, 'xml')
    srx_items = bs_obj.findAll('EXPERIMENT_PACKAGE')
    for srx_item in srx_items:
        srx_item_dict = {i: 'None' for i in COLLECT_SRR_INFO}

        # 获取SRX的信息
        srx_name = srx_item.find('PRIMARY_ID').text
        sample_info = srx_item.find('SAMPLE').attrs
        study_info = srx_item.find('STUDY').attrs
        organism = srx_item.find('SCIENTIFIC_NAME').text
        bad_choice = str(srx_item.find('EXPERIMENT'))
        srx_item_dict['SRX'] = srx_name
        srx_item_dict['GSM'] = get_attr(sample_info, 'alias', 'GSM', bad_choice)
        srx_item_dict['SRS'] = get_attr(sample_info, 'accession', 'SRS', bad_choice)
        srx_item_dict['GSE'] = get_attr(study_info, 'alias', 'GSE', bad_choice)
        srx_item_dict['SRP'] = get_attr(study_info, 'accession', 'SRP', bad_choice)
        srx_item_dict['organism'] = organism

        # 获取SRX的每个RUN即SRR信息
        srr_items = srx_item.findAll('RUN')
        for srr_item in srr_items:
            srr_item_dict = srx_item_dict.copy()

            srr_name = srr_item.find('PRIMARY_ID').text
            srr_item_dict['SRR'] = srr_name
            srr_item_dict['read'] = get_key_value(srr_item, 'Read', {}, 'index', 'count')
            srr_item_dict['quality'] = get_key_value(srr_item, 'Quality', {}, 'value', 'count')
            srr_item_dict['bases'] = get_key_value(srr_item, 'Base', {}, 'value', 'count')
            srr_item_dict['taxon'] = get_key_value(srr_item, 'taxon', {'rank': 'species'}, 'name', 'total_count')

            srr_page_dict[srr_name] = srr_item_dict
    return srr_page_dict

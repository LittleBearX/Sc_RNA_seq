__doc__ = """此模块辅助获取GSE的html中信息"""

import re
import warnings
from bs4 import BeautifulSoup
from urllib.request import urlopen

# 读取配置文件
from Sc_RNA_seq.GSE_GSM_SRR.configure import PUBMED_PREFIX
from Sc_RNA_seq.GSE_GSM_SRR.configure import COLLECT_GSE_INFO


def get_gse_info(html):
    """处理GSE的html得到gse_info"""

    if isinstance(html, bytes):
        html = html.decode('utf8')
    gse_info = get_gse_page_dict(html)
    return gse_info


# ====================辅助函数====================
def children_list(bs_obj, tag=None):
    """找一个bs_obj的标签是tag的子标签"""

    children = [x for x in bs_obj.children if str(x).strip()]
    if tag:
        children = [child for child in children if child.name == tag]
    return children


def get_gse_page_dict_warn(gse_name, bs_obj, all_loc):
    """关于处理GSE的html页面的一些建议"""

    count = 0
    for loc in all_loc:
        count += len(bs_obj.findAll('table', loc))
    if count > 3:
        warnings.warn('警告：%s的页面结构是否异于常规GSE页面结构，解析的数据可能有误！' % gse_name)


# ====================处理页面的主要函数====================
def get_gse_page_dict(html):
    """处理GSE的html的主要函数"""

    gse_page_dict = {i: 'None' for i in COLLECT_GSE_INFO}
    _gse_page_dict = gse_page_dict.copy()
    bs_obj = BeautifulSoup(html, "html.parser")
    gse_name = bs_obj.find('strong', {'class': 'acc'}).attrs['id']  # 首先在页面中找到gse的名称
    _gse_page_dict['Series'] = gse_name

    # 大多GSE的html有三大块（有些GSE可能更多），main_content，download_family，supplementary_file，标签是table
    # 我这里用于在页面中定位这三大块的三个变量如下，但是GEO页面结构可能发生改变，如超过三个块的GSE，如GSE102066
    # 我在处理这三大块的时候尽量避免上面的情况，但对未知的情况很多，所以要检查输出文件看是否正确
    main_content_loc = {'cellpadding': '2', 'cellspacing': '0', 'width': '600'}
    download_family_loc = {'cellspacing': '3', 'width': '600'}
    supplementary_file_loc = {'cellpadding': '2', 'cellspacing': '2', 'width': '600'}
    all_loc = [main_content_loc, download_family_loc, supplementary_file_loc]
    get_gse_page_dict_warn(gse_name, bs_obj, all_loc)  # 对异常页面结构进行警告

    # 处理main_content
    main_content = bs_obj.find('table', main_content_loc)
    blocks = children_list(main_content, 'tr')
    for block in blocks:
        item = children_list(block, 'td')
        if len(item) == 2:
            item_key = item[0].text
            item_value = item[1].text
            item_value = item_value.replace('\t', ';')  # 将文本中的\t转换为空格，防止保存文件时出错
            item_value = item_value.replace('\n', ';')  # 将文本中的\n转换为空格，防止保存文件时出错
            # 对特殊的项目要单独写出来
            if item_key.startswith('Platforms'):
                # 较多GPL的话有些GPL会隐藏，比如GSE38495，所以这样写可以提出隐藏的GPL
                _gse_page_dict['Platforms'] = [GPL.text.split('\n')[:2] for GPL in item[1].findAll('tr')]
            elif item_key.startswith('Experiment type'):
                _gse_page_dict['Experiment type'] = [i for i in re.split('<..>|<...>', str(item[1])) if i]
            elif item_key.startswith('Samples'):
                # 提取GSE含有的Samples数目
                _gse_page_dict['Samples'] = re.search('\(([0-9]+)\)', item_key).group(1)
            elif item_key.startswith('Organism'):
                # 有些GSE的Organism会写成Organism(s)，下面Citation同理，所以这样写
                _gse_page_dict['Organism'] = item_value
            elif item_key.startswith('Citation'):
                re_obj = re.compile('^[0-9]{8}$')  # 合法PMID格式
                _gse_page_dict['Citation'] = [i.strip() for i in item_value.split(',') if re_obj.match(i.strip())]
            elif item_key in _gse_page_dict:
                _gse_page_dict[item_key] = item_value

    # 处理download_family
    try:
        download_family = bs_obj.find('table', download_family_loc)
        blocks = children_list(download_family, 'tr')[1:]  # 第一个是标题，不要
        for block in blocks:
            if block.text.startswith('Series Matrix'):
                _gse_page_dict['Series matrix'] = block.a.attrs['href']
    except Exception:
        pass

    # 处理supplementary_file
    try:
        supplementary_file = bs_obj.find('table', supplementary_file_loc)
        blocks = supplementary_file.findAll('tr', {'valign': "top"})[1:]  # 第一个是标题，不要
        https = []
        for block in blocks:
            file_name = block.find('td').text
            http_href = 'None'
            # 注意，要得到http文本下的超链接，但是目前发现有两种情况，如GSE100058和GSE38495，注意一下
            for i in block.findAll('a'):
                if i.text == '(http)':
                    http_href = 'https://www.ncbi.nlm.nih.gov' + i.attrs['href']
            https.append([file_name, http_href])
        _gse_page_dict['File'] = https
    except Exception:
        pass

    # 最后从调用API获取PUBMED的abstract信息
    if _gse_page_dict['Citation'] != 'None':
        try:
            abstract_list = []
            for pubmed_id in _gse_page_dict['Citation'].split(','):
                pubmed_url = PUBMED_PREFIX + pubmed_id.strip()
                pubmed_html = urlopen(pubmed_url).read().decode('utf8')
                pubmed_bs_obj = BeautifulSoup(pubmed_html, 'xml')
                abstract_list.append([pubmed_id, pubmed_bs_obj.find('AbstractText').text])
            _gse_page_dict['Abstract'] = abstract_list
        except Exception:
            pass

    return _gse_page_dict

__doc__ = """这是search模块的辅助模块"""

import re
from math import ceil
from urllib.request import urlopen

# 读取配置文件
from Sc_RNA_seq.GSE_GSM_SRR.configure import N_SEARCH
from Sc_RNA_seq.GSE_GSM_SRR.configure import SEARCH_PREFIX


def get_search_result(html, n_per_page):
    """显示搜索到的结果"""

    n_result = int(re.search('<Count>([0-9]+)</Count>', html).group(1))
    n_page = int(ceil(n_result / n_per_page))
    print('搜索到%d个结果，%d个页面' % (n_result, n_page))
    return n_page


def get_page_gse_list(html):
    """获取搜索到GSE"""

    re_obj1 = re.compile('<Id>([0-9]+)</Id>')
    re_obj2 = re.compile('^2[0]+')
    # GPL、GSE、GSM的search id都是9位数，GSE是以20+开头，把容易判断失效的条件放前面可加速判断
    page_gse_list = []
    for i in re_obj1.findall(html):
        if i.startswith('20') and len(i) == 9:
            page_gse_list.append('GSE' + re_obj2.sub('', i, 1))

    return page_gse_list


def get_one_search_data(keywords, print_detail):
    """得到一个关键词的搜索结果"""

    keywords = keywords.replace(' ', '+')
    n_per_page = N_SEARCH
    if n_per_page > 50000:
        print('一次搜索的条目太多了')
        print('可修改：configure/N_SEARCH')
        exit()
    url = SEARCH_PREFIX + keywords
    html = urlopen(url).read().decode('utf8')
    n_page = get_search_result(html, n_per_page)

    n_finished = 0
    gse_list = []
    for page in range(n_page):
        if print_detail:
            print('正在搜索第%d个页面' % (page + 1))
        url = '%s%s&retstart=%d&retmax=%d' % (SEARCH_PREFIX, keywords, n_finished, n_per_page)
        html = urlopen(url).read().decode('utf8')
        page_gse_list = get_page_gse_list(html)
        gse_list.extend(page_gse_list)
        n_finished += n_per_page
    print('搜索结果经过过滤，搜索到%d个GSE' % (len(gse_list)))

    return gse_list

__doc__ = """此模块辅助获取SRR的信息"""

from bs4 import BeautifulSoup

# 读取配置文件
from Sc_RNA_seq.GSE_GSM_SRR.configure import COLLECT_PUBMED_INFO


def get_pubmed_info(html):
    """处理pubmed的html得到pubmed_info"""

    if isinstance(html, bytes):
        html = html.decode('utf8')
    pubmed_info = get_pubmed_page_list(html)
    return pubmed_info


def get_pubmed_page_list(html):
    """处理html得到里面的pubmed信息"""

    pubmed_page_dict = {}

    bs_obj = BeautifulSoup(html, 'xml')
    pubmed_items = bs_obj.findAll('PubmedArticle')
    for pubmed_item in pubmed_items:
        pubmed_item_dict = {i: 'None' for i in COLLECT_PUBMED_INFO}

        # 获取pubmed的信息
        pubmed_name = pubmed_item.find('PMID').text.strip()
        pubmed_item_dict['PMID'] = pubmed_name
        if pubmed_item.find('ELocationID', {'EIdType': 'doi'}):
            pubmed_item_dict['DOI'] = pubmed_item.find('ELocationID', {'EIdType': 'doi'}).text.strip()
        journal_info = pubmed_item.find('Journal')
        if journal_info.find('PubDate'):
            pubmed_item_dict['PubDate'] = journal_info.find('PubDate').text.replace('\n', ' ').strip()
        if journal_info.find('Title'):
            pubmed_item_dict['Title'] = journal_info.find('Title').text.strip()
        if pubmed_item.find('ArticleTitle'):
            pubmed_item_dict['ArticleTitle'] = pubmed_item.find('ArticleTitle').text.strip()
        if pubmed_item.find('AbstractText'):
            pubmed_item_dict['AbstractText'] = pubmed_item.find('AbstractText').text.strip()

        author_item_list = pubmed_item.findAll('Author')
        author_list = []
        for author_item in author_item_list:
            try:
                name = author_item.find('LastName').text.strip() + ' ' + author_item.find('ForeName').text.strip()
                if author_item.find('Affiliation'):
                    affiliation = author_item.find('Affiliation').text.replace('\t', ' ').replace('\n', ' ').strip()
                else:
                    affiliation = 'None'
                author_list.append([name, affiliation])
            except Exception:
                pass
        pubmed_item_dict['Author'] = author_list
        pubmed_page_dict[pubmed_name] = pubmed_item_dict

    return pubmed_page_dict

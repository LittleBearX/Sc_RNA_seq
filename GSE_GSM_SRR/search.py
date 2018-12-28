__doc__ = """
            INPUT：KEYWORDS
            PROCESS：利用API来获取GEO搜索结果
            OUTPUT：GSE
            """

import sys
from Sc_RNA_seq.helper import log
from Sc_RNA_seq.GSE_GSM_SRR.search_helper import get_one_search_data

# 读取配置文件
from Sc_RNA_seq.GSE_GSM_SRR.configure import PRINT_DETAIL


@log
def get_keywords_list(file_name):
    """从文件中读取keywords"""

    print('读取KEYWORDS：%s' % file_name)
    with open(file_name, 'r', encoding='utf8') as f:
        keywords_list = [i.strip().lower() for i in f if i.strip()]  # 把搜索关键词转换为小写
    # 去除重复的关键词，为了保持顺序不用set
    new_keywords_list = []
    for keywords in keywords_list:
        if keywords not in new_keywords_list:
            new_keywords_list.append(keywords)
    print('读取到%d个Keywords' % len(new_keywords_list))
    if len(new_keywords_list) == 0:
        exit()
    return new_keywords_list


@log
def get_gse_list(keywords_list, print_detail=PRINT_DETAIL):
    """得到所有关键词的搜索结果"""

    print('开始搜索')
    gse_list = []
    for i, keywords in enumerate(keywords_list, 1):
        print('正在搜索第%d个关键词：%s' % (i, keywords))
        gse_list.extend(get_one_search_data(keywords, print_detail))
        print('====================')
    gse_list = list(set(gse_list))  # 删除不同关键词搜到重复的结果
    print('一共搜索到%d个GSE' % (len(gse_list)))
    return gse_list


@log
def save_gse_list(gse_list, file_name):
    """将gse_list保存"""

    print('保存数据：%s' % file_name)
    with open(file_name, 'w', encoding='utf8') as f:
        write_data = '\n'.join(gse_list)
        f.write(write_data)
    print('已保存')


def main(input_file, output_file):
    """主函数"""

    keywords_list = get_keywords_list(input_file)  # INPUT：KEYWORDS
    gse_list = get_gse_list(keywords_list)  # PROCESS：利用API来获取GEO搜索结果
    save_gse_list(gse_list, output_file)  # OUTPUT：GSE


if __name__ == '__main__':
    # INPUT and OUTPUT
    INPUT = sys.argv[1]
    OUTPUT = './GSE.txt'
    main(INPUT, OUTPUT)

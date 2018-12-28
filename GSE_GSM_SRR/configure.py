__doc__ = """GSE_GSM_SRR的配置文件"""

# ====================常用参数====================
SEP = '\t'  # 文件保存时的分隔符号
PRINT_DETAIL = True  # 是否显示处理的详细过程
N_SEARCH = 10000  # 利用API每次搜索的条目个数，建议不要超过50000，容易请求错误
N_REQUEST = 10  # 利用API每次请求SRX的个数，建议不要超过100，容易请求错误
N_CRAWLER = 30  # 线程数、进程数或者协程数，越大越容易把带宽跑满，但操作系统调度负担大
SLEEP = 3  # 爬取一个网页后等待几秒，再爬下个网页，防止请求网页很多时被NCBI封ip
METHOD = 'multithreading'  # 串行（serial），多进程（multiprocessing），多线程（multithreading），异步IO（asyncio）
# 爬取速度 multiprocessing > multithreading > asyncio > serial，消耗资源顺序正好相反

# ====================以下参数一般保持默认就好了====================
GSE_PREFIX = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='  # GEO数据库中的GSE的url前缀(deprecated)
SEARCH_PREFIX = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term='  # GEO搜索API的前缀
PUBMED_PREFIX = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&rettype=abstract&id='  # PUBMED的前缀
REQUEST_PREFIX = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id='  # SRR请求API的前缀
# SEARCH_PREFIX = 'https://www.ncbi.nlm.nih.gov/gds/?term='  # GEO搜索界面的前缀(deprecated)

# 搜集GSE的主要信息
COLLECT_GSE_INFO = ['Series', 'Status', 'Title', 'Organism', 'Experiment type', 'Summary', 'Overall design',
                    'Contributor(s)', 'Citation', 'Submission date', 'Last update date', 'Contact name', 'E-mail',
                    'Phone', 'Organization name', 'Department', 'Lab', 'Street address', 'City', 'State/province',
                    'ZIP/Postal code', 'Country', 'Platforms', 'Samples', 'BioProject', 'SRA', 'Series matrix', 'File']

# 搜集GSM的主要信息
COLLECT_GSM_INFO = ['title', 'geo_accession', 'source_name', 'organism', 'taxid', 'characteristics', 'molecule',
                    'platform_id', 'instrument_model', 'extract_protocol', 'growth_protocol', 'data_processing',
                    'library_selection', 'library_source', 'library_strategy', 'BioSample', 'SRA',
                    'supplementary_file']

# 搜集SRR的主要信息
COLLECT_SRR_INFO = ['SRR', 'SRX', 'GSM', 'SRS', 'GSE', 'SRP', 'organism', 'read', 'quality', 'bases', 'taxon']

# 搜集PUBMED的主要信息
COLLECT_PUBMED_INFO = ['PMID', 'DOI', 'PubDate', 'Title', 'ArticleTitle', 'AbstractText', 'Author']

__doc__ = """爬虫的辅助模块"""

import time
from urllib.request import urlopen


def get_html(urls, func, print_detail, timeout=5, try_times=3, sleep=2):
    """
    爬取url，并处理html，返回处理后数据，可并行
    :parameter urls: 要处理的urls
    :parameter func: 处理html的函数
    :parameter print_detail: 是否打印处理详细信息
    :parameter timeout: 尝试连接时间，timeout时间内连接不上就断开当前url请求
    :parameter try_times: 尝试连接次数，try_times次数连接不上就放弃连接
    :parameter sleep: 爬取页面的间隔时间
    :return urls_dict: key是每个url，value是func处理url的html的结果
    """

    urls_dict = {}
    target = urls
    n_try = 0
    while target and (n_try < try_times):
        n_try += 1
        errors = []
        for url in target:
            time.sleep(sleep)
            try:
                request = urlopen(url, timeout=timeout)
                html = request.read()
                urls_dict[url] = func(html)
                if print_detail:
                    print('√ %s' % url)
            except Exception:
                errors.append(url)
        target = errors
    if target:
        print('爬取失败！', target)
    return urls_dict

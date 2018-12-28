__doc__ = """此脚本实现串行、多进程、多线程、异步IO爬取urls"""

import time
import asyncio
from math import ceil
from functools import partial
from Sc_RNA_seq.helper import add_color
from Sc_RNA_seq.GSE_GSM_SRR.crawler_helper import get_html


# ====================辅助函数====================
def distribute_urls(urls, n_crawler):
    """用于给多线程或多进程爬虫分发任务"""

    n_urls_per_crawler = int(ceil(len(urls) / n_crawler))
    split = range(0, len(urls), n_urls_per_crawler)
    urls_per_crawler = [urls[i: i + n_urls_per_crawler] for i in split]
    n_crawler = len(urls_per_crawler)
    return n_crawler, urls_per_crawler


def chunked_http_client(n_crawler):
    """定义一个协程"""

    import aiohttp
    semaphore = asyncio.Semaphore(n_crawler)

    @asyncio.coroutine
    def http_get(url):
        nonlocal semaphore
        with (yield from semaphore):
            response = yield from aiohttp.ClientSession().get(url)
            data = yield from response.content.read()
            yield from response.wait_for_close()
        return data

    return http_get


def get_urls_dict_warn(urls, method, n_crawler, sleep):
    """关于爬取数据的一些建议"""

    if n_crawler > 100:
        print(add_color('不建议n_crawler超过100！', 'red'))
        print(add_color('可修改：（configure/N_CRAWLER)', 'red'))
    if sleep < 0.2:
        print(add_color('爬取页面间隔时间太短，小心被封ip!', 'red'))
        print('可修改：（configure/SLEEP)')
    if len(urls) > 5000:
        print(add_color('有%d个网页待爬取，建议分批进行！' % len(urls), 'red'))
    if len(urls) > 1000 and method == 'serial':
        print(add_color('有%d个网页待爬取，建议使用并行或并发方法！' % len(urls), 'red'))
        print(add_color('可修改：（configure/METHOD)', 'red'))
    if method != 'serial' and not isinstance(urls, list):
        print(add_color('urls应该是个有序的数据结构，并支持切片操作！'), "red")


# ====================串行、多进程、多线程、异步IO爬虫函数====================
def serial_handle(urls, func, print_detail, sleep):
    """串行爬取urls"""

    urls_dict = get_html(urls, func, print_detail, sleep)
    return urls_dict


def multiprocessing_handle(urls, func, n_crawler, print_detail, sleep):
    """多进程爬取urls"""

    from multiprocessing import Pool
    n_crawler, urls_per_crawler = distribute_urls(urls, n_crawler)
    pool = Pool(processes=n_crawler)
    new_get_html = partial(get_html, func=func, print_detail=print_detail, sleep=sleep)
    result = pool.map(new_get_html, urls_per_crawler)
    urls_dict = {key: every_dict[key] for every_dict in result for key in every_dict}
    return urls_dict


def multithreading_handle(urls, func, n_crawler, print_detail, sleep):
    """多线程爬取urls"""

    from multiprocessing.dummy import Pool
    n_crawler, urls_per_crawler = distribute_urls(urls, n_crawler)
    pool = Pool(processes=n_crawler)
    new_get_html = partial(get_html, func=func, print_detail=print_detail, sleep=sleep)
    result = pool.map(new_get_html, urls_per_crawler)
    urls_dict = {key: every_dict[key] for every_dict in result for key in every_dict}
    return urls_dict


def asyncio_handle(urls, func, n_crawler, print_detail, sleep):
    """异步IO爬取urls"""

    http_client = chunked_http_client(n_crawler)
    tasks = [http_client(url) for url in urls]
    urls_dict = {}
    for i, future in enumerate(asyncio.as_completed(tasks)):
        html = yield from future
        urls_dict[urls[i]] = func(html)
        if print_detail:
            print('√ %s' % urls[i])
        time.sleep(sleep)
    return urls_dict


# ====================主函数====================
def get_urls_dict(urls, func, method='serial', n_crawler=20, print_detail=True, sleep=0):
    """
    串行、多进程、多线程、异步IO爬虫
    :parameter urls: 要处理的urls
    :parameter func: 处理html的函数
    :parameter method: 串行（serial），多进程（multiprocessing），多线程（multithreading），异步IO（asyncio）
    :parameter n_crawler: 在multiprocessing中是进程数，multithreading中是线程数，asyncio中是协程数
    :parameter print_detail: 是否打印处理详细信息
    :parameter sleep: 爬取页面的间隔时间
    """

    get_urls_dict_warn(urls, method=method, n_crawler=n_crawler, sleep=sleep)
    if method == 'serial':
        print('使用串行爬取%d个网页' % len(urls))
        urls_dict = serial_handle(urls, func=func, print_detail=print_detail, sleep=sleep)
    elif method == 'multiprocessing':
        print('使用%d个进程爬取%d个网页' % (n_crawler, len(urls)))
        urls_dict = multiprocessing_handle(urls, func=func, n_crawler=n_crawler, print_detail=print_detail, sleep=sleep)
    elif method == 'multithreading':
        print('使用%d个线程爬取%d个网页' % (n_crawler, len(urls)))
        urls_dict = multithreading_handle(urls, func=func, n_crawler=n_crawler, print_detail=print_detail, sleep=sleep)
    elif method == 'asyncio':
        print('使用%d个协程爬取%d个网页' % (n_crawler, len(urls)))
        loop = asyncio.get_event_loop()
        urls_dict = loop.run_until_complete(asyncio_handle(urls, func=func, n_crawler=3,
                                                           print_detail=print_detail, sleep=sleep))
    else:
        raise TypeError('method必须是{serial, multiprocessing, multithreading, asyncio}其中之一')
    return urls_dict

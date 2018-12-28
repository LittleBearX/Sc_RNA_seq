__doc__ = """
            INPUT：SRR_INFO
            PROCESS：处理GSE的html文件得到GSE的主要信息
            OUTPUT：SRR_DATA
            """

import os
import re
import sys
import time
from functools import partial
from multiprocessing import Queue
from multiprocessing import Process
from Sc_RNA_seq.helper import log
from Sc_RNA_seq.helper import make_dir
from Sc_RNA_seq.helper import add_color
from Sc_RNA_seq.SRR_master.srr_helper import check
from Sc_RNA_seq.SRR_master.srr_helper import producer
from Sc_RNA_seq.SRR_master.srr_helper import consumer
from Sc_RNA_seq.SRR_master.srr_helper import wget_srr

# 读取配置文件
from Sc_RNA_seq.SRR_master.configure import CHECK
from Sc_RNA_seq.SRR_master.configure import N_DOWNLOAD
from Sc_RNA_seq.SRR_master.configure import PRINT_DETAIL


def check_first():
    """检查操作系统，硬盘，网络"""

    foo = input('是否需要检查环境：（y/n）')
    if foo != 'y':
        return
    result = check('operator_system')
    input('\nproceed')
    if result is not None:
        if 'Windows' in result['operator_system']:
            print('====================')
            from Sc_RNA_seq.SRR_master.configure import S
            s = bytes([(i + 100) % 256 for i in list(S)]).decode('utf8').split('\n')
            for i in s:
                print(i)
                time.sleep(2)
            foo = input('>>>是否继续？(y/n)')
            if foo != 'y':
                exit()
    check('device')
    input('\nproceed')
    check('network')
    input('\nproceed')


@log
def get_srr_list(file_name):
    """从文件中读取SRR"""

    print('读取SRR：%s' % file_name)
    re_obj = re.compile('SRR[0-9]+')  # 合法的SRR格式
    srr_list = list()
    with open(file_name, 'r', encoding='utf8') as f:
        for line in f:
            foo = re_obj.search(line.strip())
            if foo:
                srr = foo.group()
                if srr not in srr_list:
                    srr_list.append(srr)

    srr_list = list(srr_list)
    print('读取到合法SRR格式%d个' % len(srr_list))
    if len(srr_list) == 0:
        exit()
    return srr_list


@log
def download_srr(srr_list, output_dir):
    """使用生产者消费者模型下载数据"""

    print('开始下载')
    make_dir(output_dir)  # 创建文件夹

    # 启动生产者
    queue_size = 100
    queue = Queue(queue_size)
    pro = Process(target=producer, args=(srr_list, queue, PRINT_DETAIL, queue_size))
    pro.start()

    try:
        # 启动消费者队列
        consumer_list = []
        new_consumer = partial(consumer, path=output_dir)
        for consumer_name in range(1, N_DOWNLOAD + 1):
            con = Process(target=new_consumer, args=(consumer_name, queue, PRINT_DETAIL, wget_srr))
            con.start()
            consumer_list.append(con)
        # 等待消费者执行完毕
        for con in consumer_list:
            con.join()
    except KeyboardInterrupt:
        print(add_color('如果是Windows下的话按两下Ctrl+C, 父进程和子进程全部马上结束', 'red'))
        print(add_color('如果是Linux下的话按两下Ctrl+C，等待子进程完成最后一个任务才会退出', 'red'))

    # 输出错误报告
    finished = {i.split('.')[0] for i in os.listdir(output_dir)}
    error = set(srr_list) - finished
    if error:
        print('%s个SRR没有下载' % add_color(len(error), 'red'))
        print(add_color(error, 'red'))
    else:
        print(add_color('所有SRR都被下载', 'green'))
    # 子进程必须要外部停止才能关掉
    print('按Ctrl+C退出')


def main(input_file, output_dir):
    """主函数"""

    if CHECK:
        check_first()  # 检测电脑的各项指标
    srr_list = get_srr_list(input_file)  # INPUT：SRR_INFO
    download_srr(srr_list, output_dir)  # OUTPUT：SRR_DATA


if __name__ == '__main__':
    # INPUT and OUTPUT
    INPUT = sys.argv[1]
    OUTPUT = './SRR_data'
    main(INPUT, OUTPUT)

__doc__ = """GSE_master的辅助函数"""

import os
from math import ceil
from Sc_RNA_seq.helper import add_color
from Sc_RNA_seq.GSE_master.gse_worker import Project


def distribute_gse(gse_list, n_worker):
    """用于给多进程分发任务"""

    n_gse_per_worker = int(ceil(len(gse_list) / n_worker))
    split = range(0, len(gse_list), n_gse_per_worker)
    gse_per_worker = [gse_list[i: i + n_gse_per_worker] for i in split]
    n_worker = len(gse_per_worker)
    return n_worker, gse_per_worker


def gse_handle(sub_gse_list, gse_data):
    """单个进程对传入的GSE列表进行处理"""

    pid = str(os.getpid())
    for count, gse in enumerate(sub_gse_list, 1):
        prefix = f'进程: {pid}\t任务: {count}/{len(sub_gse_list)}\t'
        gse_file = os.path.abspath(os.path.join(gse_data, gse, 'matrix.csv'))
        if not os.path.exists(gse_file):
            text = prefix + '文件不存在: ' + gse_file
            print(add_color(text, 'red'))

        else:
            file_size = '%.2fM' % (os.path.getsize(gse_file)/(10**6))
            text = prefix + '开始处理: ' + gse + '\t大小: ' + file_size
            print(add_color(text, 'yellow'))
            try:
                project = Project(gse_file)
                project.compute_all_steps()
                text = prefix + '处理完毕: ' + gse
                print(add_color(text, 'green'))
            except Exception as e:
                text = prefix + '处理错误: ' + gse
                print(add_color(text, 'red'))
                print(add_color(str(e), 'red'))
    text = f'进程: {pid}已结束所有任务！'
    print(add_color(text, 'green'))

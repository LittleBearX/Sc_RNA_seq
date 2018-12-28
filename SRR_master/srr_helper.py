__doc__ = """SRR_master的辅助函数"""

import os
import psutil
import platform
from math import ceil
from Sc_RNA_seq.helper import add_color

# 读取配置文件
from Sc_RNA_seq.SRR_master.configure import PRINT_DETAIL
from Sc_RNA_seq.SRR_master.configure import DOWNLOAD_PREFIX
from Sc_RNA_seq.SRR_master.configure import PRINT_WGET_DETAIL


def wget_srr(srr, print_detail, kwargs):
    """在命令行执行下载srr"""

    path = kwargs['path']
    # 如果当前目录已存在文件，则不用下载
    if os.path.exists(os.path.join(path, srr + '.sra')):
        return

    # 拼接形成下载srr的url
    url_split = [DOWNLOAD_PREFIX, srr[:6], srr, srr + '.sra']
    url = '/'.join(url_split)

    try:
        # -c是断点续传
        command_list = ['wget -c', str(url)]
        if not PRINT_WGET_DETAIL:
            command_list.append('-q')
        command_list.extend(['-P', str(path)])
        command = ' '.join(command_list)
        info = os.system(command)
        if print_detail:
            if info:
                text = add_color('× ' + url, 'red')
                # 删掉下载错误的文件
                remove_file(os.path.join(path, srr + '.sra'))
            else:
                text = add_color('√ ' + url, 'green')
            print(text)
    except KeyboardInterrupt:
        # 这条try-except只会在Windows下会执行
        # Windows下子进程先捕获到KeyboardInterrupt，然后传给父进程
        # 而Linux是父进程先捕获到KeyboardInterrupt，然后子进程变成孤儿进程继续运行
        # 我在这里用try-except专门处理Windows下的情况，在consumer函数中专门处理Linux的情况

        # 删掉中断下载时那个没下载完的文件
        remove_file(os.path.join(path, srr + '.sra'))
        raise KeyboardInterrupt


def remove_file(target):
    """对于下载错误、处理错误、中间产生的文件进行删除"""

    if os.path.exists(target):
        os.remove(target)


def check(target):
    """检测操作系统、cpu、内存、硬盘空间、网络"""

    try:
        print('====================')
        result = {}
        if target == 'operator_system':
            print('检查操作系统：')
            operator_system = platform.platform()
            text = add_color(operator_system, 'green')
            result['operator_system'] = operator_system
        elif target == 'cpu':
            print('检查CPU：')
            true_cpu = psutil.cpu_count(logical=False)
            logical_cpu = psutil.cpu_count(logical=True)
            text = '物理核数：%s  逻辑核数：%s' % (add_color(str(true_cpu), 'green'), add_color(str(logical_cpu), 'green'))
            result['true_cpu'] = true_cpu
            result['logical_cpu'] = logical_cpu
        elif target == 'memory':
            print('检查内存：')
            size = psutil.virtual_memory()
            free = round(size.free / 10 ** 9, 3)
            used = round(size.used / 10 ** 9, 3)
            text = '内存   free: %s  used: %s' % (add_color(str(free) + 'G', 'green'), add_color(str(used) + 'G', 'red'))
            result['used'] = used
            result['free'] = free
        elif target == 'device':
            print('检查硬盘：')
            print(add_color('在Linux下结果可能不准，在命令行中输入df -h查看硬盘', 'red'))
            all_devices = psutil.disk_partitions()
            text_list = []
            for device in all_devices:
                size = psutil.disk_usage(device.device)
                free = add_color(str(round(size.free / 10 ** 9, 3)) + 'G', 'green')
                used = add_color(str(round(size.used / 10 ** 9, 3)) + 'G', 'red')
                text_list.append('%s   free: %s  used: %s' % (device.device, free, used))
            text = '\n'.join(text_list)
        elif target == 'network':
            print('检查网络：')
            print(add_color('Linux下如果不能自行停止请按Ctrl+C', 'red'))
            url = 'www.baidu.com'
            connect = os.system('ping %s' % url)
            if connect:
                text = add_color('%s 连接失败' % url, 'red')
            else:
                text = add_color('%s 连接成功' % url, 'green')
            result['connect'] = connect
        else:
            text = "target must be in {operator_system, cpu, memory, device, network}"

        print(text)
        return result
    except Exception:
        text = '无法检查当前操作系统的%s' % target
        print(add_color(text, 'red'))

        return None


def producer(task, queue, print_detail, queue_size):
    """将任务放进队列，在进程间共享"""

    targets = task[:]
    n_task = len(targets)
    n_put = - queue_size
    n_rest = queue_size + len(task)

    while True:
        n_rest -= 1
        n_put += 1
        if targets:
            queue.put(targets.pop())

        else:
            # 队列里的任务分发完了，就向队列中发送终止信号
            queue.put('close')

        # 输出已分发任务个数和已完成的任务
        if 0 < n_put <= n_task and print_detail:
            text = add_color('分发第%d个任务，剩余%d个任务' % (n_put, n_rest), 'green')
            print(text)


def consumer(consumer_name, queue, print_detail, fun, **kwargs):
    """消耗队列里的任务"""

    father_id = os.getppid()

    # 当父进程被杀死后，getppid将会得到1，所以子进程的while循环也会被终止
    # 我在这里专门处理Linux下的情况，在wget_url函数中专门处理Windows的情况
    while os.getppid() == father_id:
        target = queue.get()
        if target != 'close':
            if print_detail:
                consumer_name = str(consumer_name)
                # text = add_color('消费者'+consumer_name+' get：'+target, 'yellow')
                # print(text)
                pass
            fun(target, print_detail, kwargs)
        else:
            # 如果接受到了终止信号，则停止进程
            break


def distribute_srr(finished_srr, n_worker):
    """用于给多进程分发任务"""

    n_srr_per_worker = int(ceil(len(finished_srr) / n_worker))
    split = range(0, len(finished_srr), n_srr_per_worker)
    srr_per_worker = [finished_srr[i: i + n_srr_per_worker] for i in split]
    n_worker = len(srr_per_worker)
    return n_worker, srr_per_worker


def integration_worker(sub_finished_srr, all_result_path):
    """对分配的SRR子任务进行读取"""

    pid = os.getpid()
    sub_gse_data_dict = {}
    for count, srr in enumerate(sub_finished_srr, 1):
        if count % 300 == 0 and PRINT_DETAIL:
            print(f'进程: {pid} 完成: {count}/{len(sub_finished_srr)}')
        try:
            with open(os.path.join(all_result_path, srr)) as f:
                # 读掉开始两行
                f.readline()
                f.readline()
                expression = [int(line.strip().split('\t')[-1]) for line in f]
                sub_gse_data_dict[srr] = expression
        except Exception as e:
            text = '处理错误: ' + srr
            print(add_color(text, 'red'))
            print(add_color(str(e), 'red'))
    return sub_gse_data_dict

__doc__ = """
            INPUT：SRR_INFO、SRR_DATA
            PROCESS：预先分配任务
            OUTPUT：生成可在天河二号或本地服务器上运行的脚本SRR_SCRIPT
            """

import os
import sys
import shutil
import getopt
import pandas as pd
from Sc_RNA_seq.helper import log
from Sc_RNA_seq.helper import make_dir
from Sc_RNA_seq.helper import add_color
from Sc_RNA_seq.SRR_master.srr_helper import check

# 读取配置文件
from Sc_RNA_seq.SRR_master.configure import SEP
from Sc_RNA_seq.SRR_master.configure import STAR
from Sc_RNA_seq.SRR_master.configure import CHECK
from Sc_RNA_seq.SRR_master.configure import GENOME
from Sc_RNA_seq.SRR_master.configure import TIANHE
from Sc_RNA_seq.SRR_master.configure import CACHE_LOAD
from Sc_RNA_seq.SRR_master.configure import PYTHON_PATH
from Sc_RNA_seq.SRR_master.configure import PRINT_DETAIL
from Sc_RNA_seq.SRR_master.configure import GENOME_INDEX
from Sc_RNA_seq.SRR_master.configure import FASTERQ_DUMP
from Sc_RNA_seq.SRR_master.configure import FEATURE_COUNTS
from Sc_RNA_seq.SRR_master.configure import CPU_PER_WORKER
from Sc_RNA_seq.SRR_master.configure import TASK_PER_WORKER
from Sc_RNA_seq.SRR_master.configure import TIANHE_NODE_NAME


@log
def check_fist():
    """检查软件和数据是否正确"""

    print('正在检查软件和数据路径是否有误')
    all_right = True
    # 检查数据和软件
    for name, value in zip(['GENOME', 'GENOME_INDEX', 'FASTERQ_DUMP', 'STAR', 'FEATURE_COUNTS'],
                           [GENOME, GENOME_INDEX, FASTERQ_DUMP, STAR, FEATURE_COUNTS]):
        if not os.path.exists(value):
            print(add_color('×\t' + name, 'red'), end='\t')
            print(add_color('可修改：configure/' + name, 'red'))
            all_right = False
        else:
            print(add_color('√\t' + name, 'green'))

    # 单独检查python及版本
    foo = os.popen(PYTHON_PATH + ' -V').read().strip()
    if foo == '' or foo.split()[1].startswith('2'):
        print(add_color('×\tPYTHON_PATH', 'red'), end='\t')
        print(add_color('可修改：configure/PYTHON_PATH', 'red'))
        all_right = False
    else:
        print(add_color('√\tPYTHON_PATH', 'green'))

    if not all_right:
        exit()


@log
def check_last(srr_data):
    """检查操作系统，CPU，内存"""

    foo = input('是否需要检查环境：（y/n）')
    if foo != 'y':
        return
    # 检查操作系统
    result = check('operator_system')
    if result is not None:
        if 'Windows' in result['operator_system']:
            print(add_color('Windows下将不能运行此程序生成的脚本，请切换的Linux环境下！', 'red'))
            foo = input('继续？（y/n）')
            if foo != 'y':
                exit()

    # 检查CPU
    result = check('cpu')
    if result is not None:
        n_work = len(os.listdir(os.path.join(srr_data, 'distribution')))
        if TIANHE:
            print('在天河二号上运行，预计需要%d个节点' % n_work)
            print('如果大于60节点，要用BIOJOB分区的节点进行计算')
        else:
            expect_cpu = CPU_PER_WORKER * n_work
            print('预计需要物理核数：%s' % str(expect_cpu))
            if result['true_cpu'] < expect_cpu:
                print('如果CPU核数不够任务将来回切换，开销很大')
                print('可修改：（configure/CPU_PER_WORKER)')
            input('\nproceed')

    # 检查内存
    result = check('memory')
    if result is not None:
        cache_load_organism = {i.replace(' ', '_') for i in CACHE_LOAD}
        if cache_load_organism:
            print('%s基因组STAR_Index保留在内存以加速处理' % str(cache_load_organism))
            print('每个哺乳动物的Index大约30G')
            if result['free'] < len(cache_load_organism) * 30 + 30:
                print('如果内存不够将使用虚拟内存，速度会较慢')
                print('可修改：（configure/CACHE_LOAD)')
        input('\nproceed')

    # 检查硬盘
    check('device')
    print('检查最大的几个SRR文件')
    data_path = os.path.join(srr_data, 'data')
    result_path = os.path.join(srr_data, 'result')
    data_size = {}
    finished_srr = {srr + '.sra' for walk in os.walk(result_path) for srr in walk[2] if srr.startswith('SRR')}
    all_srr = list({srr for srr in os.listdir(data_path) if srr.endswith('.sra')} - finished_srr)
    for srr in all_srr:
        data_size[srr] = os.path.getsize(os.path.join(data_path, srr)) / (10 ** 9)
    # 输出data中最大的十个文件（不计算那些已经在result里的了）
    maxsize_data = sorted(data_size.keys(), key=lambda x: data_size[x], reverse=True)[:10]
    print('======================')
    for srr in maxsize_data:
        print('%s ： %.2fG' % (srr, data_size[srr]))
    print('======================')
    print('临时文件fastq：srr转换为fastq格式后大小扩大约4倍')
    print('临时文件sam：fastq转换为sam格式后大小扩大1至3倍不等')
    print('最终文件feature：典型的人类的一个样本大概29M，小鼠为19M，其它物种较小')
    print('多个任务同时处理时将会有多个临时文件，请确保硬盘空间足够')
    input('\nproceed')


@log
def make_srr_dir(srr_data):
    """创建文件的结构"""

    # 先检查SRR_data目录下是否有文件，再检查SRR_data/data目录下是否有文件
    srr_list = [i for i in os.listdir(srr_data) if i.endswith('.sra')]
    if len(srr_list) == 0:
        data_path = os.path.join(srr_data, 'data')
        srr_list = [i for i in os.listdir(data_path) if i.endswith('.sra')]
        if len(srr_list) == 0:
            print('%s和%s：文件夹下一个SRR的文件都没有！' % (srr_data, data_path))
            exit()
    for dir_name in ['data', 'script', 'distribution', 'temp', 'result']:
        dir_path = os.path.join(srr_data, dir_name)
        # 如下三个文件如果存在，必须删了重建，否则里面已经存在的东西会影响处理
        # result当然不能删掉重建，不然前面的处理就白费了
        if dir_name in ['script', 'distribution', 'temp']:
            if os.path.exists(dir_path):
                os.system('rm -r ' + dir_path)
        print('创建文件夹：%s' % os.path.abspath(dir_path))
        make_dir(dir_path)


@log
def get_srr_dict(srr_file, sep=SEP):
    """读取srr：organism的字典"""

    print('读取SRR_INFO：%s' % srr_file)
    with open(srr_file, 'r') as f:
        srr_info = pd.read_csv(f, sep=sep)
    if srr_info.shape[0] == 0:
        print('%s：文件里什么信息都没有！' % srr_file)
        exit()
    srr_dict = dict(zip(srr_info['SRR'], srr_info['organism']))

    return srr_dict


@log
def distribute_srr(srr_dict, srr_data, print_detail=PRINT_DETAIL):
    """为worker分配工作"""

    data_path = os.path.join(srr_data, 'data')
    result_path = os.path.join(srr_data, 'result')

    # 只处理后缀是.sra的文件
    all_srr = [srr for srr in os.listdir(srr_data) if srr.endswith('.sra')]
    # 将文件移动到data目录下
    print('正在移动文件')
    for count, srr in enumerate(all_srr, 1):
        from_path = os.path.join(srr_data, srr)
        to_path = os.path.join(data_path, srr)
        shutil.move(from_path, to_path)
        if print_detail and count % 300 == 0:
            print('处理进度：%d/%d' % (count, len(all_srr)))

    # 对于在result中的，即已经处理好的srr不再处理来，先剔除掉，做差集
    finished_srr = {srr + '.sra' for walk in os.walk(result_path) for srr in walk[2] if srr.startswith('SRR')}
    all_srr = list({srr for srr in os.listdir(data_path) if srr.endswith('.sra')} - finished_srr)
    print('正在分配%d个SRR文件' % len(all_srr))
    file_name_count = {}
    file_list = {}
    error = []
    for count, srr in enumerate(all_srr, 1):

        if print_detail and count % 300 == 0:
            print('处理进度：%d/%d' % (count, len(all_srr)))

        srr_name = srr.split('.')[0]
        try:
            # 从srr_dict中得到SRR所属的物种，注意要将空格改为下划线
            organism = srr_dict[srr_name].replace(' ', '_')
        except KeyError:
            # 未找到所属物种的srr文件将不会经过后续处理
            error.append(srr_name)
            continue

        # 这段代码在给每个worker分配任务，每个worker处理的任务来自同一个物种
        if organism not in file_name_count:
            file_name_count[organism] = 1
            file_list[organism] = []
        if len(file_list[organism]) >= TASK_PER_WORKER:
            # 当任务计数大于task_per_worker时，分配给下一个worker
            # 将任务记录在distribution的文件里
            file_name = organism + '_' + str(file_name_count[organism])
            file_name = os.path.join(srr_data, 'distribution', file_name)
            with open(file_name, 'w') as f:
                f.write('\n'.join(file_list[organism]))

            file_name_count[organism] += 1  # 文件名计数加一
            file_list[organism] = []  # 列表清空
        file_list[organism].append(srr)

    # 将最后的不足task_per_worker个的任务分配给worker
    for organism, srr_list in file_list.items():
        if srr_list:
            file_name = organism + '_' + str(file_name_count[organism])
            file_name = os.path.join(srr_data, 'distribution', file_name)
            with open(file_name, 'w') as f:
                f.write('\n'.join(srr_list))

    if error:
        text = '以下SRR文件未找到所属物种，请将信息添加到SRR_info\n%s' % str(error)
        print(add_color(text, 'red'))
    text = '分配任务完毕，在%s' % os.path.abspath(os.path.join(srr_data, 'distribution'))
    print(add_color(text, 'green'))


@log
def generate_script(srr_data, output_file):
    """生成处理的脚本"""

    print('正在批量生成处理脚本')
    # 对每个worker生成一个脚本
    script_header = '#!/bin/bash\n'
    script_template = PYTHON_PATH + ' {worker_file} {task}'
    script_list = []
    script_path = os.path.abspath(os.path.join(srr_data, 'script'))
    distribution_path = os.path.abspath(os.path.join(srr_data, 'distribution'))
    if TIANHE:
        worker_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'srr_worker_tianhe.py')
    else:
        worker_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'srr_worker.py')
    for i in os.listdir(distribution_path):
        task = os.path.join(distribution_path, i)
        file_name = 'Run_' + i + '.sh'  # 一个sh脚本运行一个distribution里面的文件，即一个task
        file_path = os.path.join(script_path, file_name)
        with open(file_path, 'w') as f:
            f.write(script_header)
            script_body = script_template.replace('{worker_file}', worker_file).replace('{task}', task)
            f.write(script_body)
        os.system('chmod 755 ' + file_path)
        script_list.append(file_path)
    text = '所有脚本生成完毕：%s' % os.path.abspath(os.path.join(srr_data, 'script'))
    print(add_color(text, 'green'))

    # 生成一个总脚本
    with open(output_file, 'w') as f:
        f.write(script_header)
        for script in script_list:
            if TIANHE:
                command = 'yhbatch -N 1 -n 1 -p %s %s' % (TIANHE_NODE_NAME, script)
            else:
                command = script
            f.write(command + ' &\n')
    os.system('chmod 755 ' + output_file)
    text = '总脚本生成完毕：%s' % os.path.abspath(output_file)
    print(add_color(text, 'green'))


def main(srr_file, srr_data, output_file):
    """主函数"""

    check_fist()  # 检查软件和数据
    make_srr_dir(srr_data)  # 创建文件结构
    srr_dict = get_srr_dict(srr_file)  # INPUT：从SRR_INFO读取信息
    distribute_srr(srr_dict, srr_data)  # PROCESS：预先分配任务
    generate_script(srr_data, output_file)  # OUTPUT：生成可在天河二号或本地服务器上运行的脚本SRR_SCRIPT
    if CHECK:
        check_last(srr_data)  # 检查电脑各项指标


if __name__ == '__main__':
    options, _ = getopt.getopt(sys.argv[1:], '', ['file=', 'data='])
    options_dict = dict(options)
    miss_key = {'--file', '--data'} - options_dict.keys()
    if miss_key:
        print('缺失了参数：%s' % str(miss_key))
        exit()

    SRR_INFO = options_dict['--file']  # INPUT: SRR_INFO
    SRR_DATA = options_dict['--data']  # INPUT: SRR_DATA
    OUTPUT = './SRR_script.sh'  # OUTPUT: SRR_SCRIPT
    main(SRR_INFO, SRR_DATA, OUTPUT)

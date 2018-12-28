__doc__ = """
            INPUT：SRR_list
            PROCESS：在天河二号上，利用sratoolkit、STAR、featureCounts处理SRR文件
            """

import os
import sys
from Sc_RNA_seq.helper import make_dir
from Sc_RNA_seq.helper import add_color
from Sc_RNA_seq.SRR_master.srr_helper import remove_file

# 读取配置文件
from Sc_RNA_seq.SRR_master.configure import STAR
from Sc_RNA_seq.SRR_master.configure import TIANHE
from Sc_RNA_seq.SRR_master.configure import GENOME
from Sc_RNA_seq.SRR_master.configure import REMOVE
from Sc_RNA_seq.SRR_master.configure import CACHE_LOAD
from Sc_RNA_seq.SRR_master.configure import GENOME_INDEX
from Sc_RNA_seq.SRR_master.configure import FASTERQ_DUMP
from Sc_RNA_seq.SRR_master.configure import FEATURE_COUNTS
from Sc_RNA_seq.SRR_master.configure import CPU_PER_WORKER
from Sc_RNA_seq.SRR_master.configure import TIANHE_CPU_PER_WORKER


def first_work(task_file):
    """定义各种参数，并加进全局变量"""

    global task_list, data_path, temp_path, result_path
    global organism, organism_genome, organism_genome_gtf, organism_genome_index
    global fasterq_dump, star, feature_counts
    global cpu_use, cache_load, remove

    # ====================定位到各种文件夹的路径====================
    distribution_path, task_name = os.path.split(task_file)
    srr_data = os.path.dirname(distribution_path)
    data_path = os.path.join(srr_data, 'data')
    temp_path = os.path.join(srr_data, 'temp', task_name)
    result_path = os.path.join(srr_data, 'result', task_name)
    make_dir(temp_path)
    make_dir(result_path)
    finished = {i+'.sra' for walk in os.walk(os.path.join(srr_data, 'result')) for i in walk[2] if i.startswith('SRR')}
    task_list = [i.strip() for i in open(task_file) if i not in finished]  # 在task_list中删掉那些已经处理完的
    # ====================定位到基因组信息的文件====================
    organism = '_'.join(task_name.split('_')[:-1])
    organism_genome = os.path.join(GENOME, organism)
    organism_genome_index = os.path.join(GENOME_INDEX, organism)

    for file in os.listdir(organism_genome):
        if file.endswith('.gtf'):
            organism_genome_gtf = os.path.join(organism_genome, file)
            break
    else:
        print('missed gtf file of %s: %s' % (organism, organism_genome))
        exit()
    # ====================定位到软件的位置====================
    star = STAR
    fasterq_dump = FASTERQ_DUMP
    feature_counts = FEATURE_COUNTS
    # ====================定义其他参数====================
    if TIANHE:
        cpu_use = TIANHE_CPU_PER_WORKER
    else:
        cpu_use = CPU_PER_WORKER
    cache_load_organism = {i.replace(' ', '_') for i in CACHE_LOAD}
    if organism in cache_load_organism:
        cache_load = True
    else:
        cache_load = False
    remove = REMOVE


def final_work():
    """做一些事后处理工作"""

    # ====================从内存中释放物种的STAR_Index====================
    if cache_load:
        command_list = [star, '--genomeLoad', 'Remove', '--genomeDir', organism_genome_index]
        command = ' '.join(command_list)
        os.system(command)
        text = 'release memory of %s Index' % organism
        print(add_color(text, 'yellow'))

    # ====================删除一些在当前目录下无缘无故出来的文件====================
    remove_file('./Aligned.out.sam')
    remove_file('./Log.out')
    remove_file('./Log.progress.out')
    os.system('rm -r _STARtmp')


def match_paired_fastq(fastq_list):
    """处理双端测序数据可能reads数目不匹配问题"""

    if len(fastq_list) == 2:
        file1 = fastq_list[0]
        file2 = fastq_list[1]
        file1_name = os.path.basename(file1)
        file2_name = os.path.basename(file2)
        file1_rows = int(os.popen(f'wc -l {file1}').read().split()[0])
        file2_rows = int(os.popen(f'wc -l {file2}').read().split()[0])
        if file1_rows != file2_rows:
            print(file1, file2, 'the rows of two paired-end files is not consistent, pre-process before handle')
            # 读取read首行，然后取交集
            file1_reads = os.popen(f'grep ^@SRR {file1}').read().strip().split('\n')
            file2_reads = os.popen(f'grep ^@SRR {file2}').read().strip().split('\n')
            shared_reads = set(file1_reads) & set(file2_reads)
            new_file1_name = 'new_' + file1_name
            new_file2_name = 'new_' + file2_name
            new_file1 = os.path.join(temp_path, new_file1_name)
            new_file2 = os.path.join(temp_path, new_file2_name)

            new_f = open(new_file1, 'w')
            data = []
            with open(file1) as f:
                for line, read in enumerate(f):
                    if read.strip() in shared_reads:
                        data.append(read.strip())
                        data.append(f.readline().strip())
                        data.append(f.readline().strip())
                        data.append(f.readline().strip())
                        if len(data) > 100000:
                            # 权衡内存和CPU，所以100000次进行一次IO
                            new_f.write('\n'.join(data)+'\n')
                            new_f.flush()
                            data = []
                    else:
                        f.readline()
                        f.readline()
                        f.readline()
            new_f.write('\n'.join(data)+'\n')
            new_f.close()

            new_f = open(new_file2, 'w')
            data = []
            with open(file2) as f:
                for line, read in enumerate(f):
                    if read.strip() in shared_reads:
                        data.append(read.strip())
                        data.append(f.readline().strip())
                        data.append(f.readline().strip())
                        data.append(f.readline().strip())
                        if len(data) > 100000:
                            new_f.write('\n'.join(data)+'\n')
                            new_f.flush()
                            data = []
                    else:
                        f.readline()
                        f.readline()
                        f.readline()
            new_f.write('\n'.join(data)+'\n')
            new_f.close()

            remove_file(file1)
            remove_file(file2)


def srr_handle(srr):
    """利用sratoolkit、STAR、featureCounts处理SRR文件"""

    # ====================查看文件是否未找到====================
    srr_file = os.path.join(data_path, srr)
    if not os.path.exists(srr_file):
        return 1, 'missed file%s' % srr_file

    # ====================运行sratoolkit====================
    # 运行sratoolkit，利用fasterq_dump将sra文件转化为fastq文件
    command_list = [fasterq_dump, '--split-3', '-e', str(cpu_use), srr_file, '-O', temp_path]
    command = ' '.join(command_list)
    foo = os.system(command)
    if foo:
        return foo, command
    # ====================运行STAR====================
    # 运行STAR，利用star将fastq文件转化为sam文件
    # 单端测序只有一个fastq，双端测序有两个fastq文件
    fastq_list = [os.path.join(temp_path, i) for i in os.listdir(temp_path) if i.endswith('fastq')]
    match_paired_fastq(fastq_list)
    fastq_list = [os.path.join(temp_path, i) for i in os.listdir(temp_path) if i.endswith('fastq')]
    command_list = [star, '--runThreadN', str(cpu_use), '--genomeDir',
                    organism_genome_index, '--outFileNamePrefix', temp_path+os.sep, '--readFilesIn']
    command_list.extend(fastq_list)
    if cache_load:
        command_list.extend(['--genomeLoad', 'LoadAndKeep'])
    command = ' '.join(command_list)
    foo = os.system(command)
    if foo:
        return foo, command
    # ====================运行featureCounts====================
    # 运行featureCounts，利用feature_counts将sam文件转化为最终的文件
    srr_name = srr.split('.')[0]
    command_list = [feature_counts, '-T', str(cpu_use), '-O', '-M', '-Q', '30', '-p', '-a', organism_genome_gtf,
                    os.path.join(temp_path, 'Aligned.out.sam'), '-o', os.path.join(result_path, srr_name)]
    command = ' '.join(command_list)
    foo = os.system(command)
    if foo:
        return foo, command
    # ====================删除summary文件====================
    remove_file(os.path.join(result_path, srr_name + '.summary'))
    # 删除原始sra文件，如果之前步骤没有处理成功，就不会删掉sra文件，方便调试
    if remove:
        remove_file(os.path.join(data_path, srr))

    # 返回0既是处理正确了
    return 0, ''


def main_work(task_file):
    """处理一个任务列表"""

    first_work(task_file)
    # ====================开始处理task_list里面的所有任务====================
    # 开始处理task_list里面的所有任务
    error = []
    for count, srr in enumerate(task_list, 1):
        srr_file = os.path.join(data_path, srr)
        file_size = round(os.path.getsize(srr_file) / 10 ** 9, 2)
        text = 'start to handle: %s  --  %s  size: %sG (%d/%d)' % (organism, srr, str(file_size), count, len(task_list))
        print(add_color(text, 'yellow'))

        foo, command = srr_handle(srr)
        if foo:
            text = 'error in %s --  %s\n command: %s' % (organism, srr, command)
            print(add_color(text, 'red'))
            error.append(srr)
        else:
            text = 'success in : %s --  %s' % (organism, srr)
            print(add_color(text, 'green'))

        # 把中间产生的临时文件及时删了，防止影响下步处理
        if len(os.listdir(temp_path)) > 0:
            os.system('rm -r ' + os.path.join(temp_path, '*'))
    # ====================输出这个进程处理的报告====================
    text = 'task: %s finished' % task_file
    print(add_color(text, 'green'))

    if error:
        text = '%s do not handle correctly' % str(error)
        print(add_color(text, 'red'))
    else:
        text = 'all srr have finished correctly'
        print(add_color(text, 'green'))
    # ========================================
    final_work()


if __name__ == '__main__':
    INPUT = sys.argv[1]  # INPUT：SRR_list
    main_work(INPUT)

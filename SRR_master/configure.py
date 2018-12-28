__doc__ = """SRR_master的配置文件"""

# ====================常用参数====================
SEP = '\t'  # 文件保存时的分隔符号，表达矩阵是用csv格式储存，不会因为因为这个参数而改变
PRINT_DETAIL = True  # 是否打印处理的详细信息
PRINT_WGET_DETAIL = True  # 是否打印wget的详细信息
CHECK = True  # 是否检查环境配置
REMOVE = True  # 是否删除处理过的原始SRR文件
THRESHOLD = 0.8  # 一个GSE中的SRR处理完成度大于这个阈值才会整合成表达矩阵
VALUE = 'TPM'  # 可以是RAW, RPKM, TPM

assert 0 < THRESHOLD <= 1
assert VALUE in ['RAW', 'RPKM', 'TPM']

# ====================软件和数据路径====================
PYTHON_PATH = 'python3'  # 执行Python的命令，如果不在环境变量中要给出绝对路径，python2不行
FASTERQ_DUMP = '/public/zhengll/xiongjh/software/sratoolkit-2.9.2/bin/fasterq-dump'  # fasterq-dump命令，必须给出绝对路径
STAR = '/public/zhengll/xiongjh/software/STAR-2.6.1a/bin/Linux_x86_64/STAR'  # STAR命令，必须给出绝对路径
FEATURE_COUNTS = '/public/zhengll/xiongjh/software/subread-1.6.2/bin/featureCounts'  # featureCounts命令，必须给出绝对路径
GENOME = '/public/zhengll/xiongjh/genome'  # 基因组的fasta和gtf文件的路径，必须给出绝对路径
GENOME_INDEX = '/public/zhengll/xiongjh/genomeIndex'  # 基因组的STAR_Index的路径，必须给出绝对路径

# # 以下是我放在天河二号上面的路径
# PYTHON_PATH = 'python3.6'  # 执行Python的命令，如果不在环境变量中要给出绝对路径
# FASTERQ_DUMP = '/NSFCGZ/nsfc2015_589/xiongjh/software/sratoolkit-2.9.2/bin/fasterq-dump'
# STAR = '/NSFCGZ/nsfc2015_589/xiongjh/software/STAR-2.6.1a/bin/Linux_x86_64/STAR'
# FEATURE_COUNTS = '/NSFCGZ/nsfc2015_589/xiongjh/software/subread-1.6.2/bin/featureCounts'
# GENOME = '/NSFCGZ/nsfc2015_589/xiongjh/genome'
# GENOME_INDEX = '/NSFCGZ/nsfc2015_589/xiongjh/genomeIndex'

# ====================进程配置====================
N_DOWNLOAD = 3  # 开启下载的进程数，越大越容易把带宽跑满，但是越消耗资源
# 将如下物种的STAR_Index加载入内存以加速处理，每个哺乳动物的Index大约30G
CACHE_LOAD = ['Homo_sapiens', 'Mus_musculus']
TASK_PER_WORKER = 100  # 每个进程处理的样本数
CPU_PER_WORKER = 10  # 每个进程用到的CPU
# 关于TASK_PER_WORKER、CPU_PER_WORKER设置的一些建议
# 处理SRR的进程数目会随TASK_PER_WORKER改变而改变，但可以用总共SRR个数除以TASK_PER_WORKER粗略估计
# 一般来说，考虑操作系统开关进程所需要的开销，一般一个含有24个CPU的进程串行处理100个文件，
# 速度要比6个进程，每个进程含有4个CPU，并行的处理100个相同文件要慢
N_INTEGRATION = 5  # 用于整合SRR成表达矩阵用的进程数目

# ====================天河二号配置====================
TIANHE = False  # 是否在天河二号上运行
TIANHE_NODE_NAME = 'nsfc3'  # 在天河二号上运行时用的节点名，如果小于60个节点，则用nsfc，如果大于60个节点，用BIOJOB
TIANHE_CPU_PER_WORKER = 24  # 在天河二号上设置了一个节点跑一个进程，一个节点刚好有24个核，设置为24就好了

# ====================以下参数一般保持默认就好了====================
DOWNLOAD_PREFIX = 'ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR'  # ftp下载SRR的前缀

# 白云山高，珠江水长
S = b'\x82' \
    b'?\x1c' \
    b'\x82Q\'' \
    b'\x81$L\x80' \
    b'V"\xf3\x05\n' \
    b'\x00\x0b\x13\x0f' \
    b'\x82/)\x80Y8\x83O' \
    b'W\x83W;\xbd\xa6\x84' \
    b'[5\x80TF\x83D\'\x81V+' \
    b'\x81-4\x81Z$\x82#.\x8b' \
    b'X(\x82N=\x828%\x81"5\xf3' \
    b'\x05\n\x00\x0b\x13\x0f\x80T' \
    b'\'\x84[,\x84=(\x836 \x822U\x82' \
    b'O1\xa6\x80Y"\x824K\x83D\'\x81V+' \
    b'\x81-4\x81L+\x81/A\x824K\x81Z$\x81' \
    b'2 \x84%K\x836 \x8bX(\x8315\x80V"\x80T' \
    b'\x1c\x80TF\x81L+\x82+,\x83@V\x8bX6\xa6' \
    b'\x83Y-\x80T&\x828%\xf3\x05\n\x00\x0b\x13' \
    b'\x0f\x80T\'\x84[,\x84=(\x836 \x13\x03\x01' \
    b'\x10\xca\x01\x14\x01\x81+K\x80WA\x80T\'\x84Y' \
    b'Y\xa6\x82$2\x84\x1c!\xf3\x05\n\x00\x0b\x13\x0f' \
    b'\x836 \xec\x0b\x13\x01\x0e\xef\x04\x01\x08\x08\x85' \
    b'#(\x859>\x828%\x84#F\x81TB\x828%\x13\x03\x01\x10\xa6' \
    b'\x81+K\x80WA\x81$E\x830D\x0f\x11\xfe\x0c\x0e\x0b\xff\x01' \
    b'\x0f\x0f\x82D=\x8193\x84L\x1f\x830D\xff\t\x00\x85#(\x859>' \
    b'\x836 \x81-Y\x80W@\xa6\x85\x1c6\x84[#\x820U\x81"5\x0f\x0e\x0e' \
    b'\xfb\x04\x01\x08\x0c\x01\x0e\xca\x0c\x15\x80T\'\x836 \x00\x0b\x13' \
    b'\n\x08\x0b\xfd\x00\xfb\x11\x0e\x08\x81#Y\x821L\xa6\x81LM\x81+K\x80WA\x81' \
    b'8D\xf3\x05\n\x00\x0b\x13\x0f\x80T\'\x84[,\x84=(\xbd\xa6\x83D\'\x81V+\x81-4\x81' \
    b'L+\x81/A\x84B\x1d\x80X-\x81\x1d#\x81*W\x80V"\x8bX(\x83A9\x80Y<\x82 %\x81[G\xbd\x8c;4&'

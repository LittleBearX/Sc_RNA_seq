__doc__ = """GSE_master的配置文件"""

SEP = '\t'  # 文件保存时的分隔符号，有些结果是用csv格式储存，不会因为因为这个参数而改变
N_WORKER = 10  # 处理GSE文件的进程数目
PRINT_DETAIL = True  # 是否打印处理的详细信息
# ModularityOptimizer程序的位置
MODULARITY_JAR = '/home/zhengll/xiongjh/workplace/Sc_RNA_seq/GSE_master/ModularityOptimizer.jar'

__doc__ = """此脚本为辅助脚本"""

import os
import time
import gzip


def log(func):
    """辅助打印分割线和记录时间"""

    line = '◆┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈◆ 华丽的分割线 ◆┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈◆'
    split_line = add_color(line, 'purple')

    def decorator(*args, **kwargs):
        t1 = time.time()
        print('\n%s' % split_line)
        result = func(*args, **kwargs)
        t2 = time.time()
        print('耗时：%.2fs' % (t2 - t1))
        return result

    return decorator


def make_dir(path):
    """创建文件夹"""

    try:
        os.makedirs(path)
    except FileExistsError:
        pass


def unpack_gz(file_name):
    """解压gz压缩，直接返回数据，不保存到文件"""

    with gzip.GzipFile(file_name) as gz_file:
        data = gz_file.read()
    if isinstance(data, bytes):
        data = data.decode('utf8')
    return data


def add_color(text, color='r'):
    """为字符串带上颜色（Linux）"""

    color_palette = {'black': 30, 'red': 31, 'green': 32, 'yellow': 33, 'blue': 34, 'purple': 35, 'cyan': 36,
                     'white': 37, 'r': 31, 'g': 32, 'y': 33, 'b': 34, 'p': 35, 'c': 36, 'w': 37}

    if color in color_palette:
        text = str(text)
        color_mode = '1;%dm' % color_palette[color]
        color_text = '\033[%s%s\033[0m' % (color_mode, text)
        return color_text
    else:
        print('color必须是%s' % set(color_palette.keys()))
        return text

__doc__ = """分析表达矩阵的主要算法"""

import os
import sys
import time
import warnings
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn import manifold
import matplotlib.pyplot as plt
from scipy.stats import ranksums
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors

# 读取配置文件
from Sc_RNA_seq.GSE_master.configure import PRINT_DETAIL
from Sc_RNA_seq.GSE_master.configure import MODULARITY_JAR
plt.style.use('seaborn')
warnings.filterwarnings("ignore")


def time_log(func, print_detail=PRINT_DETAIL):
    """记录每一步耗时"""

    def wrapper(self, *args, **kw):
        if print_detail:
            time1 = time.time()
            func(self, *args, **kw)
            time2 = time.time()
            print(f'{self.name}: {func.__name__}', end='')
            print('\tTime: %9.3fs' % (time2 - time1))
        else:
            func(self, *args, **kw)

    return wrapper


class Project:
    """分析表达矩阵"""

    __slots__ = ['path', 'name', 'raw_data', 'data', 'cells_numi', 'genes_use', 'cells_use', 'genes_num', 'cells_num',
                 'data_norm', 'data_log1p', 'data_scale', 'genes_high_var', 'cells_pca', 'corr_matrix', 'snn_matrix',
                 'cells_ident', 'cells_use_sort_by_ident', 'cells_tsne_2D', 'cells_tsne_3D', 'all_marker_genes',
                 'heat_map']

    def __init__(self, file_path):
        """初始化并读入原始数据"""

        self.path = os.path.dirname(file_path)
        self.name = os.path.split(self.path)[1]
        with open(file_path) as f:
            self.raw_data = pd.read_csv(f, index_col=0)
        self.data = None
        self.cells_numi = None
        self.genes_use, self.cells_use = None, None
        self.genes_num, self.cells_num = None, None
        self.data_norm = None
        self.data_log1p = None
        self.data_scale = None
        self.genes_high_var = None
        self.cells_pca = None
        self.corr_matrix = None
        self.snn_matrix = None
        self.cells_ident = None
        self.cells_use_sort_by_ident = None
        self.cells_tsne_2D = None
        self.cells_tsne_3D = None
        self.all_marker_genes = None
        self.heat_map = None

    @time_log
    def check(self):
        """检查数据和参数是否出错"""

        genes_num, cells_num = self.raw_data.shape
        if genes_num < 30:
            raise TypeError(f'{self.name}: 基因太少了，只有{genes_num}个!')
        if cells_num < 30:
            raise TypeError(f'{self.name}: 细胞太少了，只有{cells_num}个!')
        if len(set(self.raw_data.index)) < genes_num:
            raise TypeError(f'{self.name}: 基因名有重复!')
        if self.raw_data.min().min() < 0:
            raise TypeError(f'{self.name}: 矩阵中有负数!')
        if not os.path.exists(MODULARITY_JAR):
            raise ModuleNotFoundError(f'{os.path.abspath(MODULARITY_JAR)}: 没有找到ModularityOptimizer.jar！')

    @time_log
    def filter_data(self, min_genes_percent=0.001, min_cells_percent=0.005, max_cells_percent=0.995):
        """
        对低表达细胞和在大多细胞没有表达的基因进行过滤
        :param min_genes_percent: 在大多细胞中没有表达的基因将被过滤
        :param min_cells_percent: 表达基因数较少的细胞将被过滤
        :param max_cells_percent: 表达基因数较多的细胞将被过滤
        """

        genes_num, cells_num = self.raw_data.shape
        data_non_zero = self.raw_data > 0.00001
        cells_count = data_non_zero.sum(axis=0)  # 细胞表达多少个基因
        genes_count = data_non_zero.sum(axis=1)  # 基因在多少个细胞中表达
        cells_numi = self.raw_data.sum(axis=0)  # 细胞表达多少条RNA
        min_cells_rank, max_cells_rank = int(cells_num * min_cells_percent), int(cells_num * max_cells_percent)
        min_cells_count, max_cells_count = cells_count.sort_values()[[min_cells_rank, max_cells_rank]]
        min_genes_count = max(int(cells_num * min_genes_percent), 1)
        # 根据参数过滤
        self.cells_use = self.raw_data.columns[(cells_count > min_cells_count) & (cells_count < max_cells_count)]
        self.genes_use = self.raw_data.index[genes_count > min_genes_count]
        # 过滤后更新数据
        self.data = self.raw_data.loc[self.genes_use, self.cells_use]
        self.cells_numi = cells_numi[self.cells_use]
        self.genes_num, self.cells_num = self.data.shape

    @time_log
    def normalize_data(self, scale_factor=10000):
        """
        对细胞表达数据进行初步标准化
        :param scale_factor: 缩放比例
        """

        self.data_norm = (self.data * scale_factor / self.cells_numi)
        self.data_log1p = np.log1p(self.data_norm)

    @time_log
    def find_variable_genes(self, bin_nums=20, x_low_cutoff=0.0125, x_high_cutoff=3, y_low_cutoff=0.5):
        """
        找出高变异基因
        :param bin_nums: 根据基因表达水平对基因进行分bins
        :param x_low_cutoff: 不处理表达量低于此值的基因
        :param x_high_cutoff: 不处理表达量高于此值的基因
        :param y_low_cutoff: genes_bin_var高于此值的基因认为是高变异基因
        """

        foo = self.data_norm.mean(axis=1)
        genes_mean = np.log1p(foo)
        genes_dispersion = np.log(self.data_norm.var(axis=1, ddof=1) / foo).fillna(0)  # 基因表达方差/基因表达均值
        genes_bin = pd.cut(genes_mean, bin_nums)  # 按基因在不同细胞中平均表达水平分bins
        bin_mean = genes_dispersion.groupby(genes_bin).mean().fillna(0)  # 每个bin中基因的dispersion的平均
        bin_std = genes_dispersion.groupby(genes_bin).std(ddof=1).fillna(0)  # 每个bin中基因的dispersion的方差
        genes_bin_var = ((genes_dispersion - bin_mean[genes_bin].get_values()) / bin_std[genes_bin].get_values())
        genes_bin_var.fillna(0)
        foo = (x_low_cutoff < genes_mean) & (genes_mean < x_high_cutoff) & (genes_bin_var > y_low_cutoff)
        self.genes_high_var = self.genes_use[foo]

    @time_log
    def scale_data(self, scale_max=10):
        """
        用线性回归消除细胞整体表达高低对表达量的影响，并进一步标准化数据，方差为1，均值为0
        :param scale_max: 标准化后数据不能超过这个值
        """

        foo_list = []
        x = np.mat(np.row_stack((self.cells_numi, np.ones(self.cells_num)))).T
        x_inverse = x.I
        for i in range(self.genes_num):
            y = np.mat(np.array(self.data_log1p.iloc[i, :])).T
            coef = x_inverse * y  # 最小二乘求系数
            residuals = y - x * coef
            foo = (residuals - np.mean(residuals)) / np.std(residuals, ddof=1)  # 标准化，方差为1，均值为0
            foo[foo > scale_max] = scale_max
            foo_list.append(foo)
        self.data_scale = pd.DataFrame(np.hstack(foo_list).T, index=self.genes_use, columns=self.cells_use)
        self.data_scale.fillna(0, inplace=True)

    @time_log
    def run_pca(self, only_hvg=True, n="auto", corr=True):
        """
        PCA降维
        :param only_hvg: True的话，只以高变异基因为特征进行降维
        :param n: 主成分个数，'auto'的话，根据数据规模自动确定值
        :param corr: 是否计算相关系数矩阵
        """

        if only_hvg:
            cells_high_var = self.data_scale.loc[self.genes_high_var, :]
            x = cells_high_var.T
        else:
            x = self.data_scale.T
        if n == 'auto':
            n = int(np.log2(self.cells_num / 10000 + 1) * 14 + 8)
        pca = PCA(n_components=n)
        pca.fit(x)
        pca_embeddings = pca.transform(x)[:, :n]
        self.cells_pca = {'variance_ratio': pca.explained_variance_ratio_, 'embeddings': pca_embeddings}
        if corr:
            self.corr_matrix = np.corrcoef(self.cells_pca['embeddings'])

    @time_log
    def find_clusters(self, k=30, do_fast=0, snn=True, resolution='auto', modularity=1, algorithm=1, n_start=100,
                      n_iter=10, random_seed=0):
        """
        利用SNN距离，进行社团聚类算法
        :param k: 算法寻找k个近邻
        :param do_fast: 是否加速计算
        :param snn: 是否计算snn矩阵
        :param resolution: 粒度，越大分的越细，'auto'的话，根据数据规模自动确定值
        :param modularity: 最大模块度
        :param algorithm: 算法类型
        :param n_start: 进行多少次运算
        :param n_iter: 每次运算迭代多少次
        :param random_seed: 随机种子
        """

        # =============计算SNN相似度=============
        # 用kd树计算k最近邻
        k = min(k, self.cells_num)
        knn = NearestNeighbors(algorithm="kd_tree").fit(self.cells_pca['embeddings'])
        knn_matrix = knn.kneighbors(self.cells_pca['embeddings'], n_neighbors=k, return_distance=False)
        knn_set = [set(i) for i in knn_matrix.tolist()]
        edge_data = []

        # 完整算法，复杂度是O(细胞数*细胞数)
        if do_fast == 0:
            for i in range(self.cells_num):
                i_set = knn_set[i]
                for j in range(i + 1, self.cells_num):
                    count = len(i_set & knn_set[j])
                    if count > 3:
                        edge_value = (count / (k * 2 - count))
                        edge_data.append((i, j, edge_value))
        # 近似算法，保留85%精确度，复杂度是O(k*k*细胞数)
        elif do_fast == 1:
            knn_check_list = []
            for i in range(self.cells_num):
                neighbors = knn_matrix[i]
                neighbors_neighbors = set()
                for neighbor in neighbors:
                    for neighbor_neighbor in knn_matrix[neighbor]:
                        if neighbor_neighbor > i:
                            neighbors_neighbors.add(neighbor_neighbor)
                knn_check_list.append(neighbors_neighbors)
            for i in range(self.cells_num):
                i_set = knn_set[i]
                for j in knn_check_list[i]:
                    count = len(i_set & knn_set[j])
                    if count > 3:
                        edge_value = (count / (k * 2 - count))
                        edge_data.append((i, j, edge_value))
        # 更一步近似，保留40%精度，复杂度是O(k*细胞数)
        elif do_fast == 2:
            for i in range(self.cells_num):
                i_set = knn_set[i]
                check_set = {_i for _i in i_set if _i > i}
                for j in check_set:
                    count = len(i_set & knn_set[j])
                    if count > 3:
                        edge_value = (count / (k * 2 - count))
                        edge_data.append((i, j, edge_value))
        else:
            raise ValueError('Parameter: do_fast must be in {0,1,2}')

        if snn:
            snn_matrix = np.diag([1.] * self.cells_num)
            for i, j, edge_value in edge_data:
                snn_matrix[i][j] = snn_matrix[j][i] = edge_value
            self.snn_matrix = snn_matrix

        # =============社团发现算法=============
        edge_file = os.path.join(self.path, "edge_data.txt")
        cluster_file = os.path.join(self.path, "cluster.txt")
        with open(edge_file, "w", encoding="utf8") as file:
            edge_data = [f'{i}\t{j}\t{edge_value}' for i, j, edge_value in edge_data]
            file.write('\n'.join(edge_data))
        if resolution == 'auto':
            resolution = np.log10(self.cells_num / 5000 + 1) * 1.6 + 0.3
        command_list = ["java -jar", MODULARITY_JAR, edge_file, cluster_file, str(modularity), str(resolution),
                        str(algorithm), str(n_start), str(n_iter), str(random_seed), '1']
        command = ' '.join(command_list)
        os.popen(command)
        while True:
            try:
                with open(cluster_file, "r", encoding="utf8") as file:
                    cells_ident = [int(i) for i in file]
                    break
            except FileNotFoundError:
                time.sleep(0.5)

        os.remove(edge_file)
        os.remove(cluster_file)

        # 将单独成一类的细胞根据近邻细胞类别进行重新划分
        single_idents = [i for i in set(cells_ident) if cells_ident.count(i) == 1]
        for i in single_idents:
            ident_index = cells_ident.index(i)
            ident_knn = knn_matrix[ident_index][1:]
            new_ident = np.bincount(np.array(cells_ident)[ident_knn]).argmax()
            cells_ident[ident_index] = new_ident

        self.cells_ident = np.array(cells_ident)
        ident_dict = dict(zip(self.cells_use, self.cells_ident))
        self.cells_use_sort_by_ident = sorted(self.cells_use, key=lambda x: ident_dict[x])

    @time_log
    def run_tsne(self, n_iter=1000, three=True):
        """
        tSNE可视化
        :param n_iter: 迭代次数
        :param three: 是否三维可视化
        """

        tsne = manifold.TSNE(n_components=2, random_state=0, n_iter=n_iter, angle=0.5)
        self.cells_tsne_2D = tsne.fit_transform(self.cells_pca['embeddings'])
        if three:
            tsne = manifold.TSNE(n_components=3, random_state=0, n_iter=n_iter, angle=0.5)
            self.cells_tsne_3D = tsne.fit_transform(self.cells_pca['embeddings'])

    @time_log
    def find_all_marker_genes(self, min_pct=0.25, min_diff_pct=0.0, log_threshold=0.25, only_pos=True, min_pvalue=0.01):
        """
        one vs rest，秩和检验找不同类的差异表达基因
        :param min_pct: 低于这个比例的细胞数表达的基因被排除
        :param min_diff_pct: 细胞数表达比例数之差低于这个的基因被排除
        :param log_threshold: 细胞数表达量之差低于这个的基因被排除
        :param only_pos: 是否只找每个类中上调表达的基因
        :param min_pvalue: 秩和检验p-value经过Bonferroni修正后小于此的基因为marker基因
        """

        all_idents = set(self.cells_ident)
        if len(all_idents) > 1:
            marker_genes_list = []
            for i in all_idents:
                foo = (self.cells_ident == i)
                ident1, ident2 = (self.cells_use[foo], self.cells_use[~foo])
                ident1_num, ident2_num = (len(ident1), len(ident2))
                ident1_count = (self.data_norm.loc[:, ident1] > 0).sum(axis=1) / ident1_num
                ident2_count = (self.data_norm.loc[:, ident2] > 0).sum(axis=1) / ident2_num
                ident_count_diff = ident1_count - ident2_count
                new_genes_use = self.genes_use[
                    ((ident1_count > min_pct) | (ident2_count > min_pct)) & (ident_count_diff.abs() > min_diff_pct)]
                ident1_data = self.data_norm.loc[new_genes_use, ident1]
                ident2_data = self.data_norm.loc[new_genes_use, ident2]
                ident1_mean = np.log1p(np.mean(ident1_data, axis=1))
                ident2_mean = np.log1p(np.mean(ident2_data, axis=1))
                ident_mean_diff = ident1_mean - ident2_mean
                if only_pos:
                    new_genes_use = np.array(new_genes_use[ident_mean_diff > log_threshold])
                else:
                    new_genes_use = np.array(new_genes_use[ident_mean_diff.abs() > log_threshold])
                ident1_data = self.data_log1p.loc[new_genes_use, ident1]
                ident2_data = self.data_log1p.loc[new_genes_use, ident2]
                pvalue = np.array([ranksums(ident1_data.iloc[i, :], ident2_data.iloc[i, :]).pvalue
                                   for i in range(len(new_genes_use))])
                adj_pvalue = pvalue * self.genes_num
                filter_pvalue = adj_pvalue < min_pvalue
                pvalue = pvalue[filter_pvalue]
                adj_pvalue = adj_pvalue[filter_pvalue]
                new_genes_use = new_genes_use[filter_pvalue]
                ident_data = [[i] * len(new_genes_use), new_genes_use, pvalue, adj_pvalue,
                              ident1_count[new_genes_use], ident2_count[new_genes_use], ident_count_diff[new_genes_use],
                              ident1_mean[new_genes_use], ident2_mean[new_genes_use], ident_mean_diff[new_genes_use]]
                ident_data = np.column_stack(ident_data)
                sorter = np.lexsort([-ident_mean_diff[new_genes_use], adj_pvalue])
                marker_genes_list.append(ident_data[sorter])
            columns = ['cluster', 'gene', 'p_value', 'adjusted_p_value', 'count.1', 'count.2', 'count.diff', 'mean.1',
                       'mean.2', 'mean.diff']
            all_marker_genes = pd.DataFrame(np.row_stack(marker_genes_list), columns=columns)
            if all_marker_genes.shape[0] > 0:
                self.all_marker_genes = all_marker_genes

    @time_log
    def run_heat_map(self, top=10, sort_by='mean.diff'):
        """
        计算heap_map
        :param top: 每个类中与其他类中差异最大几个基因
        :param sort_by: 根据这个来判断差异
        """
        if (self.all_marker_genes is not None) and (self.all_marker_genes.shape[0] > 0):
            # 根据sort_by选择每个类别细胞中和其他类别细胞表达差别最大的基因
            heatmap_genes = pd.concat([i[1].sort_values(sort_by, ascending=False)['gene'][:top]
                                       for i in self.all_marker_genes.groupby('cluster')])
            heat_map = self.data_scale.loc[heatmap_genes, self.cells_use_sort_by_ident]
            self.heat_map = heat_map

    @time_log
    def save_picture(self, pca=True, tsne=True, snn_matrix=True, corr_matrix=True, heat_map=True, dpi=300):
        """保存图片"""

        idents_count = len(set(self.cells_ident))
        color = sns.husl_palette(idents_count, h=0, s=0.75)

        if pca and self.cells_pca is not None:
            # 画前两个主成分
            x = self.cells_pca['embeddings'][:, 0]
            y = self.cells_pca['embeddings'][:, 1]
            plt.figure(figsize=(8, 8))
            plt.title('PCA', size=14)
            for i in range(idents_count):
                mask = self.cells_ident == i
                ident_x = x[mask].tolist()
                ident_y = y[mask].tolist()
                plt.scatter(ident_x, ident_y, s=5, color=color[i], alpha=0.8)
            plt.xlabel('PC-1')
            plt.ylabel('PC-2')
            bar = [plt.plot([], 'o', markersize=8, color=color[i], alpha=0.8)[0] for i in range(idents_count)]
            ax = plt.legend(bar, range(idents_count), loc='upper right')
            ax.set_title(r'# Cells idents')
            plt.savefig(os.path.join(self.path, 'PCA.png'), dpi=dpi)

        if tsne and self.cells_tsne_2D is not None:
            # 画二维tsne
            x = self.cells_tsne_2D[:, 0]
            y = self.cells_tsne_2D[:, 1]
            plt.figure(figsize=(8, 8))
            plt.title('tSNE', size=14)
            for i in range(idents_count):
                mask = self.cells_ident == i
                ident_x = x[mask].tolist()
                ident_y = y[mask].tolist()
                plt.scatter(ident_x, ident_y, s=5, color=color[i], alpha=0.8)
            plt.xlabel('tSNE-1')
            plt.ylabel('tSNE-2')
            bar = [plt.plot([], 'o', markersize=8, color=color[i], alpha=0.8)[0] for i in range(idents_count)]
            ax = plt.legend(bar, range(idents_count), loc='upper right')
            ax.set_title(r'# Cells idents')
            plt.savefig(os.path.join(self.path, 'tSNE.png'), dpi=dpi)

        if corr_matrix and self.corr_matrix is not None:
            plt.figure(figsize=(8, 8))
            plt.title("Correlation coefficient")
            foo = pd.DataFrame(self.corr_matrix, index=self.cells_use, columns=self.cells_use)
            sns.heatmap(foo.loc[self.cells_use_sort_by_ident, self.cells_use_sort_by_ident], cmap='afmhot',
                        yticklabels=False, xticklabels=False)
            plt.savefig(os.path.join(self.path, 'corr.png'), dpi=dpi)

        if snn_matrix and self.snn_matrix is not None:
            plt.figure(figsize=(8, 8))
            plt.title("SNN")
            foo = pd.DataFrame(self.snn_matrix, index=self.cells_use, columns=self.cells_use)
            sns.heatmap(foo.loc[self.cells_use_sort_by_ident, self.cells_use_sort_by_ident], cmap='afmhot',
                        yticklabels=False, xticklabels=False)
            plt.savefig(os.path.join(self.path, 'snn.png'), dpi=dpi)

        if heat_map and self.heat_map is not None:
            plt.figure(figsize=(12, 8))
            plt.title('Heat map')
            sns.heatmap(self.heat_map, cmap='afmhot', yticklabels=False, xticklabels=False)
            plt.savefig(os.path.join(self.path, 'heat_map.png'), dpi=dpi)

    @time_log
    def save_data(self, hvg_genes=True, scale_data=True, corr_matrix=True, snn_matrix=True, pca=True, tsne=True,
                  marker_genes=True, heat_map=True, precision=3):
        """保存数据"""

        if hvg_genes and self.genes_high_var is not None:
            with open(os.path.join(self.path, 'hvg_genes.csv'), 'w', encoding='utf8') as f:
                f.write('\n'.join(self.genes_high_var))

        if scale_data and self.data_scale is not None:
            with open(os.path.join(self.path, 'scale_data.csv'), 'w', encoding='utf8') as f:
                self.data_scale.round(precision).to_csv(f)

        if corr_matrix and self.corr_matrix is not None:
            with open(os.path.join(self.path, 'corr.csv'), 'w', encoding='utf8') as f:
                foo = pd.DataFrame(self.corr_matrix, index=self.cells_use, columns=self.cells_use)
                foo.loc[self.cells_use_sort_by_ident, self.cells_use_sort_by_ident].round(precision).to_csv(f)

        if snn_matrix and self.snn_matrix is not None:
            with open(os.path.join(self.path, 'snn.csv'), 'w', encoding='utf8') as f:
                foo = pd.DataFrame(self.snn_matrix, index=self.cells_use, columns=self.cells_use)
                foo.loc[self.cells_use_sort_by_ident, self.cells_use_sort_by_ident].round(precision).to_csv(f)

        if pca and self.cells_pca is not None:
            with open(os.path.join(self.path, 'pca.csv'), 'w', encoding='utf8') as f:
                foo = pd.DataFrame(np.row_stack([self.cells_pca['variance_ratio'], self.cells_pca['embeddings']]),
                                   index=['variance_ratio'] + self.cells_use.tolist(),
                                   columns=['PC' + str(i + 1) for i in range(len(self.cells_pca['variance_ratio']))])
                foo['cluster'] = [sum(self.cells_pca['variance_ratio'])] + self.cells_ident.tolist()
                foo.to_csv(f)

        if tsne:
            if self.cells_tsne_2D is not None:
                with open(os.path.join(self.path, 'tsne2D.csv'), 'w', encoding='utf8') as f:
                    pd.DataFrame(np.hstack([self.cells_ident.reshape([-1, 1]), self.cells_tsne_2D]),
                                 index=self.cells_use, columns=['cluster', 'x', 'y']).to_csv(f)
            if self.cells_tsne_3D is not None:
                with open(os.path.join(self.path, 'tsne3D.csv'), 'w', encoding='utf8') as f:
                    pd.DataFrame(np.hstack([self.cells_ident.reshape([-1, 1]), self.cells_tsne_3D]),
                                 index=self.cells_use, columns=['cluster', 'x', 'y', 'z']).to_csv(f)

        if marker_genes and self.all_marker_genes is not None:
            with open(os.path.join(self.path, 'marker_genes.csv'), 'w', encoding='utf8') as f:
                self.all_marker_genes.to_csv(f, index=False)

        if heat_map and self.heat_map is not None:
            with open(os.path.join(self.path, 'heat_map.csv'), 'w', encoding='utf8') as f:
                self.heat_map.to_csv(f)

    def compute_all_steps(self):
        """计算所有步骤"""

        self.check()
        self.filter_data()
        self.normalize_data()
        self.find_variable_genes()
        self.scale_data()
        self.run_pca()
        self.find_clusters()
        self.run_tsne()
        self.find_all_marker_genes()
        self.run_heat_map()
        self.save_picture()
        self.save_data()


if __name__ == '__main__':
    INPUT = sys.argv[1]  # INPUT：matrix.csv
    project = Project(INPUT)
    project.compute_all_steps()

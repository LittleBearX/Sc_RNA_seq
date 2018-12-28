import numpy as np
import pandas as pd
from sklearn.metrics import mean_absolute_error


def test_new_feature(train_x, train_y, test_x, test_y, model, fun, features_list,
                     corr_threshold=0, promotion_threshold=0.00001):
    """
    用fun组合features里的特征得到新特征，然后加入新特征测试预测效果是否变好
    :param train_x: 训练集特征dataframe
    :param train_y: 训练集目标向量
    :param test_x: 测试集特征dataframe
    :param test_y: 测试集特征dataframe
    :param model: 用于测试新特征的机器学习模型
    :param fun: 用于组合特征的函数
    :param features_list: 用于组合的特征列表
    :param corr_threshold: 新特征和目标的值的相关系数阈值
    :param promotion_threshold: 效果提升阈值
    """

    print('函数：', fun.__name__, '\t特征：', features_list)

    new_train_feature = fun(train_x, features_list)  # 组合得到新特征
    corr = np.corrcoef(new_train_feature, train_y)[1, 0]  # 新特征和目标的相关系数
    print('新特征和目标的相关系数：%f' % corr)

    # 相关系数的绝对值如果大于阈值，运行机器学习模型进行预测
    if abs(corr) > corr_threshold:
        # 用旧特征预测
        model.fit(train_x, train_y)
        old_predict_y = model.predict(test_x)
        old_mae = mean_absolute_error(test_y, old_predict_y)

        # 为训练集和测试集加入新特征，然后预测
        new_train_x = train_x.copy()
        new_train_x['new_feature'] = new_train_feature
        new_test_feature = fun(test_x, features_list)
        new_test_x = test_x.copy()
        new_test_x['new_feature'] = new_test_feature

        # 用新特征预测
        model.fit(new_train_x, train_y)
        new_predict_y = model.predict(new_test_x)
        new_mae = mean_absolute_error(test_y, new_predict_y)

        print('用旧特征预测mae：%f\n新特征预测mae：%f' % (old_mae, new_mae))
        promotion = old_mae - new_mae

        if promotion > promotion_threshold:
            print('\U0001F60A 效果得到了提升！')
            return [fun.__name__, features_list, promotion]
        else:
            return None
    else:
        print('新特征和目标的相关系数低于阈值：%f' % corr_threshold)
        return None


def multi_process(task_list, train_x, train_y, test_x, test_y, model, fun,
                  corr_threshold=0, promotion_threshold=0.00001):
    foo_list = []
    for features_list in task_list:
        result = test_new_feature(train_x, train_y, test_x, test_y, model, fun, features_list,
                                  corr_threshold=corr_threshold, promotion_threshold=promotion_threshold)
        if result:
            foo_list.append(result)
    return foo_list

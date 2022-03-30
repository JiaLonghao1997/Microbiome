"""
reproduce of the result of classification
study to study  transfer
leave one study out
using GridSearchCV method to search for the best C for LR
"""

import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GridSearchCV
import warnings
from sklearn.svm import SVC

tag = 'species'
Study = ['FR-CRC', 'AT-CRC', 'CN-CRC', 'US-CRC', 'DE-CRC']
Study1 = ['CN-CRC']
log_n0 = 1e-5
sd_min_q = 10
kernel = 'linear'

parameters_stst = {'FR-CRC': 1, 'AT-CRC': 1, 'CN-CRC': 1, 'US-CRC': 1, 'DE-CRC': 1}

warnings.filterwarnings("ignore")
species_url = 'data//' + tag + '//feat_rel_crc.tsv'
meta_url = 'data//meta//meta_crc.tsv'
species = pd.read_csv(species_url, sep='\t', header=0)
meta = pd.read_csv(meta_url, sep='\t', header=0)

# row name
sample_name = meta['Sample_ID'].tolist()
# predict_matrix store the prediction result
predict_matrix = pd.DataFrame(np.zeros(shape=(species.shape[1], 6)), index=sample_name,
                              columns=['FR-CRC', 'AT-CRC', 'CN-CRC', 'US-CRC', 'DE-CRC', 'LOSO'])
accuracy_matrix = pd.DataFrame(np.zeros(shape=(7, 5)),
                               index=['FR-CRC', 'AT-CRC', 'CN-CRC', 'US-CRC', 'DE-CRC', 'LOSO', 'LOSO_CV']
                               , columns=['FR-CRC', 'AT-CRC', 'CN-CRC', 'US-CRC', 'DE-CRC'])
recall_matrix = pd.DataFrame(np.zeros(shape=(7, 5)),
                             index=['FR-CRC', 'AT-CRC', 'CN-CRC', 'US-CRC', 'DE-CRC', 'LOSO', 'LOSO_CV']
                             , columns=['FR-CRC', 'AT-CRC', 'CN-CRC', 'US-CRC', 'DE-CRC'])

"""
train model for study to study transfer
"""
for study in Study:
    print('training on {} with Study to Study Transfer'.format(study))

    train_data = pd.read_csv('pre_data//species_stst//{}//{}_data.csv'.format(study, study), header=0, index_col=0)
    train_sample_id = train_data._stat_axis.values.tolist()  ###行名是sample_id

    y = train_data['Label']
    x = train_data.drop(columns='Label')
    x = x.values
    y = y.values
    train_pre_matrix = pd.DataFrame(np.zeros(shape=(train_data.shape[0], 10)),
                                    index=train_sample_id, columns=['1', '2', '3', '4', '5', '6', '7', '8', '9', '10'])

    right_num_study = 0
    right_num_positive = 0
    # parameters = {'C':[0.1,1,5,10,100,1000]}
    for i in range(10):
        skf = StratifiedKFold(n_splits=10)
        for train_index, test_index in skf.split(x, y):
            train_x = x[train_index]
            train_y = y[train_index]
            test_x = x[test_index]
            test_y = y[test_index]

            """
                  remove features with std=0
                  reference: https://www.cnblogs.com/pinard/p/9032759.html
            """
            std_preprocess = np.std(train_x, axis=0)
            train_x = train_x[:, std_preprocess != 0]
            train_x = np.log10(train_x + log_n0)
            mean = np.mean(train_x, axis=0)
            std = np.std(train_x, axis=0)
            q = np.percentile(std, 10)
            train_x = (train_x - mean) / (std + q)

            test_x = test_x[:, std_preprocess != 0]
            test_x = np.log10(test_x + log_n0)
            test_x = (test_x - mean) / (std + q)

            # lr = LogisticRegression(penalty='l1',solver='liblinear')
            # clf = GridSearchCV(lr,parameters,cv=5,n_jobs=-1,pre_dispatch=5)
            # clf.fit(train_x,train_y)
            # lr = LogisticRegression(penalty='l1',solver='liblinear',C=parameters_stst[study])
            # lr.fit(train_x,train_y)

            svc = SVC(kernel=kernel, C=parameters_stst[study], probability=True)
            svc.fit(train_x, train_y)

            proba = svc.predict_proba(test_x)[:, 1]
            train_sample_id = np.array(train_sample_id)
            #print(svc.predict_proba(test_x))
            ##predict_proba返回的是一个 n 行 k 列的数组， 第 i 行 第 j 列上的数值是模型预测 第 i 个预测样本为某个标签的概率，并且每一行的概率和为1。
            test_sample = train_sample_id[test_index]
            train_pre_matrix.loc[test_sample, '%d' % (i + 1)] = proba
            right_num_study += np.sum(test_y == svc.predict(test_x))
            pre_y = svc.predict(test_x)

            for m in range(len(test_x)):
                if test_y[m] == 1 and test_y[m] == pre_y[m]:
                    right_num_positive += 1

            for study_test in Study:
                if study_test == study:
                    continue

                else:
                    test_stst = pd.read_csv('pre_data//species_stst//{}//{}_data.csv'.format(study_test, study_test),
                                            header=0, index_col=0)
                    test_stst_index = test_stst._stat_axis.values.tolist()
                    test_stst_y = test_stst['Label']
                    test_stst_x = test_stst.drop(columns='Label')
                    test_stst_x = test_stst_x.values
                    test_stst_x = test_stst_x[:, std_preprocess != 0]
                    test_stst_x = np.log10(test_stst_x + log_n0)
                    test_stst_x = (test_stst_x - mean) / (std + q)
                    proba_stst = svc.predict_proba(test_stst_x)[:, 1] / 100
                    predict_matrix.loc[test_stst_index, study] += proba_stst
                    score = svc.score(test_stst_x, test_stst_y)
                    predict_y = svc.predict(test_stst_x)
                    num = 0
                    print(len(test_stst_x))
                    for k in range(len(test_stst_x)):
                        if test_stst_y[k] == 1 and predict_y[k] == 1:
                            num += 1
                    recall = num / np.sum(test_stst_y == 1)
                    # recall = np.sum(predict_y == 1) / np.sum(test_stst_y == 1)
                    accuracy_matrix.loc[study, study_test] += score / 100
                    recall_matrix.loc[study, study_test] += recall / 100

        train_pre_mean = train_pre_matrix.mean(1)
        predict_matrix.loc[train_pre_matrix._stat_axis.values.tolist(), study] = train_pre_mean
        accuracy_matrix.loc[study, study] += right_num_study / (len(x) * 10)
        recall_matrix.loc[study, study] += right_num_positive / (np.sum(y == 1) * 10)
        right_num_study = 0
        right_num_positive = 0

"""
LOSO
"""
for study in Study:
    print('testing on study {} with LOSO'.format(study))
    parameters_loso = {'FR-CRC': 10, 'AT-CRC': 10, 'CN-CRC': 0.1, 'US-CRC': 1, 'DE-CRC': 10}
    train_data = pd.read_csv('pre_data//species_loso//{}//train_data.csv'.format(study), header=0, index_col=0)
    train_sample_id = train_data._stat_axis.values.tolist()
    y = train_data['Label']
    x = train_data.drop(columns='Label')
    x = x.values
    y = y.values
    #parameters = {'C': [0.1,  1, 5, 10, 100,1000]}

    for i in range(10):
        right_num_study = 0
        right_num_positive = 0
        skf = StratifiedKFold(n_splits=10)
        for train_index,test_index in skf.split(x,y):
            train_x = x[train_index]
            train_y = y[train_index]
            std_preprocess = np.std(train_x, axis=0)
            train_x = train_x[:, std_preprocess != 0]
            train_x = np.log10(train_x + log_n0)
            mean = np.mean(train_x, axis=0)
            std = np.std(train_x, axis=0)
            q = np.percentile(std, 10)
            train_x = (train_x - mean) / (std + q)

            # lr = LogisticRegression(penalty='l1',solver='liblinear')
            # clf = GridSearchCV(lr,parameters,cv=5,n_jobs=-1,pre_dispatch=5)
            # clf.fit(train_x,train_y)
            #lr = LogisticRegression(penalty='l1',solver='liblinear',C=0.1)
            #lr.fit(train_x,train_y)
            #svc = SVC(kernel=kernel,C=1,probability=True)
            #svc.fit(train_x,train_y)
            lr = LogisticRegression(penalty='l1', solver='liblinear', n_jobs=-1, C=parameters_stst[study])
            lr.fit(train_x, train_y)

            test_x = x[test_index]
            test_y = y[test_index]
            test_x = test_x[:, std_preprocess != 0]
            test_x = np.log10(test_x + log_n0)
            test_x = (test_x - mean) / (std + q)
            lr.score(test_x, test_y)

            right_num_study += np.sum(test_y == lr.predict(test_x))

            pre_y = lr.predict(test_x)
            for m in range(len(test_x)):
                if test_y[m] == 1 and test_y[m] == pre_y[m]:
                    right_num_positive += 1

            test_data = pd.read_csv('pre_data//species_loso//{}//test_data.csv'.format(study), header=0, index_col=0)
            test_sample_id = test_data._stat_axis.values.tolist()
            test_loso_y = test_data['Label']
            test_loso_y = test_loso_y.values
            test_loso_x = test_data.drop(columns='Label')
            test_loso_x = test_loso_x.values
            test_loso_x = test_loso_x[:, std_preprocess != 0]
            test_loso_x = np.log10(test_loso_x + log_n0)
            test_loso_x = (test_loso_x - mean) / (std + q)
            svc = SVC(kernel='linear',C=parameters_loso[study],probability=True)
            svc.fit(train_x,train_y)
            predict_y = svc.predict(test_loso_x)
            num = 0
            for m in range(len(test_loso_x)):
                if test_loso_y[m] == 1 and predict_y[m] == 1:
                    num += 1
            recall = num / np.sum(test_loso_y == 1)
            proba_loso = svc.predict_proba(test_loso_x)[:, 1] / 100
            predict_matrix.loc[test_sample_id, 'LOSO'] += proba_loso
            recall_matrix.loc['LOSO',study] += recall / 100
            score = lr.score(test_loso_x, test_loso_y)
            accuracy_matrix.loc['LOSO', study] += score / 100
        accuracy_matrix.loc['LOSO_CV', study] += right_num_study / (len(x) * 10)
        recall_matrix.loc['LOSO_CV',study] += right_num_positive / (np.sum(y==1) * 10)

predict_matrix.to_csv('predict_matrix_svm.csv')
accuracy_matrix.to_csv('accuracy_matrix_svm.csv')
recall_matrix.to_csv('recall_matrix_svm.csv')


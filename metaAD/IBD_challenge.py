import pandas as pd
import numpy as np
import os
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.svm import SVC
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import accuracy_score,precision_score,recall_score

def data_preprocess(data, feature="function"):
    """
    :param data: row with samples and columns with features; relative abundance
    :param feature: function profile or species profile
    :return:
    """
    # convert to relative abundance
    data_rel = (data/data.sum(axis=0))
    data = data_rel
    if feature == "function":
        func_cutoff = 1e-6
        func_pesudo = 1e-9
        filter_data = data.loc[:, (data > func_cutoff).sum(axis=0) >= 3]
        feature_list = [i for i in range(data.shape[1]) if (data.values[:,i]>func_cutoff).sum(axis=0) >= 3 ]
        print("filtered featrues: " + str(data.shape[1] - filter_data.shape[1]))
        # log10-transformed
        log_data = np.log10(filter_data.astype(np.float64) + func_pesudo)
        # z-score normalization
        scale = StandardScaler()
        z_data = scale.fit_transform(log_data)

        return z_data,feature_list
    elif feature == "species":
        tax_cutoff = 1e-3
        tax_pesudo = 1e-5
        filter_data = data.loc[:, (data > tax_cutoff).sum(axis=0) >= 3]
        feature_list = [i for i in range(data.shape[1]) if (data.values[:,i]>tax_cutoff).sum(axis=0) >= 3 ]
        print("filtered featrues: " + str(data.shape[1] - filter_data.shape[1]))
        # log10-transformed
        log_data = np.log10(filter_data.astype(np.float64) + tax_pesudo)
        # z-score normalization
        scale = StandardScaler()
        z_data = scale.fit_transform(log_data)
        return z_data,feature_list


def train_cv(workdir, data, group, model, type='function',cv=10,repeat=10):
    """
    model:'svm','randomforest','lasso'
    :return:
    """

    # data_matrix = 'TrainingHe_TaxonomyAbundance_matrix.txt'
    # label = 'Class_labels_He.txt'
    data = pd.read_csv(os.path.join(workdir, data), header=0, index_col=0)
    data = data.loc[data['group'].isin(['NC', group])]
    feature_data = data.iloc[:, 3:]
    label = data.values[:,2]

    label =  np.array([1 if i == 'group' else -1 for i in label]).reshape(feature_data.shape[0],1)
    #print(label)
    feature_data,_ = data_preprocess(feature_data,feature=type)
    print('CrossValidation using {}'.format(model))
    if model == 'svm':

        # print(param)
        accuracy=[]
        precision=[]
        recall=[]
        for i in range(repeat):
            skf = StratifiedKFold(n_splits=cv)
            accuracy_result = np.zeros((len(label),1))
            for train_index, test_index in skf.split(feature_data, label):
                train_x = feature_data[train_index]
                train_y = label[train_index]
                test_x = feature_data[test_index]
                test_y = label[test_index]
                svm = SVC()
                param_grid = {'C': [0.0001, 0.001, 0.1, 1, 10, 100, 1000], 'gamma': [0.001, 0.0001],
                              'kernel': ['rbf', 'linear']}
                grid = GridSearchCV(svm, param_grid, cv=10, scoring='accuracy')
                grid.fit(train_x, train_y)
                param = grid.best_params_
                svm = SVC(kernel=param['kernel'],C=param['C'],gamma=param['gamma'])
                svm.fit(train_x,train_y)
                predict_y = svm.predict(test_x)
                #accuracy_result.append(accuracy_score(test_y,predict_y))
                accuracy_result[test_index] = np.array(predict_y).reshape(len(predict_y),1)
            #print(accuracy_result)
            accuracy.append(accuracy_score(label,accuracy_result))
            precision.append(precision_score(label,accuracy_result,average=None))
            recall.append(recall_score(label,accuracy_result,average=None))
        print("accuracy:" , np.array(accuracy).mean())
        print("precision:" , np.array(precision).mean(axis=0))
        print("recall:" , np.array(recall).mean(axis=0))

    if model == 'randomforest':

        #print(param)
        accuracy = []
        precision=[]
        recall=[]
        for i in range(repeat):
            skf = StratifiedKFold(n_splits=cv)
            accuracy_result = np.zeros((len(label),1))
            for train_index, test_index in skf.split(feature_data, label):
                train_x = feature_data[train_index]
                train_y = label[train_index]
                test_x = feature_data[test_index]
                test_y = label[test_index]
                rf = RandomForestClassifier(n_jobs=-1)
                param_grid = {'n_estimators': [10, 50, 100, 500, 1000]}
                grid = GridSearchCV(rf, param_grid, cv=10, scoring='accuracy')
                grid.fit(train_x, train_y)
                param = grid.best_params_
                rf = RandomForestClassifier(n_estimators =param['n_estimators'])
                rf.fit(train_x,train_y)
                predict_y = rf.predict(test_x)
                accuracy_result[test_index] = np.array(predict_y).reshape(len(predict_y),1)
                #accuracy_result.append(accuracy_score(test_y,predict_y))
            accuracy.append(accuracy_score(label,accuracy_result))
            precision.append(precision_score(label,accuracy_result,average=None))
            recall.append(recall_score(label,accuracy_result,average=None))
        print("accuracy:",np.array(accuracy).mean())
        print("precision:" ,np.array(precision).mean(axis=0))
        print("recall:" , np.array(recall).mean(axis=0))

    if model == 'lasso':

        #print(param)
        accuracy = []
        precision=[]
        recall=[]
        for i in range(repeat):
            skf = StratifiedKFold(n_splits=cv)
            accuracy_result = np.zeros((len(label),1))
            for train_index, test_index in skf.split(feature_data, label):
                train_x = feature_data[train_index]
                train_y = label[train_index]
                test_x = feature_data[test_index]
                test_y = label[test_index]
                lr = LogisticRegression(penalty='l1', solver='liblinear', n_jobs=-1)
                param_grid = {'C': [0.0001, 0.001, 0.1, 1, 10, 100, 1000]}
                grid = GridSearchCV(lr, param_grid, cv=10, scoring='accuracy')
                grid.fit(train_x, train_y)
                param = grid.best_params_
                lr = LogisticRegression(penalty='l1', solver='liblinear', n_jobs=-1,C=param['C'])
                lr.fit(train_x,train_y)
                predict_y = lr.predict(test_x)
                accuracy_result[test_index] = np.array(predict_y).reshape(len(predict_y),1)
                #accuracy_result.append(accuracy_score(test_y,predict_y))
            accuracy.append(accuracy_score(label,accuracy_result))
            precision.append(precision_score(label,accuracy_result,average=None))
            recall.append(recall_score(label,accuracy_result,average=None))
        print("accuracy:",np.array(accuracy).mean())
        print("precision:" ,np.array(precision).mean(axis=0))
        print("recall:" , np.array(recall).mean(axis=0))


def train_transfer(train_data,train_label,test_data,test_label,repeat=10,type='function',model='lasso',is_TCA=True):

    train_feature_data = pd.read_csv('data/testset_subchallenge2_files/' + train_data, sep='\t').T[1:]
    train_label = pd.read_csv('data/testset_subchallenge2_files/' + train_label, sep='\t').values[:, 2]
    train_label = np.array([1 if i != 'nonIBD' else -1 for i in train_label]).reshape(len(train_feature_data), 1)

    test_feature_data = pd.read_csv('data/testset_subchallenge2_files/' + test_data, sep='\t').T[1:]
    test_label = pd.read_csv('data/testset_subchallenge2_files/' + test_label, sep='\t').values[:, 2]
    test_label = np.array([1 if i != 'nonIBD' else -1 for i in test_label]).reshape(len(test_feature_data), 1)

    train_feature_data,train_feature_list = data_preprocess(train_feature_data, feature=type)
    test_feature_data = np.array(test_feature_data)[:,train_feature_list]
    if type == 'function':
        func_cutoff = 1e-6
        func_pesudo = 1e-9
        # log10-transformed
        log_data = np.log10(test_feature_data.astype(np.float64) + func_pesudo)
        # z-score normalization
        scale = StandardScaler()
        test_feature_data = scale.fit_transform(log_data)
    else:
        tax_cutoff = 1e-3
        tax_pesudo = 1e-5
        # log10-transformed
        log_data = np.log10(test_feature_data.astype(np.float64) + tax_pesudo)
        # z-score normalization
        scale = StandardScaler()
        test_feature_data = scale.fit_transform(log_data)

    if is_TCA == True:
        from TCA import TCA
        train_feature_data, test_feature_data = TCA(dim=150, kernel='linear', gamma=1, lamb=1).fit(train_feature_data.astype(np.float64)
                                                                                                   ,test_feature_data.astype(np.float64))
        scale = StandardScaler()
        # train_feature_data = scale.fit_transform(train_feature_data)
        # test_feature_data = scale.fit_transform(test_feature_data)
    print('Transfer training using {}'.format(model))
    if model == 'svm':
        svm = SVC(class_weight='balanced')
        param_grid = {'C': [0.0001, 0.001, 0.1, 1, 10, 100, 1000], 'gamma': [0.001, 0.0001],
                      'kernel': ['rbf', 'linear']}
        grid = GridSearchCV(svm, param_grid, cv=10, scoring='accuracy')
        grid.fit(train_feature_data, train_label)
        param = grid.best_params_
        print(param)
        accuracy = []
        precision = []
        recall = []
        for i in range(repeat):
            svm = SVC(kernel=param['kernel'], C=param['C'], gamma=param['gamma'],class_weight='balanced')
            svm.fit(train_feature_data, train_label)
            predict_y = svm.predict(test_feature_data)
            # accuracy_result.append(accuracy_score(test_y,predict_y))
            # print(accuracy_result)
            accuracy.append(accuracy_score(test_label, predict_y))
            precision.append(precision_score(test_label, predict_y, average=None))
            recall.append(recall_score(test_label, predict_y, average=None))
        print("accuracy:", np.array(accuracy).mean())
        print("precision:", np.array(precision).mean(axis=0))
        print("recall:", np.array(recall).mean(axis=0))

    if model == 'randomforest':
        rf = RandomForestClassifier(n_jobs=-1,class_weight='balanced')
        param_grid = {'n_estimators': [10, 50, 100,500,1000]}
        grid = GridSearchCV(rf, param_grid, cv=10, scoring='accuracy')
        grid.fit(train_feature_data,train_label)
        param = grid.best_params_
        print(param)
        accuracy = []
        precision=[]
        recall=[]
        for i in range(repeat):
            rf = RandomForestClassifier(n_estimators =param['n_estimators'],class_weight='balanced')
            rf.fit(train_feature_data,train_label)
            predict_y = rf.predict(test_feature_data)
            accuracy.append(accuracy_score(test_label,predict_y))
            precision.append(precision_score(test_label,predict_y,average=None))
            recall.append(recall_score(test_label,predict_y,average=None))
        print("accuracy:",np.array(accuracy).mean())
        print("precision:" ,np.array(precision).mean(axis=0))
        print("recall:" , np.array(recall).mean(axis=0))

    if model == 'lasso':
        lr = LogisticRegression(penalty='l1', solver='liblinear', n_jobs=-1,class_weight='balanced')
        param_grid = {'C':[ 0.0001,0.001, 0.1, 1, 10, 100, 1000]}
        grid = GridSearchCV(lr, param_grid, cv=10, scoring='accuracy')
        grid.fit(train_feature_data,train_label)
        param = grid.best_params_
        print(param)
        accuracy = []
        precision=[]
        recall=[]
        for i in range(repeat):
            lr = LogisticRegression(penalty='l1', solver='liblinear', n_jobs=-1,C=param['C'],class_weight='balanced')
            lr.fit(train_feature_data,train_label)
            predict_y = lr.predict(test_feature_data)
            accuracy.append(accuracy_score(test_label,predict_y))
            precision.append(precision_score(test_label,predict_y,average=None))
            recall.append(recall_score(test_label,predict_y,average=None))
        print("accuracy:",np.array(accuracy).mean())
        print("precision:" ,np.array(precision).mean(axis=0))
        print("recall:" , np.array(recall).mean(axis=0))




if __name__ == '__main__':
    import warnings
    warnings.filterwarnings("ignore")
    #train_cv('TrainingSchirmer_PathwayAbundance_matrix.txt','Class_labels_Schirmer.txt',model='svm',type='function')
    # print('\n\n')
    workdir = "F:\\Zhaolab2020\\gut-brain-axis\\metaAD\\01crc_meta\\meta-analysis-microbiome\\AD_species"
    train_cv(workdir, 'CHN.csv', group='AD', model='lasso',type='species')

    # print('\n\n')

    # train_cv('TrainingSchirmer_TaxonomyAbundance_matrix.txt', 'Class_labels_Schirmer.txt', model='lasso', type='species')

    # train_cv('TrainingSchirmer_TaxonomyAbundance_matrix.txt', 'Class_labels_Schirmer.txt', model='lasso', type='function')

    # train_transfer('TrainingHe_TaxonomyAbundance_matrix.txt','Class_labels_He.txt',
    #                'TrainingSchirmer_TaxonomyAbundance_matrix.txt','Class_labels_Schirmer.txt',type='species',model='svm')
    #
    #
    # train_transfer('TrainingSchirmer_TaxonomyAbundance_matrix.txt','Class_labels_Schirmer.txt',
    #              'TrainingHe_TaxonomyAbundance_matrix.txt','Class_labels_He.txt',type='species',model='svm')


import matplotlib.pyplot as plt
import numpy as np
import os

from sklearn.datasets import make_classification
from sklearn.model_selection import KFold
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_curve

def plot_kfold_auroc(X, y, img):
    kf = KFold(n_splits=10)

    tprs = []
    base_fpr = np.linspace(0, 1, 101)

    plt.figure(figsize=(5, 5))
    plt.axes().set_aspect('equal', 'datalim')

    for i, (train, test) in enumerate(kf.split(X)):
        model = LogisticRegression().fit(X[train], y[train])
        y_score = model.predict_proba(X[test])
        fpr, tpr, _ = roc_curve(y[test], y_score[:, 1])

        plt.plot(fpr, tpr, 'b', alpha=0.15)
        tpr = np.interp(base_fpr, fpr, tpr)
        tpr[0] = 0.0
        tprs.append(tpr)

    tprs = np.array(tprs)
    mean_tprs = tprs.mean(axis=0)
    std = tprs.std(axis=0)

    tprs_upper = np.minimum(mean_tprs + std, 1)
    tprs_lower = mean_tprs - std

    plt.plot(base_fpr, mean_tprs, 'b')
    plt.fill_between(base_fpr, tprs_lower, tprs_upper, color='grey', alpha=0.3)

    plt.plot([0, 1], [0, 1], 'r--')
    plt.xlim([-0.01, 1.01])
    plt.ylim([-0.01, 1.01])
    plt.ylabel('True Positive Rate')
    plt.xlabel('False Positive Rate')
    plt.show()
    plt.savefig(img)

def plot_repeat_kfold_auroc(X,y, img):
    idx = np.arange(0, len(y))

    tprs = []
    base_fpr = np.linspace(0, 1, 101)

    plt.figure(figsize=(5, 5))
    plt.axes().set_aspect('equal', 'datalim')

    for j in np.random.randint(0, high=10000, size=10):
        print(j)
        np.random.shuffle(idx)
        kf = KFold(n_splits=10, shuffle=True, random_state=j)

        for i, (train, test) in enumerate(kf.split(X)):
            model = LogisticRegression().fit(X[idx][train], y[idx][train])
            y_score = model.predict_proba(X[idx][test])
            fpr, tpr, _ = roc_curve(y[idx][test], y_score[:, 1])

            plt.plot(fpr, tpr, 'b', alpha=0.05)
            tpr = np.interp(base_fpr, fpr, tpr)
            tpr[0] = 0.0
            tprs.append(tpr)

    tprs = np.array(tprs)
    mean_tprs = tprs.mean(axis=0)
    std = tprs.std(axis=0)

    tprs_upper = np.minimum(mean_tprs + std, 1)
    tprs_lower = mean_tprs - std

    plt.plot(base_fpr, mean_tprs, 'b')
    plt.fill_between(base_fpr, tprs_lower, tprs_upper, color='grey', alpha=0.3)

    plt.plot([0, 1], [0, 1], 'r--')
    plt.xlim([-0.01, 1.01])
    plt.ylim([-0.01, 1.01])
    plt.ylabel('True Positive Rate')
    plt.xlabel('False Positive Rate')
    plt.show()
    plt.savefig(img)

def main():
    X, y = make_classification(n_samples=500, random_state=100, flip_y=0.3)
    workdir = "D:\\Zhaolab2020\\gut-brain-axis\\metaAD\\01MLmodels\\multikingdom\\testwork"
    img1 = os.path.join(workdir, "plot_kfold_auroc.png")
    img2 = os.path.join(workdir, "plot_repeat_kfold_auroc.png")
    plot_kfold_auroc(X, y, img1)
    plot_repeat_kfold_auroc(X, y, img2)

if __name__ == '__main__':
    main()

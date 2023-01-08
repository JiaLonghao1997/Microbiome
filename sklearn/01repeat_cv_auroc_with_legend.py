import matplotlib.pyplot as plt
import numpy as np
import os
from sklearn.datasets import load_iris
from sklearn.datasets import make_classification
from sklearn.model_selection import KFold, StratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn import svm
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import RocCurveDisplay

def plot_kfold_auroc(X, y, target_names, clf, img):
    #参考: https://scikit-learn.org/stable/auto_examples/model_selection/plot_roc_crossval.html

    tprs = []
    aucs = []
    mean_fpr = np.linspace(0, 1, 101)

    ##fig, ax = plt.figure(figsize=(5, 5))
    ##Figure是一个大画板，Axes是x,y坐标轴等。
    fig, ax = plt.subplots(figsize=(6, 6))  ##nrows, ncols默认是1.
    #plt.axes().set_aspect('equal', 'datalim')
    cv = StratifiedKFold(n_splits=5)

    for fold, (train, test) in enumerate(cv.split(X, y)):
        clf = clf.fit(X[train], y[train])
        y_score = clf.predict_proba(X[test])
        fpr, tpr, _ = roc_curve(y[test], y_score[:, 1])
        roc_auc = auc(fpr, tpr)

        ##
        interp_tpr = np.interp(mean_fpr, fpr, tpr)
        interp_tpr[0] = 0.0
        tprs.append(interp_tpr)
        aucs.append(roc_auc)

    ax.plot([0, 1], [0, 1], "k--", label="chance level (AUC = 0.5)")
    mean_tpr = np.mean(tprs, axis=0)
    mean_tpr[-1] = 1.0
    mean_auc = auc(mean_fpr, mean_tpr)
    std_auc = np.std(aucs)
    ax.plot(
        mean_fpr,
        mean_tpr,
        color="b",
        label=r"Mean ROC (AUC = %0.2f $\pm$ %0.2f)" % (mean_auc, std_auc),
        lw=2,
        alpha=0.8,
    )

    ax.set(
        xlim=[-0.05, 1.05],
        ylim=[-0.05, 1.05],
        xlabel="False Positive Rate",
        ylabel="True Positive Rate",
        title=f"Mean ROC curve with variability\n(Positive label '{target_names[1]}')",
    )
    ax.axis("square")
    ax.legend(loc="lower right")
    plt.show()
    fig.savefig(img)

def plot_repeat_kfold_auroc(X,y, target_names, clf, img):
    idx = np.arange(0, len(y))

    tprs = []
    aucs = []
    mean_fpr = np.linspace(0, 1, 101)

    fig, ax = plt.subplots(figsize=(6, 6))  ##nrows, ncols默认是1.

    for j in np.random.randint(0, high=10000, size=10):
        print(j)
        np.random.shuffle(idx)
        #kf = KFold(n_splits=10, shuffle=True, random_state=j)
        cv = StratifiedKFold(n_splits=6, shuffle=True, random_state=j)

        for i, (train, test) in enumerate(cv.split(X, y)):
            #model = LogisticRegression().fit(X[idx][train], y[idx][train])
            clf = clf.fit(X[train], y[train])
            y_score = clf.predict_proba(X[idx][test])
            fpr, tpr, _ = roc_curve(y[idx][test], y_score[:, 1])
            roc_auc = auc(fpr, tpr)

            ##
            interp_tpr = np.interp(mean_fpr, fpr, tpr)
            interp_tpr[0] = 0.0
            tprs.append(interp_tpr)
            aucs.append(roc_auc)

    ax.plot([0, 1], [0, 1], "k--", label="chance level (AUC = 0.5)")
    mean_tpr = np.mean(tprs, axis=0)
    mean_tpr[-1] = 1.0
    mean_auc = auc(mean_fpr, mean_tpr)
    std_auc = np.std(aucs)
    ax.plot(
        mean_fpr,
        mean_tpr,
        color="b",
        label=r"Mean ROC (AUC = %0.2f $\pm$ %0.2f)" % (mean_auc, std_auc),
        lw=2,
        alpha=0.8,
    )

    # std_tpr = np.std(tprs, axis=0)
    # tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
    # tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
    # ax.fill_between(
    #     mean_fpr,
    #     tprs_lower,
    #     tprs_upper,
    #     color="grey",
    #     alpha=0.2,
    #     label=r"$\pm$ 1 std. dev.",
    # )

    ax.set(
        xlim=[-0.05, 1.05],
        ylim=[-0.05, 1.05],
        xlabel="False Positive Rate",
        ylabel="True Positive Rate",
        title=f"Mean ROC curve with variability\n(Positive label '{target_names[1]}')",
    )
    ax.axis("square")
    ax.legend(loc="lower right")
    plt.show()
    fig.savefig(img)

def main():
    #X, y = make_classification(n_samples=500, random_state=100, flip_y=0.3)
    iris = load_iris()
    target_names = iris.target_names
    X, y = iris.data, iris.target
    X, y = X[y != 2], y[y != 2]
    n_samples, n_features = X.shape

    random_state = np.random.RandomState(0)
    X = np.concatenate([X, random_state.randn(n_samples, 200 * n_features)], axis=1)

    ##
    clf = svm.SVC(kernel="linear", probability=True, random_state=random_state)

    workdir = "D:\\Zhaolab2020\\gut-brain-axis\\metaAD\\01MLmodels\\multikingdom\\testwork"
    img1 = os.path.join(workdir, "plot_kfold_auroc.png")
    img2 = os.path.join(workdir, "plot_repeat_kfold_auroc.png")
    plot_kfold_auroc(X, y, target_names, clf, img1)
    plot_repeat_kfold_auroc(X, y, target_names, clf, img2)

if __name__ == '__main__':
    main()
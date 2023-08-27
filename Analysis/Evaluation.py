def plot_hist(hist,title):
#     plt.plot(hist.history["accuracy"])
#     plt.plot(hist.history["val_accuracy"])
    plt.plot(hist.history["loss"])
    plt.plot(hist.history["sparse_categorical_accuracy"])
    plt.plot(hist.history["val_loss"])
    plt.plot(hist.history["val_sparse_categorical_accuracy"])
#     plt.title("Training Progress")
    plt.ylabel("Accuracy/Loss")
    plt.xlabel("Epochs")
    plt.legend(["loss","sparse_categorical_accuracy","val_loss","val_sparse_categorical_accuracy"])
    #plt.show()
    plt.savefig('mi/{}.pdf'.format(title), bbox_inches='tight')
    plt.clf()





def plot_auc(y_train_le,y_score,title):
    y_test = label_binarize(y_train_le, classes=[0,1, 2, 3, 4])
    n_classes = 5
    fpr = dict()
    tpr = dict()
    roc_auc = dict()
    for i in range(n_classes):
        fpr[i], tpr[i], _ = roc_curve(y_test[:, i], y_score[:, i])
        roc_auc[i] = auc(fpr[i], tpr[i])

    # Compute micro-average ROC curve and ROC area
    fpr["micro"], tpr["micro"], _ = roc_curve(y_test.ravel(), y_score.ravel())
    roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])

    # Compute macro-average ROC curve and ROC area
    # First aggregate all false positive rates
    all_fpr = np.unique(np.concatenate([fpr[i] for i in range(n_classes)]))
    # Then interpolate all ROC curves at this points
    mean_tpr = np.zeros_like(all_fpr)
    for i in range(n_classes):
        mean_tpr += interp(all_fpr, fpr[i], tpr[i])
    # Finally average it and compute AUC
    mean_tpr /= n_classes
    fpr["macro"] = all_fpr
    tpr["macro"] = mean_tpr
    roc_auc["macro"] = auc(fpr["macro"], tpr["macro"])

    # Plot all ROC curves
    lw = 3
    plt.figure()

    plt.plot(fpr["micro"], tpr["micro"],
             label='micro-average ROC curve (area = {0:0.3f})'
                   ''.format(roc_auc["micro"]),
             color='#DEB887', linestyle=':', linewidth=2)

    plt.plot(fpr["macro"], tpr["macro"],
             label='macro-average ROC curve (area = {0:0.3f})'
                   ''.format(roc_auc["macro"]),
             color='#8B8B7A', linestyle=':', linewidth=2)

    ####################################
    colors = cycle(['#D53E4F', '#F46D43', '#B1DDA4', '#66C2A5', '#3288BD'])
    for i, color in zip(range(n_classes), colors):
        plt.plot(fpr[i], tpr[i], color=color, lw=lw,
                 label='ROC curve of class {0} (area = {1:0.3f})'
                       ''.format(i, roc_auc[i]))
    plt.plot([0, 1], [0, 1], 'k--', lw=lw)
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    # plt.title('Some extension of Receiver operating characteristic to multi-class')
    plt.title=title
    plt.legend(loc="lower right")
    plt.savefig('mi/{}.pdf'.format(title), bbox_inches='tight')
    plt.clf()


def plot_confusion_matrix(o,cm, classes, normalize=False, title='Confusion matrix', cmap=plt.cm.Blues):
    if normalize:
        cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
        print("Normalized confusion matrix")
    else:
        print('Confusion matrix, without normalization')
    plt.imshow(cm, interpolation='nearest', cmap=cmap)
    plt.colorbar()
    tick_marks = np.arange(len(classes))
    plt.xticks(tick_marks, classes, rotation=45)
    plt.yticks(tick_marks, classes)
    fmt = '.2f' if normalize else 'd'
    thresh = cm.max() / 2.
    for i, j in itertools.product(range(cm.shape[0]), range(cm.shape[1])):
        plt.text(j, i, format(cm[i, j], fmt),
                 horizontalalignment="center",
                 color="white" if cm[i, j] > thresh else "black")
    plt.tight_layout()
    plt.ylabel('True label')
    plt.xlabel('Predicted label')
    plt.savefig('mi/{}.pdf'.format(o), bbox_inches='tight')
    plt.clf()

def evaluation( y, y_train_pre, y_test,y_pred_model,o):

    train_label=pd.DataFrame(y,columns=['y'])
    train_label['y_train_pre']=y_train_pre
    
    test_label=pd.DataFrame(y_test,columns=['y_test'])
    test_label['y_test_pre']=y_pred_model
        
    print('\n>>>evaluation---wait...')
    score_train = metrics.accuracy_score(y, y_train_pre)
    score_test = metrics.accuracy_score(y_test, y_pred_model)
    
    matrix = metrics.confusion_matrix(y_test, y_pred_model)
    report = metrics.classification_report(y_test, y_pred_model)
    
    
    print('>>>train_acc\n', score_train)
    print('>>>test_acc\n', score_test)
    print('\n>>>confusion_matrix\n', matrix)
    labels_name = ['1','2','3','4','5']
    plot_confusion_matrix(o,matrix,labels_name,"confusion_matrix")
    print('\n>>>reacall\n', report)
    print('>>>end end ...')

import numpy as np
from sklearn.cluster import DBSCAN
from sklearn.svm import SVC
from sklearn.naive_bayes import GaussianNB

def main(df):

  if (df["dcode"]==2).sum() < 2:
    import sys
    msg = "Insufficient pathogenic variation for DBSCAN clustering.\n"
    sys.stderr.write(msg)
    return df

  print "\nTotal residues:     %d"%len(df)
  print "\nTotal pathogenic:   %d"%(df["dcode"]==2).sum()

  # Identify clusters of pathogenic variants (class IDs >0)
  # Pathogenic "noise" is treated like benign variation
  df.ix[df["dcode"]==2,"cluster"] = \
                    pathclust(df.ix[df["dcode"]==2,["x","y","z"]].values)
  # Label known or putative benign as class -1
  df.ix[df["dcode"].isin([0,1]),"cluster"] = -1
  print "\nDBSCAN Cluster IDs and Counts:"
  print "-1: Putative Benign"
  print " 0: Pathogenic 'Noise'"
  print ">0: Pathogenic Clusters"
  print df.groupby("cluster").apply(len)
  print "\nDBSCAN Cluster ID Predictions and Count:"
  print df.ix[df["dcode"]==2,["seqid","pos","dclass","dcode","cluster"]]

  # Subset the dataframe to training samples
  tdf = df.ix[(~df["cluster"].isnull()) & (df["cluster"]!=0)].copy()

  # Plot the decision space
  def assign_color(v):
    if int(v)==-1: 
      return "darkblue"
    elif int(v)==0: 
      return "darkred"
    elif v>0:
      return ["orange","yellow","purple"][int(v)%3]
  from mpl_toolkits.mplot3d import Axes3D
  import matplotlib.pyplot as plt
  print "\n%d training values...\n"%len(tdf)
  # Setup the figure and axes
  fig = plt.figure(figsize=(40,40))
  ax = fig.add_subplot(111,projection='3d')
  ax.set_xlabel("X Coordinate")
  ax.set_ylabel("Y Coordinate")
  ax.set_zlabel("Z Coordinate")
  ax.set_xticks(())
  ax.set_yticks(())
  plt.title("SVM Prediction Space")
  tdf["colors"] = tdf["cluster"].apply(assign_color)
  ax.scatter(tdf["x"].values,
             tdf["y"].values,
             tdf["z"].values,
             # color=tdf["colors"].values,
             c=tdf["cluster"].values,
             cmap=plt.cm.prism,
             marker='^',
             s=100)

  # C set to 0.01 for negatives, and 1.0 for positives
  w = dict((i,10) for i in tdf["cluster"].drop_duplicates().values)
  svm = SVC(C=0.5,gamma=0.01,kernel="rbf",probability=True,
              class_weight=w,decision_function_shape='ovr')
  svm = svm.fit(tdf[["x","y","z"]].values,tdf["cluster"].values)
  s   = svm.score(tdf[["x","y","z"]].values,tdf["cluster"].values)
  # svm = GaussianNB().fit(tdf[["x","y","z"]].values,tdf["cluster"].values)
  # Predict class for non-training set residues
  # pdf = df.copy()
  pdf = df[(df["cluster"].isnull()) | (df["cluster"]==0)].copy()
  pdf["cluster"] = svm.predict(pdf[["x","y","z"]].values)
  print "\nC=0.5; G=0.01; W=10; Mean Accuracy=%.3f"%s
  for name,pdfg in pdf.groupby("cluster"):
    print " % d: %4d"%(name,len(pdfg))
  # Probability of in a pathogenic cluster
  # pdf["probability"] = svm.predict_proba(pdf[["x","y","z"]].values)[:,1]
  pdf["probability"] = np.nanmax(svm.predict_proba(pdf[["x","y","z"]].values)[:,1:],axis=1)
  pdf["colors"]      = pdf["cluster"].apply(assign_color)
  ax.scatter(pdf["x"].values,
             pdf["y"].values,
             pdf["z"].values,
             # color=pdf["colors"].values,
             c=pdf["cluster"].values,
             cmap=plt.cm.prism,
             marker='o',
             s=40)
  # Save the figure to PDF
  plt.savefig("SVM_Cluster_Decision_Function.pdf")
  plt.close(fig)

  # LOO Train/Test for training data
  for idx in tdf.index:
    tempdf = tdf.drop([idx])
    # Retrain DBSCAN without this variant
    tempdf.ix[tempdf["dcode"]==2,"cluster"] = \
              pathclust(tempdf.ix[tempdf["dcode"]==2,["x","y","z"]].values)
    tempdf.ix[tempdf["dcode"].isin([0,1]),"cluster"] = -1
    # Retrain the SVM with the DBSCAN result
    tempdf = tempdf[(~tempdf["cluster"].isnull()) & (tempdf["cluster"]!=0)]
    svm = svm.fit(tempdf[["x","y","z"]].values,tempdf["cluster"].values)
    # Predict variant class
    probs = svm.predict_proba(tdf.ix[idx,["x","y","z"]].values.reshape(1,-1))[0,1:]
    if probs.ndim<2:
      probs = probs.reshape(probs.shape[0],1)
    tdf.ix[idx,"probability"] = np.nanmax(probs,axis=1)

  # ndf = df.ix[df["cluster"]==0].copy()
  # ndf["probability"] = 0.0
  # Recombine training, test, and noise
  df = tdf.append(pdf)
  # df = pdf

  print df.groupby("cluster").apply(len)

  with open("prediction.attr","wb") as fout:
    fout.write("#model\tseqid\tchain\tpathogenic\n")
    fout.write("attribute: prediction\n")
    fout.write("match mode: 1-to-1\n")
    fout.write("recipient: residues\n")
    for _,r in df.iterrows():
      fout.write("\t:%d\t%d\n"%(r["seqid"],r["cluster"]))
  print ""
  with open("probability.attr","wb") as fout:
    fout.write("#model\tseqid\tchain\tprobability\n")
    fout.write("attribute: probability\n")
    fout.write("match mode: 1-to-1\n")
    fout.write("recipient: residues\n")
    for _,r in df.iterrows():
      fout.write("\t:%d\t%.4f\n"%(r["seqid"],r["probability"]))
  return df,svm

def pathclust(v):
  """ v is a numpy matrix of pathogenic feature vectors """
  db = DBSCAN(eps=20.,min_samples=3.).fit(v)
  return db.labels_.astype(int)+1
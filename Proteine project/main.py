import pandas as pd
from pandas import *
from Bio import pairwise2
import io
import Bio
from Bio import Align
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
import difflib
from Bio.pairwise2 import format_alignment

df = pd.read_csv(r'C:\Users\hsain\PycharmProjects\Proteine project\Datasetprotein.csv', low_memory=False, )

df = df[df["macromoleculeType"].str.contains("Protein", "Protein#DNA") == True]

df.fillna(0)
# print(df.describe())
# print(df)

df1 = pd.read_csv(r'C:\Users\hsain\PycharmProjects\Proteine project\Reference sequence.csv', low_memory=False, )
#print(df1)

X1 = df['sequence'].tolist()
X2 = df['structureId'].tolist()
#del X1[-1]
# print(X1)
control = "GHMEQTHRAIFRFVPRHEDELELEVDDPLLLELQAEDYWYEAYNMRTGARGVFPAYYAIEVTK"
# print(X2)

sequences = []

for item in X1:
    i = str(item)
    if i.isalpha() and len(i) < 100:
        sequences.append(item)
        continue
# print difflib.get_close_matches(
# print(duplicated), [X1])

# alignments = pairwise2.align.globalxx("RAIFRFVPRHE",  "NGEEHEQTHRAIFRFVPRHEDELELEVDDPLLVELQAEDYWYEAYNMRTGARGVFPAYYAIEVTKEPEHMA")
# from Bio.pairwise2 import format_alignment
# print(format_alignment(*alignments[0]))


aligner = Align.PairwiseAligner()

r = []

del sequences[-1]


for sequence in sequences:
    #print(sequence)
    m = len(control) * 0.6
    score = aligner.align(control, str(sequence))

    if score.score > m:
        #a = pairwise2.align.globalxx(control, sequence)
        a = pairwise2.align.globalms(control, sequence, 2, -1, -.5, -.1)
        b = format_alignment(*a[0])
        #print(type(b))
        r.append(b)
        #print(score.score)
        print(format_alignment(*a[0]))

for item in r:
    print(item)


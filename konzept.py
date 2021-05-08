import pandas as pd
import numpy as np
import os
from msa import msa

result_file = "result.csv"

lm_eps = {"AN":36.64, "DCM": 8.93}
sigma = [4.94, 1.26]
c0 = [1, 0.002, 0.00007]
cs = []
base = []

for c in c0:
    expo = np.floor(np.log10(c))-2
    x = 10 ** (expo)
    a = np.linspace(c-x, c+x,num=21, endpoint=True)
    a = np.round(a, int(2-expo))
    cs.extend(a)
    for val in a:
        base.append([c, val])
base = pd.DataFrame(base, columns=["base", "c"]).set_index("c")

values=[]
for (lm,eps) in lm_eps.items():
    for rhomol in cs:
        wert = msa(rhomol, eps, sigma)
        values.append([lm, rhomol, wert])

test = pd.DataFrame(values, columns=["Lösungsmittel", "c", "Ergebnis"])

spalten = [base]
for lm in test["Lösungsmittel"].unique():
    an = test[test["Lösungsmittel"]==lm]
    an = an[["c", "Ergebnis"]] # Beschränkt die Spalten
    an = an.set_index("c")
    an.columns=[lm] # Benennt die Spalte Ergebnis um
    spalten.append(an)

result = pd.concat(spalten, axis=1)
print(result)
result.to_csv(result_file, sep=";", decimal=",") #, float_format="%.10E")

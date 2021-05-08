import pandas as pd
import numpy as np
import os
from msa import msa

result_file = "result.csv"

lm_eps = {"AN":36.64, "DCM": 8.93}
sigma = [4.94, 1.26]
c0 = [1, 0.002, 0.00007]
cs = []

for c in c0:
    expo = np.floor(np.log10(c))-2
    x = 10 ** (expo)
    a = np.linspace(c-x, c+x,num=21, endpoint=True)
    a = np.round(a, int(2-expo))
    cs.extend(a)


values=[]
for (lm,eps) in lm_eps.items():
    for rhomol in cs:
        wert = msa(rhomol, eps, sigma)
        values.append([lm, rhomol, wert])

# values=[]
# # Gehe rekursiv durch alle Ordner in "experimente"
# for lm in os.scandir("experimente"):
#     if lm.is_dir():
#         for c in os.scandir(lm.path):
#             if c.is_dir:
#                 # Öffne die Ergebnis-Datei und schreibe die Werte in values
#                 file = open(c.path+"/values.8", 'r') # values.8 noch anpassen!
#                 val = float(list(file)[3]) # richtige Zeile (ab 0 zählen)
#                 values.append([lm.name, c.name, val])

test = pd.DataFrame(values, columns=["Lösungsmittel", "c", "Ergebnis"])
print(test)

spalten = []
for lm in test["Lösungsmittel"].unique():
    an = test[test["Lösungsmittel"]==lm]
    an = an[["c", "Ergebnis"]] # Beschränkt die Spalten
    an = an.set_index("c")
    an.columns=[lm] # Benennt die Spalte Ergebnis um
    spalten.append(an)

result = pd.concat(spalten, axis=1)
print(result)
result.to_csv(result_file, sep=";", decimal=",") #, float_format="%.10E")

# -*- coding: utf-8 -*-
""" 
Created on Sat May  1 09:58:59 2021

@author: Timo
"""
import os


Salz = input('agpf? (Y/n)')

if Salz == 'n':
    print('Du musst den Radius in der Datei test.2 an dein Salz anpassen')           
 
else:
    LM = input('Lösugsmittel?')
    epsilonr = input('Permittivität des Lösungsmittel?')
    c = input ('Konzentration? (Format: #.###...d-#)')
    
    if os.path.exists(str(LM))==True:
        os.chdir(str(LM))
    else:
        os.system('mkdir '+str(LM))
        os.chdir(str(LM))

    if os.path.exists(str(c)+"mM") == True: #Wenn Ordner schon existiert wurde die Rechnung schon durchgeführt. Ebene zurück gehen und Programm beenden. sonst Ordner anlegen und zurück..
        print("Rechnung schon mal durchgeführt.")
        os.chdir('..')
        exit()

    else:
        os.system('mkdir '+str(c)+"mM")
        os.chdir('..')
    
    os.system('cp test2.f '+str(LM)+'/'+str(c)+'mM')
    os.system('cp msa '+str(LM)+'/'+str(c)+'mM')
    os.chdir(str(LM)+'/'+str(c)+'mM')
    os.system('mv test2.f agpf'+str(LM)+str(c)+'.f')
    
    os.system('gfortran( -o msa '+ str(LM)+str(c)+'.f')
    os.system("./msa "+ str(c) +" " + str(epsilonr))

    Ergebnis=open('fort.8','r')
    Ergebnis.readline()
    Ergebnis.readline()
    Ergebnis.readline()
    Ergebnis.readline()
    Ergebnis.readline()
    mu=str(Ergebnis.readline())
    mu=mu.strip()
    mu=mu.replace(".",",")
    Ergebnis.close()
    print(mu)
    os.chdir('..')
    os.chdir('..')
    
   
    os.system("echo '" + str(LM) + ";" + str(c) + ";" + mu + "\n' >> ausgabe.csv")
    
   # file = open("ausgabe.csv", "a")
   # file.write(str(mu))
   # file.close()
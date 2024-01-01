
import pandas as pd
from tkinter.filedialog import askopenfilename
#from selenium import webdriver
#import chromedriver_binary
#import time
from tkinter import Tk,Text,mainloop,END,Button
from selenium.webdriver.firefox.service import Service
from selenium.webdriver.firefox.options import Options
from selenium.webdriver.firefox.firefox_binary import FirefoxBinary
import time
import os
from selenium import webdriver

options=Options()

torpath="C:\\Tor Browser\\Browser\\TorBrowser\\Tor\\tor.exe"
torprofpath="C:\\Tor Browser\\Browser\\TorBrowser\\Data\\Browser\\profile.default"
firepath="C:\\Tor Browser\\Browser\\firefox.exe"
geckopath="C:\\WebDrivers\\geckodriver.exe"
options.set_preference('network.proxy.type', 1)
options.set_preference('network.proxy.socks', '127.0.0.1')
options.set_preference('network.proxy.socks_port', 9050)
options.set_preference("network.proxy.socks_remote_dns", False)
options.set_preference('profile', torprofpath)
os.popen(torpath)
binary = FirefoxBinary(firepath)
service = Service(geckopath)
driver = webdriver.Firefox(service=service,options=options,firefox_binary=binary)



dfile="C:\\Users\\domin\\sciebo2\\Promotion\\Wissen\\Doktorarbeit\\Quellen_zum_Checken.xlsx"
data=pd.read_excel(dfile)

Number=data["Number"]
References=data["References"]


#binary = FirefoxBinary(firepath)
#service = Service(geckopath)
#driver = webdriver.Chrome("C:\\Program Files (x86)\\Google\\Chrome\\Application\\chrome.exe")#(service=service,options=options,firefox_binary=binary)
#driver = webdriver.Chrome()
url= "https://scholar.google.com/scholar?hl=de&as_sdt=0%2C5&q="
# enter the publication in google scholarb
# geb leertaste zu Prozent 20

   
count=0
for i,val in enumerate (References):
    count+=1
    driver.get(url+str(References[i]))
    root = Tk()
    #root.geometry('%dx%d+%d+%d' % (500,400,200,200))
    #T = Text(root)
    def GotYou():
        with open("Target"+str(i)+".txt", "w") as f:
            f.write(References[i])
    T = Button(root, text=References[i], command=GotYou)
    T.pack()
    T.focus_force()
    mainloop()
    
    if count>5:
        driver.close()
        #time.sleep(5)
        driver = webdriver.Firefox(service=service,options=options,firefox_binary=binary)
        count=0
#[driver(url+str(val)) for i,val in enumerate (References)]

#Show text box with publication




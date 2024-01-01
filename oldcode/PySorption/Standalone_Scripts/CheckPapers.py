



from selenium import webdriver
from selenium.webdriver.firefox.service import Service
from selenium.webdriver.firefox.options import Options
from selenium.webdriver.firefox.firefox_binary import FirefoxBinary
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.action_chains import ActionChains
from numpy.random import rand
import time
import os
def CallPapers(url1,url2):
    options=Options()
        #from selenium.webdriver.firefox.firefox_profile import FirefoxProfile
    #from selenium.webdriver.firefox.firefox_binary import FirefoxBinary
    #os.environ['MOZ_HEADLESS'] = '1'
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
    lb=1
    ub=5
    random1=lb+(ub-lb)*rand()
    driver.get(url1)
    time.sleep(random1)
    driver.get(url2)
    time.sleep(random1)
    driver.quit()
url1="https://www.mdpi.com/1999-4923/14/6/1181"
url2="https://www.mdpi.com/1999-4923/14/6/1181/pdf?version=1654149202"
url3="https://www.mdpi.com/2077-0375/12/4/434"
url4="https://www.mdpi.com/2077-0375/12/4/434/pdf?version=1652081937"
url5="https://www.mdpi.com/1999-4923/14/9/1897"
url6="https://www.mdpi.com/1999-4923/14/9/1897/pdf"
for i in range(3):
    CallPapers(url1,url2)
    CallPapers(url3,url4)
    CallPapers(url5,url6)
    time.sleep(10)

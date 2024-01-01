import datetime as dat
import matplotlib.pyplot as plt
import pandas_datareader as pdr
import numpy as np

plt.style.use("dark_background")
ma_1=30
ma_2=100

start=dat.datetime.now() -dat.timedelta(days=365*3)
end= dat.datetime.now()

data=pdr.DataReader("GLD","yahoo",start,end)
day= np.arange(1,len(data)+1)
data["day"]=day
#data.drop(columns=["Adj Close","Volume"],inplace=True)

#data=pdr.get_data_yahoo('GLD')
data[f'SMA_{ma_1}'] = data['Adj Close'].rolling(window=ma_1).mean().shift()
data[f'SMA_{ma_2}'] = data['Adj Close'].rolling(window=ma_2).mean().shift()



data=data.iloc[ma_2:]



buy_signals=[]
sell_signals=[]
trigger=0

for x in range(len(data)):
    if data[f'SMA_{ma_1}'].iloc[x]>data[f'SMA_{ma_2}'].iloc[x] and trigger!=1:
        buy_signals.append(data["Adj Close"].iloc[x])
        sell_signals.append(float("nan"))
        trigger=1
    elif data[f'SMA_{ma_1}'].iloc[x]<data[f'SMA_{ma_2}'].iloc[x] and trigger!=-1:
        sell_signals.append(data["Adj Close"].iloc[x])
        buy_signals.append(float("nan"))
        trigger=-1
    else:
        sell_signals.append(float("nan"))
        buy_signals.append(float("nan"))
data["Buy Signals"] = buy_signals
data["Sell Signals"] = sell_signals


data["return"]=np.log(data["Adj Close"]).diff()
data["signal"]=np.where(data[f'SMA_{ma_1}']>data[f'SMA_{ma_2}'],1,0)
data["signal"]=np.where(data[f'SMA_{ma_1}']<data[f'SMA_{ma_2}'],-1,data["signal"])
data["system_return"]=data["signal"]*data["return"]
data["entry"]=data.signal.diff()

fig1,ax1=plt.subplots()
fig2,ax2=plt.subplots()
ax1.plot(data['Adj Close'],label="Share Price", alpha=0.5)
ax1.plot(data[f'SMA_{ma_1}'],label=f'SMA_{ma_1}', color="orange", linestyle="--")
ax1.plot(data[f'SMA_{ma_2}'],label=f'SMA_{ma_2}', color="purple", linestyle="--")
ax1.scatter(data.index,data["Buy Signals"],label="Buy Signal",marker="^",color="#00ff00",lw=3)
ax1.scatter(data.index,data["Sell Signals"],label="Sell Signal",marker="v",color="#ff0000",lw=3)
plt.legend()
ax2.plot(np.exp(data["return"]).cumprod(),label='Buy/Hold', color="green")
ax2.plot(np.exp(data["system_return"]).cumprod(),label='System', color="blue")
plt.legend()

Price=data['Adj Close']


#mu=0.3

n=700

T=1

M=10

S0=Price[0]
#S0=100

#sigma=0.3

dt=T/n

logPrice=np.log(Price)
dlogPrice=np.diff(logPrice)
mu=np.mean(dlogPrice)/dt
sigma=np.sqrt(np.std(dlogPrice))

#St=np.exp((mu-sigma**2/2)*dt+sigma*np.random.normal(0,np.sqrt(dt),size=(M,n)).T)
St=np.exp(mu*dt+sigma*np.random.normal(0,np.sqrt(dt),size=(M,n)).T)



St=np.vstack([np.ones(M),St])

St=S0* St.cumprod(axis=0)

time=np.linspace(0,T,n+1)
tt= np.full(shape=(M,n+1),fill_value=time).T

fig3,ax3=plt.subplots()

ax3.plot(tt,St)

plt.show()





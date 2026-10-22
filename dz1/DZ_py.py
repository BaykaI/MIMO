import matplotlib.pyplot as plt
import numpy as np
def give_me_equal(num_r,num_t,crx,ctx,r):
    #print("xbckj 'ktvtynjd ", num_r*num_t)
    Ec=np.matrix([np.zeros(num_r*num_t),np.zeros(num_r*num_t)])
    for i in range(0,num_t):
        for j in range(0,num_r):
            m=num_r*(i)+j
            Ec[0,m-1]=crx[0,j]+ctx[0,i]+r[0]
            Ec[1,m-1]=crx[1,j]+ctx[1,i]+r[1]
    return Ec;
def giive_me_radiation_pattern(Ec,l):
    k=2*np.pi/l
    theta=np.linspace(-80,80,1440)
    DN_phi0=np.zeros_like(theta,dtype=complex)
    DN_phi90=np.zeros_like(theta)
    dB3=-3*np.ones_like(theta)
    M=Ec.shape[1] 
    r=np.zeros(M)
    arg=np.zeros(M)
    for i in range(0,M):
        r[i]=np.sqrt((Ec[0,i])**2+(Ec[1,i])**2)
        if (Ec[1,i]==0)and(Ec[0,i]==0):
            arg[i]=0
        elif((Ec[0,i]<0)and(Ec[1,i]==0)):
            arg[i]=np.pi
        else:   
            arg[i]=np.arctan(Ec[1,i]/Ec[0,i])
            

 
    for i in range (0,len(theta)):
        DN=0+0j;
        for jc in range(0,M):
            DN=DN+np.exp(1j*k*r[jc]*np.sin(theta[i]*np.pi/180)*np.cos(-arg[jc]))
        DN_phi0[i]=DN;
    
    plt.figure(figsize=(10, 4))
    plt.plot(theta, 10*np.log10((np.abs(DN_phi0)/np.abs(DN_phi0[int(np.size(theta)/2)]))**2), label="Диаграмма направленности")
    plt.plot(theta, dB3, color="red", label="Уровень -3дБ")
    plt.legend()
    plt.title("Диаграмма направленности")
    plt.xlabel("Угол, град.")
    plt.ylabel("Усиление, дБ")
    plt.grid(True)
    plt.show()
lyambda=1
flag=int(input('Введите 0 для ручного ввода и 1 для считывания из файла \n'))
if flag==0:
#задание координат Пр и Прд
    Nrx=int(input());
    Ntx=int(input());
    Rx=np.matrix([np.zeros(Nrx),np.zeros(Nrx)])
    Tx=np.matrix([np.zeros(Ntx),np.zeros(Ntx)])
    r0=np.zeros(2)
    r0[0]=float(input("x координата вектора смещения \n"))
    r0[1]=float(input("y координата вектора смещения \n"))

    for i in range (0,Nrx):
        print("x координата ", i+1, "-ого приёмника в длинах волн")
        Rx[0,i]=float(input())
        print("y координата ", i+1, "-ого приёмника в длинах волн")
        Rx[1,i]=float(input())
    for i in range (0,Ntx):
        print("x координата ", i+1, "-ого передатчика в длинах волн")
        Tx[0,i]=float(input())
        print("y координата ", i+1, "-ого передатчика в длинах волн")
        Tx[1,i]=float(input())
else:
    data = np.loadtxt(str(input('Укажите путь и расширения для файла данных \n')))   
    Rx_x = data[:, 0]
    Rx_y = data[:, 1]
    Tx_x = data[:, 2]
    Tx_y = data[:, 3]
    Tx_x=Tx_x[~np.isnan(Tx_x)]
    Tx_y=Tx_y[~np.isnan(Tx_y)]
    Rx=np.vstack((Rx_x,Rx_y))
    Tx=np.vstack((Tx_x,Tx_y))
    r0=np.zeros(2)
    r0[0]=float(input("x координата вектора смещения \n"))
    r0[1]=float(input("y координата вектора смещения \n"))

Equal_arr=give_me_equal(Rx.shape[1],Tx.shape[1],Rx,Tx,r0)
#Строим график
Rx_x = np.ravel(Rx[0])   
Rx_y = np.ravel(Rx[1])
Tx_x = np.ravel(Tx[0])
Tx_y = np.ravel(Tx[1])
Eq_x  = np.ravel(Equal_arr[0])
Eq_y  = np.ravel(Equal_arr[1])
plt.figure()#figsize=(10, 4))
plt.scatter(Rx_x, Rx_y, marker='o', color='red', label='Приёмники')
plt.scatter(Tx_x, Tx_y, marker='s', color='blue', label='Передатчики')
plt.scatter(Eq_x, Eq_y, marker='^', color='green', label='Эквивалент')
plt.title("Расстановка элементов")
plt.xlabel("Координаты элементов по оси ОХ в длинах волн")
plt.ylabel("Координаты элементов по оси ОУ в длинах волн")
#plt.grid(True)
plt.legend()
plt.show()
giive_me_radiation_pattern(Equal_arr,lyambda)

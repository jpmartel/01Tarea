"""
PREGUNTA 1
Programa que obtiene el espectro del Sol de 'sun_AM0.dat' 
y grafica flujo (en cgs) vs longitud de onda (micrómetros)
"""
import numpy as np
import matplotlib.pyplot as plt
datos=np.loadtxt('sun_AM0.dat') #abre archivo y obtiene arreglo [L_onda, Flujo]
londa=datos[:, 0]*0.001 #se obtienen longitudes de onda y convierte a nanómetros
flujo=datos[:, 1]*10000000*1000/10000 #se obtiennen flujos y se convierte a cgs
plt.xlim(0,3) #grafico hasta 3 micrómetros en longitud de onda
plt.plot(londa,flujo) #grafica
xlabel('Longitud de onda ($\mu m$)')
ylabel('flujo ($ergs$ $ s^{-1}$$cm^2$ $\mu m^{-1}$)')
title('Espectro solar')
grid(True)
show()



"""
PREGUNTA 2
Programa que calcula la integral del flujo del espectro solar
mediante la regla del trapezoide compuesta (algoritmo propio)
"""
import numpy as np
#obtiene arreglo doble con datos  de flujo y longitud de onda
#y Se obtienen las variables en arreglos independientes 
datos=np.loadtxt('sun_AM0.dat')
londa=datos[:, 0]*0.001 #en micrómetros
flujo=datos[:, 1]*10000000*1000/10000 #en cgs 
def trapecio_p2(londa,flujo):
    """
    función que recibe arreglo de longitud de onda y arreglo de flujo
    y entrega valor integral integrada numéricamente por método 
    trapezoides compuesto
    """
    t=0
    for i in range(len(londa)-1):
        h=(londa[i+1]-londa[i])
        t+=h*(flujo[i]+flujo[i+1])/2.
    return t
    
T=trapecio_p2(londa,flujo) #muestra el valor de la integral en ergs-1 s-1 cm-2
#cálculo luminosidad en ergs s-1
Luminosidad=T*4*np.pi*(1.49598*10**(13))**2
print Luminosidad



"""
PREGUNTA 3
Programa integra numéricamente función de Plank con una tolerancia escogida

"""
import numpy as np
from astropy import constants as cte
#valor analítico de la función dentro de la integral
valor=(np.pi)**4/15.
def f(y):
    """
    función a integrar (P sin las constantes)
    """
    return (np.tan(y)**5+np.tan(y)**3)/(np.exp(np.tan(y))-1)

def trapecio_p3(f,londa,h):
    """
    Recibe función a integrar, arreglo de valores para longitud de onda,
    y paso de integración
    Entrega valor para integral de f(x) calculado numéricamente con método de 
    trapezoides compuesto
    """
    fa=6 #primer valor obtenido con L'hopital
    s=h*(fa+f(londa[1]))/2. #trapezoide simple para primer intervalo
    for i in range(1,len(londa)-1): # calcula integral contrapezoide simple para cada intervalo
        h=(londa[i+1]-londa[i]) #paso integración       
        s=s+h*(f(londa[i])+f(londa[i+1]))/2.
    return s
tol=0.001 #tolerancia escogida
def trapecio_iterado(f,tol):  
    """
    Recibe función a integrar y tolerancia
    Entrega valor para integral de f(x) calculado numéricamente con método de 
    trapezoides compuesto que cumpla con la tolerancia escogida
    """
    a=0; b=np.pi/2 #valores extremos
    #trapezoide simple antes de aumentar número de intervalos
    h=(b-a)
    londa=np.linspace(a,b,2)
    I=trapecio_p3(f,londa,h)
    n=2
    while abs(valor-I)>=tol: #se aumenta número de intervalos hasta alcanzar precision deseada
        h=(b-a)/(2*(n-1))
        londa=np.linspace(a,b,n)
        I=trapecio_p3(f,londa,h)
        n=2**n
    return I
I=trapecio_iterado(f,tol)
print I
print valor
#ahora se optiene valor P al multiplicar I por las constantes
T=5776  #temperatura efectiva del Sol en kelvin
ctes=2*np.pi*cte.h.cgs/(cte.c.cgs)**2*(cte.k_B.cgs*T/cte.h.cgs)**4
P=ctes*I
print P
L=3.84185530786e+33 #luminosidad solar calculada en pregunta 2 en ergs s-1
R=np.sqrt(L/(4*np.pi*P))#Radio efectivo
print R



"""
PREGUNTA 4
scrip que calcula y muestra:
1-integral espectro solar de pregunta 1 con sp.integrate.trapz de scimpy
2-integral de de pregunta 2 con valor analítico pi^4/15 con sp.integrate.quad de scimpy
3-tiempo de ejecución integral algoritmo propio pregunta 2
4-tiempo de ejecución integral de módulo scimpy pregunta 2
5-tiempo de ejecución integral algoritmo propio pregunta 3
6-tiempo de ejecución integral de módulo scimpy pregunta 3
"""

mgc = get_ipython().magic
print sp.integrate.trapz(flujo,londa)
print sp.integrate.quad(f,0,np.pi/2)[0]*ctes
print mgc('%timeit trapecio_p2(londa,flujo)')
print mgc('%timeit sp.integrate.trapz(flujo,londa)')
print mgc('%timeit trapecio_iterado(f,tol)')
print mgc('%timeit sp.integrate.quad(f,0,np.pi/2)')
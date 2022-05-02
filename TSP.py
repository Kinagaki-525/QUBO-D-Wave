from dimod import sampleset
from pyqubo import solve_qubo, Array, Constraint, Placeholder, solve_qubo, LogEncInteger
import dimod
import math
import numpy as np
from dwave.system.composites import EmbeddingComposite
from dwave.system.samplers import DWaveSampler
from scipy.spatial import distance

N=5 

x = Array.create('x', (N, N), 'BINARY')#'BINARY' or 'SPIN'

Q=np.array([[1000,312,211,269,224],
            [312,1000,499,491,229],
            [211,499,1000,302,380],
            [269,491,302,1000,189],
            [224,229,380,189,1000]])

#Aichi-Tokyo 312km, Aichi-Osaka 211km, Aichi-Kanazawa 269km, Aichi-Matsumoto 224m
#Tokyo-Osaka 499km, Tokyo-Kanazawa 491km, Tokyo-Matsumoto 229km
#Osaka-Kanazawa 302km, Osaka-Matsumoto 380km
#Kanazawa-Matsumoto 189km

cost=0 
for i in range(N):
    for j in range(N):
        for k in range(N-1):
            cost=cost+Q[i][j]*x[i][k]*x[j][k+1]
for i in range(N):
    for j in range(N):
        cost=cost+Q[i][j]*x[i][N-1]*x[j][0]
cost=cost


#Constraints 1
constr_1=0
a=np.sum(x,axis=1)
for i in range(N):
    constr_1+=(a[i]-1)**2 
constr_1=constr_1

#Constraints 2
constr_2=0
b=np.sum(x,axis=0)
for i in range(N):
    constr_2+=(b[i]-1)**2
constr_2=constr_2

cost_func=cost+Placeholder('lam')*Constraint(constr_1, label='constr_1')+Placeholder('lam')*Constraint(constr_2, label='constr_2')
model=cost_func.compile()

lam=500

feed_dict={'lam':lam} #penalty
qubo, offset = model.to_qubo(feed_dict=feed_dict)

D_OBJECT = EmbeddingComposite(DWaveSampler(endpoint='https://cloud.dwavesys.com/sapi',
                                           token='', #your APItoken
                                           solver='DW_2000Q_6')) #Using Solver
count=1

resu=np.zeros((count))

responseFA = D_OBJECT.sample_qubo(qubo, num_reads=1, annealing_time=1)

    #anneal_schedule=[(0,0), (10,0.5), (110,0.5), (120,1)] or,annealing_time=1(μs)

result=responseFA.record[0][0].reshape(N,N)

c=10000
l=0
for (sample, energy, occurrence, _) in responseFA.data():
    if c>energy:
        l=occurrence
        c=energy
    elif c==energy:
        l=l+occurrence
    #print('Sample:{} Energy:{} Occurrence:{}'.format(list(sample.values()), energy, occurrence))
#print(l)#the number of Optimal solution

dista=0 #distance
for i in range(N):
    for j in range(N):
        for k in range(N-1):
            dista=dista+Q[i][j]*result[i][k]*result[j][k+1]
for i in range(N):
    for j in range(N):
        dista=dista+Q[i][j]*result[i][N-1]*result[j][0]

    #Constraints 1
a=np.sum(result,axis=1)
for i in range(N):
    dista+=((a[i]-1)**2 )*lam
        
    
    #Constraints 2
b=np.sum(result,axis=0)
for i in range(N):
    dista+=((b[i]-1)**2)*lam
        
print("distance:"+str(dista)+"km")
print(result)
print("Vertical:ctiy, Horizontal:order")
print("Visit the city when 1")
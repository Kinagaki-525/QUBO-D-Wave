#ナーススケジューリング問題（5人5日）をD-Waveで解く

from dimod import sampleset
from pyqubo import solve_qubo, Array, Constraint, Placeholder, solve_qubo, LogEncInteger
import dimod
import math
import numpy as np
from dwave.system.composites import EmbeddingComposite
from dwave.system.samplers import DWaveSampler
from scipy.spatial import distance
import random

N=5 #看護師の数
D=5 #総日数
J=3.5 #ペナルティ

x = Array.create('x', (N, D), 'BINARY')#'BINARY' or 'SPIN'


H=0
for n in range(N):
    for i in range(D-1):
        H=H+J*x[n][i]*x[n][i+1]

#h_1(n):nは人の識別番号、出勤回数の調整のためのパラメータ
def h_1(i):
    if i==0:
        h_1=1
    elif i==1:
        h_1=1
    elif i==2:
        h_1=2
    elif i==3:
        h_1=2
    elif i==4:
        h_1=3
    return h_1

#h_2(d):dは日付、その日の忙しさを考慮したパラメータ
def h_2(i):
    if i==0:
        h_2=1
    elif i==1:
        h_2=2
    elif i==2:
        h_2=1
    elif i==3:
        h_2=2
    elif i==4:
        h_2=1
    return h_2

#E(n):看護師nの労働力
def E(i):
    if i==0:
        E=1
    elif i==1:
        E=2
    elif i==2:
        E=2
    elif i==3:
        E=3
    elif i==4:
        E=3
    return E

#W(d):d日に必要なパラメータ
def W(i):
    if i==0:
        W=1
    elif i==1:
        W=2
    elif i==2:
        W=2
    elif i==3:
        W=3
    elif i==4:
        W=3
    return W

#F(n):看護師nの希望出勤数
def F(i):
    if i==0:
        F=1
    elif i==1:
        F=2
    elif i==2:
        F=2
    elif i==3:
        F=3
    elif i==4:
        F=3
    return F

lam=100

constr_1=0
for d in range(D):
    P=0
    for n in range(N):
        P=P+E(n)*x[n][d]
    constr_1=constr_1+(P-W(d))**2

constr_2=0
for n in range(N):
    Q=0
    for d in range(D):
        Q=Q+h_1(n)*h_2(d)*x[n][d]
    constr_2=constr_2+(Q-F(n))**2

H_func=H+Placeholder('lam')*Constraint(constr_1, label='constr_1')+Placeholder('lam')*Constraint(constr_2, label='constr_2')
model=H_func.compile()

feed_dict={'lam':lam} #ペナルティ関数
qubo, offset = model.to_qubo(feed_dict=feed_dict)

D_OBJECT = EmbeddingComposite(DWaveSampler(endpoint='https://cloud.dwavesys.com/sapi',
                                           token='', #あなたのAPIトークンを記入
                                           solver='DW_2000Q_6')) #利用可能なソルバーを指定
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
    elif c==energy: #=だと代入、==だと真偽
        l=l+occurrence
    #print('Sample:{} Energy:{} Occurrence:{}'.format(list(sample.values()), energy, occurrence))
#print(l)#最適解の数

print(result)
print("Vertical:nurse identification number, Horizontal:days")
print("Comes to work when 1")

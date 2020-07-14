""" 29 June 2020 """


import numpy as np
import math
from qiskit import *
from qiskit import Aer
from qiskit.circuit import Gate
from qiskit.visualization import plot_histogram
import statistics
import matplotlib

#set simulator
simulator = Aer.get_backend('qasm_simulator')

#coupling (dimensionless: G = ag)
G = 1.0

#parameter (W = aw)
W = 0.5
theta = 0

j = G*G/(4*W)

#fermion mass
m = 0.1
m_0 = m + 0.5

#time steps

#test vals
"""T = 1
M = 1*T"""
#real run
#T = 150
#M = 10*T
del_T = 0.1

#number of measurements/shots
SHOTS = 10**6
#repeats
REPEATS = 1

#number of qubits
N=4

print(" number of qubits: ", N,
      "\n coupling g = : ", G,
      "\n w = ", W,
      "\n j = ", j,
      "\n theta = ", theta,
      "\n init fermion mass m_0 = ", m_0,
      "\n actual fermion mass m = ", m,
      "\n delta_T = ", del_T,
      "\n no. of repeats = ", REPEATS,
      "\n no. of measurements = ", SHOTS,
      "\n"
     )
     
     
#____________________________________________________________
#time dependent constant functions

def h_X_Y(n,t,M):
    n = int(n)
    t = float(t)
    frac = t/M
    W_t = W*frac
    m_t = ((1 - frac)*m_0) + frac*m
    theta_t = theta*frac
    h = W_t - ((-1)**n)*m_t*math.sin(theta_t)/2
    return h

def h(n,t,M):
    n = int(n)
    t = float(t)
    frac = t/M
    m_t = (1 - frac)*m_0 + frac*m
    theta_t = theta*frac
    h = m_t*math.cos(theta_t)
    h = ((-1)**n)*h
    return h


#____________________________________________________________
#function to find VEV of mass operator

def VEV_m(M,R):

    VEV = []

    for r in range(R):
        """setup circuit"""
        qr = QuantumRegister(N)
        cr = ClassicalRegister(N)
        qc = QuantumCircuit(qr,cr)

        """set known vacuum"""
        for n in range (N):
            if (n+1)%2==0:
                qc.x(qr[n])

        """improved Suzuki-Trotter time evolution"""
        for t in range (M):
            for n in range (N-1):
                qc.rz(-math.pi/2, qr[n])
                qc.rz(-math.pi/2, qr[n+1])
                qc.cx(qr[n], qr[n+1])
                qc.rx(h_X_Y(n+1,t+1,M)*del_T/2, qr[n])
                qc.cx(qr[n], qr[n+1])
                qc.rz(math.pi/2, qr[n])
                qc.rz(math.pi/2, qr[n+1])                       #U_y/2


            for n in range (N-1):
                qc.cx(qr[n], qr[n+1])
                qc.rx(h_X_Y(n+1,t+1,M)*del_T/2, qr[n])
                qc.cx(qr[n], qr[n+1])                          #U_x/2



            for n in range (N):
                qc.rz(h(n+1,t+1,M)*del_T, qr[n])
            for n in range (1,N,2):
                for k in range (n):
                    qc.rz(-j*del_T, qr[k])                     #U_z


            for l in range (N-1):
                for k in range (l):
                    qc.cx(qr[k], qr[l])
                    qc.rz(j*del_T, qr[l])
                    qc.cx(qr[k], qr[l])                         #U_zz



            for n in range (N-1):
                qc.cx(qr[n], qr[n+1])
                qc.rx(h_X_Y(n+1,t+1,M)*del_T/2, qr[n])
                qc.cx(qr[n], qr[n+1])                           #U_x/2



            for n in range (N-1):
                qc.rz(-math.pi/2, qr[n])
                qc.rz(-math.pi/2, qr[n+1])
                qc.cx(qr[n], qr[n+1])
                qc.rx(h_X_Y(n+1,t+1,M)*del_T/2, qr[n])
                qc.cx(qr[n], qr[n+1])
                qc.rz(math.pi/2, qr[n])
                qc.rz(math.pi/2, qr[n+1])                          #U_y/2

        #qc.draw()
        #decomposed_circ = qc.decompose() # Does not modify original circuit
        #decomposed_circ.draw()

        """Output here is the time evolved vacuum"""

        """measure new vacuum"""

        qc.measure(qr,cr)
        result = execute(qc,backend = simulator,shots=SHOTS).result()
        #result.get_counts(qc)
        #plot_histogram(result.get_counts(qc))


        """label states using equivalent binary value, save probabilities as a list"""

        prob_dict = result.get_counts(qc)
        prob_states = np.zeros((2**N))
        for state in prob_dict.keys():
            num = int(state,2)
            prob_states[num] = prob_dict.get(state)/SHOTS
        #prob_states

        """measure VEV of mass operator"""

        vev = 0.0
        for state in prob_dict.keys():
            s = [int(x) for x in state]
            state_cbits = list(reversed(s))
            #print(state_cbits)
            op = 0
            for n in range(N):
                op = op + ((-1)**((n+1)+state_cbits[n]))
            #print(op)
            vev = vev + op*prob_states[int(state,2)]
            
            
        VEV.append(vev*W/N)
        del qc,qr,cr
        
    #res = [statistics.mean(VEV), statistics.stdev(VEV)]
    #return res
    return VEV


#____________________________________________________________
#"main" loop
T = list(range(150,250,10))
VEV_list = []
VEV_err = []

for time in T:
    #[vev_val, vev_err] = VEV_m(int(time/del_T),REPEATS)
    [vev_val] = VEV_m(int(time/del_T),REPEATS)
    VEV_list.append(vev_val)
    #VEV_err.append(vev_err)
    #print(time, ": completed, VEV = ", vev_val, ", error with ", REPEATS, " sets = ", vev_err)
    print(time, ": completed. VEV = ", vev_val)
    
print("T = ", T,"\n VEV = ",VEV_list,"\n avg VEV = ",statistics.mean(VEV_list))
    
#plotting
#import matplotlib.pyplot as plt
#plt.errorbar(T, VEV_list, yerr=VEV_err, fmt="o")
#plt.plot(T,VEV_list)
#plt.xlabel("t")
#plt.ylabel("VEV")
#plt.show()




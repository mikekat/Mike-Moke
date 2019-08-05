import random
import numpy as np
import numpy.linalg as la
import math						

def py_ang(v1, v2): 
   # Returns the angle in radians between np-vectors 'v1' and 'v2'   
    cosang = np.dot(v1, v2)
    sinang = la.norm(np.cross(v1, v2))
    return np.arctan2(sinang, cosang)

def is_in(v, M) : #tests if np-vector v is in the list of np-vectors M
    return tuple(v) in [tuple(x) for x in M]  

def removed(a,A): #removes np-vector from the list of np-vectors A
    y=[]
    for i in range(len(A)): 
        if  np.any(A[i]!=a) :
            y.append(A[i])
    return(y)


class Autopoiesis:
    def __init__(self, n, K): #takes as input the length of the matrix and the disintegration constant
        self.n=n
        self.K=K
        self.x=[]
        for i in range(self.n):
            self.x.append([])
        M=np.array([(0,-1),(-1,-1),(-1,0),(-1,1),(0,1),(1,1),(1,0),(1,-1)])
        for i in range(self.n):
            for j in range(self.n):
                y=list(input('Give: ' + str(i) + ',' + str(j) +' (K,S,H,L,BL) :')) #asks for the states in each position 
                if len(y)==2:
                    y[1]=M[int(y[1])]
                if len(y)==3:
                    y[1]=M[int(y[1])]
                    y[2]=M[int(y[2])]
                self.x[i].append(y)
            
    def R_p(self): #returns random direction from (West,North,East, South)
        return np.array(random.choice(((0,-1),(-1,0),(0,1),(1,0)))) 
                        
    def R_s(self, state, pos) : #creates a list of directions (West,North,East, South) from pos that are on the given state
                                #and returns one at random
                             
        M_1=(pos[0]-1,pos[1])
        M_2=(pos[0]+1,pos[1])
        M_3=(pos[0],pos[1]-1)
        M_4=(pos[0],pos[1]+1)
        D=[]
        for p in {M_1,M_2,M_3,M_4}  :
            if p[0] in range(self.n) and p[1] in range(self.n):
                if state != 'B' :
                    if self.x[p[0]][p[1]][0]==state :
                        D.append((p[0]-pos[0],p[1]-pos[1]))
                elif self.B(p) :
                    D.append((p[0]-pos[0],p[1]-pos[1]))
                    
        if D!=[]:
            return np.array(random.choice(D))
        else:
            return (-self.n-1,-self.n -1)

    def reorder(self, u, v):  #u=[position_1,...,position_m] , v=[f(1),...,f(n)] , f is a permutation of {1,...,n} 
        temp=[u[i] for i in v]
        A=np.array(self.x)
        for i in range(len(u)):
            self.x[u[i][0]][u[i][1]] = [z for z in A[temp[i][0]][temp[i][1]] ]
    
    #pos is a 2-element-np.array
    #each method tests if position pos is at some state (H,S,C,L,B,B1,B2)
    def H(self,pos):
        if pos[0] in range(self.n) and pos[1] in range(self.n)  :
            if self.x[pos[0]][pos[1]][0]=='H':
                return 1
            else:
                return 0
    def S(self,pos):
         if pos[0] in range(self.n) and pos[1] in range(self.n)  :
            if self.x[pos[0]][pos[1]][0]=='S'  :
                return 1
            else:
                return 0
    def L(self,pos):
        if pos[0] in range(self.n) and pos[1] in range(self.n)  :
            if self.x[pos[0]][pos[1]][0]=='L' and len(self.x[pos[0]][pos[1]])==1  :
                return 1
            else:
                return 0
    def C(self,pos):
        if pos[0] in range(self.n) and pos[1] in range(self.n)  :
            if self.x[pos[0]][pos[1]][0]=='K'   :
                return 1
            else:
                return 0
    def B1(self,pos):
        if pos[0] in range(self.n) and pos[1] in range(self.n)  :
            if len(self.x[pos[0]][pos[1]])==2 :
                return 1
            else:
                return 0
    
    def B2(self,pos):
        if pos[0] in range(self.n) and pos[1] in range(self.n)  :
            if len(self.x[pos[0]][pos[1]])==3 :
                return 1
            else:
                return 0
    def B(self,pos):
        if self.B1(pos) or self.B2(pos):
            return 1
        else:
            return 0
    
    def motion1(self):
        H_list=[]
        for i in range(self.n):
            for j in range(self.n):
                pos=np.array((i,j))
                if self.H(pos) :
                    H_list.append(pos)
        random.shuffle(H_list)
        y=[]
        for pos in H_list :
            R_r=self.R_p()
            if self.S(pos+R_r) or self.L(pos+R_r) or self.C(pos+R_r):
                u=[pos,pos+R_r]
                v=[1,0]
                self.reorder(u,v)
                if self.L(pos):
                    y.append(pos)
            if self.B(pos+R_r) and self.S(pos+2*R_r) :
                u=[pos,pos+R_r,pos+2*R_r]
                v=[2,1,0]
                self.reorder(u,v)
        for pos in y :
            self.bonding(pos)
    
    def motion2(self):
        L_list=[]
        for i in range(self.n):
            for j in range(self.n):
                pos=np.array((i,j))
                if self.L(pos) :
                    L_list.append(pos)
        random.shuffle(L_list)
        y=[]
        for pos in L_list :
            R_r=self.R_p()
            R_H=self.R_s('H',pos+R_r)
            R_B=self.R_s('B',pos+R_r)
            if self.S(pos+R_r) and self.H(pos+R_r+R_H) :
                u=[pos, pos+R_r, pos+R_r+R_H]
                v=[2,0,1]
                self.reorder(u,v)
                y.append(pos+R_r)
            elif self.S(pos+R_r) and self.B(pos+R_r+R_B) and self.H(pos+R_r+2*R_B) :
                u=[pos, pos+R_r, pos+R_r+R_B, pos+R_r+2*R_B]
                v=[3,0,2,1]
                self.reorder(u,v)
                y.append(pos+R_r)
            elif self.S(pos+R_r) or self.H(pos+R_r):
                u=[pos, pos+R_r]
                v=[1,0]
                self.reorder(u,v)
                y.append(pos+R_r)
        for p in y :
            self.bonding(pos)
    
    def motion3(self):
        C_list=[]
        for i in range(self.n):
            for j in range(self.n):
                pos = np.array((i,j))
                if self.C(pos) :
                    C_list.append(pos)
        random.shuffle(C_list)
        y=[]
        for pos in C_list :
            R_r=self.R_p()
            R_H=self.R_s('H',pos+2*R_r)
            R_B=self.R_s('B',pos+2*R_r)
            if self.L(pos+R_r) :
                if self.S(pos+2*R_r):
                    if self.H(pos+2*R_r+R_H):
                        u=[pos,pos+R_r,pos+2*R_r,pos+2*R_r+R_H]
                        v=[3,0,1,2]
                        self.reorder(u,v)
                        y.append(pos+2*R_r)
                    elif self.B(pos+2*R_r+R_B) and self.H(pos+2*R_r+2*R_B):
                        u=[pos,pos+R_r,pos+2*R_r,pos+2*R_r+R_B,pos+2*R_r+2*R_B]
                        v=[4,0,1,3,2]
                        self.reorder(u,v)
                        y.append(pos+2*R_r)
                elif self.H(pos+2*R_r):
                    u=[pos,pos+R_r,pos+2*R_r]
                    v=[2,0,1]
                    self.reorder(u,v)
                    y.append(pos+2*R_r)
            elif self.S(pos+R_r):
                R_H=self.R_s('H',pos+R_r)
                R_B=self.R_s('B',pos+R_r)
                if self.H(pos+R_r+R_H):
                    u=[pos,pos+R_r,pos+R_r+R_H]
                    v=[2,0,1]
                elif self.B(pos+R_r+R_B) and self.H(pos+R_r+2*R_B):
                    u=[pos,pos+R_r,pos+R_r+R_B,pos+R_r+2*R_B]
                    v=[3,0,2,1]
                    self.reorder(u,v)
                else:
                    u=[pos,pos+R_r]
                    v=[1,0]
                    self.reorder(u,v)
            if  self.L(pos+R_r) :
                if self.B(pos+2*R_r) or self.L(pos+2*R_r) or self.C(pos+2*R_r):
                    u=[pos,pos+R_r,pos+2*R_r]
                    v=[1,0,2]
                    y.append(pos+R_r)
            if self.H(pos+R_r):
                u=[pos,pos+R_r]
                v=[1,0]
                self.reorder(u,v)
                
        for p in y :
            self.bonding(p)

    def find_B(self, pos):
        list=[]
        M=np.array([(0,-1),(-1,-1),(-1,0),(-1,1),(0,1),(1,1),(1,0),(1,-1)])
        for d in M :
            if self.B(pos+d):
                list.append(d)
        return list

    def rename(self, u, v): #u=[position_1,...,position_m , v=[state_1,...,state_m] 
        i=0
        for pos in u :
            if self.B(pos) :
                for d in self.x[pos[0]][pos[1]][1:] :
                    p=pos+d
                    self.x[p[0]][p[1]] = removed(-d, self.x[p[0]][p[1]])
            self.x[pos[0]][pos[1]] = [v[i]]
            i+=1
            
    def production(self):
        M=np.array([(0,-1),(-1,-1),(-1,0),(-1,1),(0,1),(1,1),(1,0),(1,-1)])
        C_list=[]
        for i in range(self.n):
            for j in range(self.n):
                pos=np.array((i,j))
                if self.C(pos) :
                    C_list.append(pos)
        random.shuffle(C_list)
        y=[]
        for pos in C_list:
            S=[]
            for direction in M:
                if self.S(pos+direction):
                    S.append(direction)
            pairs=[]
            for p in S:
                for q in S:
                    if la.norm(p-q)==1:
                        pairs.append((p,q))
            if pairs !=[]:
                (R_S1,R_S2)=random.choice(pairs)
                u=[pos+R_S1,pos+R_S2]
                v=['L','H']
                self.rename(u,v)
                y.append(pos+R_S1)
                
        for p in y:
            self.bonding(p)
           
        
    def disintegration(self):
        M=np.array(((0,-1),(-1,0),(0,1),(1,0)))
        LB_list=[]
        for i in range(self.n):
            for j in range(self.n):
                pos = np.array((i,j))
                if self.L(pos) or self.B(pos) :
                    LB_list.append(pos)
        random.shuffle(LB_list)
        for pos in LB_list:
            R_N=random.random()
            if R_N <= self.K :
                R_H=self.R_s('H',pos)
                T1=[]
                for S1 in range(4):
                    if self.S(pos+M[S1]):
                        for i in {(S1-1)%4,S1%4,(S1+1)%4}:
                            if self.H(pos+M[S1]+M[i]):
                                T1.append(M[S1])


                T2=[]
                for S2 in range(4):
                    if self.S(pos+M[S2]):
                        for i in {(S2-1)%4,S2%4,(S2+1)%4}:
                            if self.B(pos+M[S2]+M[i]) and self.H(pos+M[S2]+2*M[i]) :
                                T2.append(M[S2])
                if self.H(pos+R_H):
                    u=[pos,pos+R_H]
                    v=['S','S']
                    self.rename(u,v)
                    self.rebond(pos)
                elif T1!=[]:
                    R_T1=random.choice(T1)
                    R_H=self.R_s('H',pos+R_T1)
                    u=[pos,pos+R_T1+R_H]
                    v=['S','S']
                    self.rename(u,v)
                    self.rebond(pos)
                elif T2!=[] :
                    R_T2=random.choice(T2)
                    R_B=self.R_s('B',pos+R_T2) #!!!
                    if self.H(pos+R_T2+2*R_B):
                        u=[pos,pos+R_T2+2*R_B]
                        v=['S','S']
                        self.rename(u,v)
                        self.rebond(pos)

    def find_L(self, pos): #creates a list the positions is the Moore neighbourhood of pos that are in the L state
        list=[]
        M=np.array([(0,-1),(-1,-1),(-1,0),(-1,1),(0,1),(1,1),(1,0),(1,-1)])
        for d in M :
            if self.L(pos+d):
                list.append(pos+d)
        return list

    def find_B1(self, pos): #creates a list the positions is the Moore neighbourhood of pos that are in the B1 state
        list=[]
        M=np.array([(0,-1),(-1,-1),(-1,0),(-1,1),(0,1),(1,1),(1,0),(1,-1)])
        for d in M :
            if self.B1(pos+d):
                list.append(pos+d)
        return list

    def bondable(self, p, q): #tests if positions p,q can be bonded
        key=0
        if np.any(q!=p) and (self.L(p) or self.B1(p)) and is_in(q, self.find_B1(p)+self.find_L(p) ):
            key=1
            if self.B1(p):
                v=self.x[p[0]][p[1]][1]
                u=q-p
                if py_ang(v,u) < math.pi/2 :
                    key=0
            if self.B1(q):
                v=self.x[q[0]][q[1]][1]
                u=p-q
                if py_ang(v,u) < math.pi/2 :
                    key=0
        return key

    def bond(self,p,q) : 
        if self.bondable(p,q) :
            self.x[p[0]][p[1]].append(q-p)
            self.x[q[0]][q[1]].append(p-q)

    def bonding(self,pos) :
        N = self.find_L(pos)
        M = self.find_B1(pos)
        l=len(M)
        
        
        for i in range(l) :
            if self.bondable(pos,M[i]) :
                M.append(M[i])
        M=M[l:]
        
        if M !=[] :
            m1=random.choice(M)
            self.bond(pos,m1)
            
        M = self.find_B1(pos)
        l=len(M)
        for i in range(l) :
            if self.bondable(pos,M[i]) :
                M.append(M[i])
        M=M[l:]
        
        if M !=[] :
            m2=random.choice(M)
            self.bond(pos,m2)
            
        l=len(N)
        for i in range(l) :
            if self.bondable(pos,N[i]) :
                N.append(N[i])
        N=N[l:]
        
        if N !=[] :
            n1=random.choice(N)
            self.bond(pos,n1)
                   
        N = self.find_L(pos)
        l=len(N)
        for i in range(l) :
            if self.bondable(pos,N[i]) :
                N.append(N[i])
        N=N[l:]        
        
        if N !=[] :
            n2=random.choice(N)
            self.bond(pos,n2)

    def rebond(self,pos) :
        L=self.find_B1(pos)
        pairs=[]
        for i in range(len(L)) :
            for j in range(len(L)) :
                if self.bondable(L[i],L[j]):
                    pairs.append((L[i],L[j]))
        while pairs != [] :
            random_pair = random.choice(pairs)
            self.bond(random_pair[0],random_pair[1])
            L=self.find_B1(pos)
            pairs=[]
            for i in range(len(L))  :
                for j in range(len(L)) :
                    if self.bondable(L[i],L[j]):
                        pairs.append((L[i],L[j]))
        L = self.find_B1(pos)+self.find_L(pos)
        pairs=[]
        for i in range(len(L)) :
            for j in range(len(L)) :
                if self.bondable(L[i],L[j]):
                    pairs.append((L[i],L[j]))
        while pairs != [] :
            random_pair = random.choice(pairs)
            self.bond(random_pair[0],random_pair[1])
            L=self.find_B1(pos)+self.find_L(pos)
            pairs=[]
            for i in range(len(L)) :
                for j in range(len(L)) :
                    if self.bondable(L[i],L[j]):
                        pairs.append((L[i],L[j]))

    def __str__(self) : 
        M=np.array(((-1,1),(0,1),(1,1),(1,0)))
        D=('/','-','\\','|')
        m=[[' ' for i in range(3*self.n -2)] for j in range(3*self.n -2)]
        for i in range(self.n) :
            for j in range(self.n) :
                m[3*i][3*j] = self.x[i][j][0]
                for k in range(len(M)) :
                    if is_in(M[k], self.x[i][j] ) :
                        p=3*np.array((i,j)) + M[k]
                        q=3*np.array((i,j)) + 2*M[k]
                        m[p[0]][p[1]]=D[k]
                        m[q[0]][q[1]]=D[k]
        S='\n'.join([''.join(m[i]) for i in range(len(m))]) +'\n' +  len(m)*'-'
        return S
                                    
        
def main(t_s) :
    n=int(input('Give length : '))
    K=float(input('Give disintegration rate: '))
    A=Autopoiesis(n,K)
    print(0)
    print(A)
    for i in range(t_s):
        print(i+1)
        A.motion1()
        
        A.motion2()
       
        A.motion3()
       
        A.production()
   
        A.disintegration()
        
        print()
        print(A)

t_s=int(input('Give number of time steps: '))       
main(t_s)


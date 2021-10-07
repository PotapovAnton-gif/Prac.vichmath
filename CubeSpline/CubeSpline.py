import pandas as pd
import matplotlib.pyplot as plt
import numpy as np 


class CubicSpline():
    
    
    
    def __init__ (self, file):
        
        self.data = pd.read_csv(file)

        self.xes = self.data["Hight"].values
        self.yes = self.data["plotnost"].values


    def tridiag_matrix_alg(self, A, F):
        
        n = A.shape[0]
        
        al = np.zeros(shape=(n))
        bet = np.zeros(shape =(n))
        
        al[1] = -A[0,1]/A[0,0]
        bet[1] = F[0]/A[0,0]
        
        for i in range(1,n-1):
        
            a,b,c,f = A[i,i-1],A[i,i],A[i,i+1],F[i]
            al[i+1]=-c/(a*al[i]+b)
            bet[i+1]=(f-a*bet[i])/(a*al[i]+b)
        
        x_es = np.zeros(shape = n)
        x_es[n-1] = (F[n-1]-A[n-1,n-2]*bet[-1])/(A[n-1,n-1]+A[n-1,n-2]*al[-1])
        
        for i in range(n-2,-1, -1):
        
            x_es[i]=al[i+1]*x_es[i+1]+bet[i+1]
        
        return x_es
    
    
    def h(self, i):
        
        if i == 0:
        
            return self.xes[1]-self.xes[0]
        
        return self.xes[i]-self.xes[i-1]
    
    
    def runge_ex(self, xes):
        
        y = []
        
        for (i,x) in enumerate(self.xes):
        
            print(1/(1+x))
            y.append(1/(1+x^2))
        
        return y

    
    def interpolate(self):
       
        plt.plot(self.xes,self.yes)

        n=len(self.xes)-1

        sp_cof = np.zeros(shape=(n+1, 4))
        sp_cof[:,0] = self.yes

        F = []

        for i in range(1,n):
        
            F.append(3*( (sp_cof[i+1,0]-sp_cof[i,0])/self.h(i+1) - (sp_cof[i,0]-sp_cof[i-1,0])/(self.h(i))) )

        

        A=np.zeros(shape = (n-1,n-1))
        
        for i in range(n-1):
            
            if i!=0 and i!=n-2:
            
                A[i,i-1],A[i,i],A[i,i+1] = self.h(i+1), 2*(self.h(i+2)+self.h(i+1)),self.h(i+2)
            
            elif i==0:
            
                A[0,0],A[0,1]=4*self.h(i),self.h(i)
            
            elif i == n-2:
            
                A[n-2,n-3],A[n-2,n-2] = self.h(i),4*self.h(i)
        
        
        
        sp_cof[1:n,2] = self.tridiag_matrix_alg(A,F)
        print(sp_cof)
        
        
        for i in range(n+1):
        
            sp_cof[i,3] = (sp_cof[i,2]-sp_cof[i-1,2])/(3*self.h(i))
            sp_cof[i,1] = (sp_cof[i,0]-sp_cof[i-1,0])/self.h(i) +(2*sp_cof[i,2]+sp_cof[i-1,2]) *self.h(i)/3
        
        
        plt.ylim([0,0.2e-10])
        
        for i in range(1,n+1):
            
            [a,b,c,d] = sp_cof[i,:]
            x_s = np.linspace(self.xes[i-1],self.xes[i],100)
            y_s = a + b*(x_s-self.xes[i])+c*(x_s-self.xes[i])**2+d*(x_s-self.xes[i])**3
            plt.plot(x_s,y_s)
            

        plt.show()

if __name__ == "__main__":
    cub_spline = CubicSpline("Atmosphere_1.csv")
    cub_spline.interpolate()




 
        

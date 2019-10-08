class ran(object):   
    
#     The object class ran() generates pseudo random numbers using either the linear congrugential algorithm, or the Mersenne  
#     Twister algorithm (in the first case, using my own implementation, in the second one, by calling np.random.uniform()). 
#     To choose between one or the other, pass the string "LGC" or "Twistter" when creating an instance of the class, in example:

#         ran=ran('LGC')

#     If no string is passed, the Mersenne Twistter algorithm will be used

#     In the case of the LGC, the parameters for the multipliyer 'a', additive constant 'c' and the modulus 'm' can be choosen 
#     among 4 options: Numerical Recipes, Borland, glibc, and Delphi, the values being retrieved from the wikipedia LCG page. 

#     By default, Numerical Recipes will apply, but the user can generate pseudo random numbers using any of the others 
#     by calling the method set_params and passing one of the following strings: 'Borland', 'glibc' or 'Delphi'.

#     In example:
#         ran.set_params('Delphi')

#     The methods ran.ranu(), ran.rang() and ran.rane() return random numbers drawn from, respectively, the uniform, 
#     gaussian, or exponential distributions. 

#     If the are called without passing any argument, only one random number is returned. On the other hand, if they are called
#     passing an integer n, they return n random numbers drawn from the corresponding distribution.

#     The user can set the edges of the interval, the mean and standard deviation, and the tau value (respectively in the uniform, 
#     gaussian and exponential cases) by calling the methods set_interval(min,max), set_gauss(mu,sigma), and set_tau(tau)

#     In example:

#         ran.set_interval(3,5)
#         ran.ranu(100)

#     will generate 100 uniformly distributed numbers in the (3,5) interval
    

    def __init__(self,algorithm='Twistter'):
        self.iseed=None
        self.min,self.max=0,1
        self.avg,self.std=0,1
        self.tau=1
        self.algorithm=algorithm
        self.a,self.c,self.m,self.norm=1664525,1013904223,2**32,np.sqrt(5.5)/9999999999
        if(self.algorithm!='Twistter' and self.algorithm!='LCG'): 
            raise ValueError("Bad algorithm String passed. Please pass either 'Twistter' or 'LCG'")
    
#     Setter method for the seed. If is called without given any integer value for the seed, a seed is generated from the 
#     current seconds in the system clock. If the seed is a number with less than 7 digists, it is padded with the
#     relevant number of 0 on the right until its made out of 7 digits.

    def set_seed(self,iseed):
#         IN: iseed, integer
        self.iseed=iseed
        if(type(self.iseed)!=int):
            x=str(time.time())
            self.iseed=x[re.search('\.', x).start()+1:]
            while(len(self.iseed)<7):self.iseed=self.iseed + '0'
            self.iseed=int(self.iseed)
        if(self.algorithm=='Twistter'):
            np.random.seed(self.iseed)
        print('final seed state', type(self.iseed))
        
    def set_params(self,params):
#         IN: parameters for the Linear Congruguential Algorithm
        self.params=params
        if(self.params=='Numerical Recipes'): self.a,self.c,self.m,self.norm=1664525,1013904223,2**32,np.sqrt(5.5)/9999999999
        elif(self.params=='Borland'): self.a,self.c,self.m,self.norm=22695477,1,2**32,np.sqrt(5.5)/9999999999
        elif(self.params=='glibc'):self.a,self.c,self.m,self.norm=1103515245,12345,2**32,np.sqrt(5.5)/9999999999
        elif(self.params=='Delphi'):self.a,self.c,self.m,self.norm=134775813,1,2**32,np.sqrt(5.5)/9999999999
    
    def set_interval(self,a,b):
#         IN: values for the interval in the uniform distribution
        self.min,self.max=a,b
        
#     By calling this method the current seed in the sequence can be retreived. If no seed has been yet settled, a seed
#     is generated from the secons in the system clock and then returned.
    def get_seed(self):
#         OUT: If the LCG algorithm is being used, the previous produced random number in the sequence (i.e. current seed).
#             If the Twistter algorithm is being used, the seed that was used to initialise the generator is retrieved, as
#             because of the algebra on which the Twistter generator relies, there is no concept of "seed", but instead an
#             array that can be retrieved calling np.rand.get_state()
        if(self.iseed==None): 
            self.set_seed(self)
        return self.iseed

#     This method is internally called within the class to generate uniforme distributed random values using the LCG recurrence
#     relation. As it is an internal method, I decided to define it as private, because never will be call outside the calss.
#     The result is multiplied by a normalization factor in order to get the numbers in the interval [0,1]
    def __r(self):
        if(self.algorithm=='Twistter'):
            if(self.iseed==None): self.set_seed(self)
            return np.random.uniform(self.min,self.max)
        elif(self.algorithm=='LCG'):
            if(self.iseed==None): 
                self.set_seed(self)
            x=(self.a*self.iseed+self.c)%self.m
            self.iseed=x
            return (self.max-self.min)*x*self.norm+self.min

#     If the user calls the ranu method passing an integer n, an array of n float numbers uniformly distributed in 
#     (in the case of the LGC, roughly) the range [0,1].
    def ranu(self,n=1):
#         IN: Integer n, number of uniform distributed random numbers to be returned, by default 1
#         OUT: random number or array of random numbers drawn from uniform distribution between the parameters settled by 
#         the method set_interval(a,b)
        if(n>1): return [self.__r() for i in range(n)]
        else: return self.__r() 

#     This method is internally called within the class to generate random values drawn from the Gaussian normal distribution
#     (with mean = 0 and stdev = 1) using the Box Muller transform algorithm. Again, as it is an internal method, it is defined
#     as private.
    def __rg(self):
        W=1.1
        self.set_interval(-1,1)
        while(W>1 or W==0):
            u,v=self.__r(),self.__r()
            W=np.power(u,2)+np.power(v,2)
        A=np.sqrt(-2*np.log(W)/W)
        rg=u*A
        rg=rg*self.std+self.avg
        return rg
    
#     Bellow is a setter method for the exponential distribution parameter Tau.
    def set_tau(self,tau):
#         IN: Tau, by default one.
        self.tau=tau

#   Bellow is a setter method for the parameters of the Gaussian distribution:
    def set_gauss(self,avg=0,std=1):
# IN: avg, mean of the distribution; std, standard deviation of the distribution, by default 0 and 1 respectively, corresponding
# to the normal distribution.
        self.avg,self.std=avg,std
    
    def rang(self,n=1):
#         IN: Integer n, number of gaussian distributed random numbers to be returned, by default 1
#         OUT: random number or array of random numbers drawn from Gaussian distribution between the parameters settled by 
#         the method set_gauss(a,b)
        if(n>1): return [self.__rg() for i in range(n)]
        else: return self.__rg() 

#     The private method bellow is internally called within the class to generate random values drawn from the Exponential normal distribution
#     using the Inverse Transform Sampling approach.
    def __re(self):
        u=1.1
        self.set_interval(0,1)
        while((1-u)<0): u=self.__r()
        re=-(self.tau)*np.log(1-u)
        return re
    
    def rane(self,n=1):
#         IN: Integer n, number of exponential distributed random numbers to be returned, by default 1
#         OUT: random number or array of random numbers drawn from Exponential distribution between the parameters settled by 
#         the method set_tau(t)
        if(n>1): return [self.__re() for i in range(n)]
        else: return self.__re()
    

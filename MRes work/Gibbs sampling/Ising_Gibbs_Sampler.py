#This code performs Gibbs sampling on an Ising model (NxN graphical model)



import numpy as np 
#We attempt an easier problem first, which is that of a 2x2 system
N=2
sigma_0=np.ones(4) #initialise 
sigma_t=sigma_0
sigma_tplus1=sigma_0
iteration_no=100000 #want to do this many iterations, truly does x-1, excludes initialised case
beta=10 #set beta above the critical value, magnetisation should be almost 1
sample_holder=np.zeros((iteration_no,N**2)) #initialise matrix to hold the samples
sample_holder[0,:]=sigma_t #set first row equal to initial value
#Perform Gibbs sampling
for i in range(1,iteration_no):
    #Approach block with 1 and 4 first
    argument=4*beta*(sigma_t[1]+sigma_t[2])
    jumpprob1_4=(np.exp(argument)+1)/(2*np.cosh(argument)+2)
    one_go_to_1=np.random.uniform(0,1,1)<jumpprob1_4
    four_go_to_1=np.random.uniform(0,1,1)<jumpprob1_4
    if one_go_to_1:
        sigma_tplus1[0]=1
    else:
        sigma_tplus1[0]=-1
    if four_go_to_1:
        sigma_tplus1[3]=1
    else:
        sigma_tplus1[3]=-1
    #Approach block with 2 and 3
    argument=4*beta*(sigma_tplus1[0]+sigma_tplus1[3])
    jumpprob2_3=(np.exp(argument)+1)/(2*np.cosh(argument)+2)
    two_go_to_1=np.random.uniform(0,1,1)<jumpprob2_3
    three_go_to_1=np.random.uniform(0,1,1)<jumpprob2_3
    if two_go_to_1:
        sigma_tplus1[1]=1
    else:
        sigma_tplus1[1]=-1
    if three_go_to_1:
        sigma_tplus1[2]=1
    else:
        sigma_tplus1[2]=-1
    #Now store and update values
    sample_holder[i,:]=sigma_tplus1 #store generated values
    sigma_t=sigma_tplus1 #update values for next iteration
#Now find averages
print(sample_holder[iteration_no-10:iteration_no,:])
sample_holder=sample_holder.astype(int)
count_1=np.count_nonzero(sample_holder==1)
count_minus1=np.count_nonzero(sample_holder==-1)
print(count_1)
print(count_minus1)
expectation=np.mean(sample_holder,axis=0)
print(expectation)
m=0.25*np.sum(expectation) #This is the magnetisation
print('Magnetisation for 2x2 system subject to beta=',beta,' is:',m) #Should be near zero if beta<0.44


def GibbsSamplerIsing(N,beta,sigma_initial,iteration_no):
    sigma_current=sigma_initial #initialise vectors
    sigma_next=sigma_initial
    sample_holder=np.zeros((iteration_no,N**2)) #matrix to hold samples
    sample_holder[0,:]=sigma_current #set first row equal to initial value
    for t in range(1,iteration_no):
        if N%2==1: #in this case, N is an odd number
            for i in range(1,N**2+2,2): #over these odd numbers first
                #Now need to identify if point on edge or corner
                left=(i-1)%N==0 #if true, on left edge of grid
                right=(i+1)%N==1 #if true, on right edge of grid
                top=(i-N)<=0 #if true, on top edge of grid
                bottom=(i+N)>N**2 #if true, on bottom edge of grid

                if left:
                    if top:
                        #Point 1 here
                        argument=2*beta*(sigma_current[i]+sigma_current[i+N-1])
                    elif bottom:
                        #Point N(N-1)+1 here
                        argument=2*beta*(sigma_current[i]+sigma_current[i-N-1])
                    else:
                        #Point on left-hand side, not on corner
                        argument=2*beta*(sigma_current[i]+sigma_current[i+N-1]+sigma_current[i-N-1])
                elif right:
                    if top:
                        #Point N here
                        argument=2*beta*(sigma_current[i-2]+sigma_current[i+N-1])
                    elif bottom:
                        #Point N^2 here
                        argument=2*beta*(sigma_current[i-2]+sigma_current[i-N-1])
                    else:
                        #Point on right-hand side, not corner
                        argument=2*beta*(sigma_current[i-2]+sigma_current[i+N-1]+sigma_current[i-N-1])
                elif top:
                    #Top side, not corners
                    argument=2*beta*(sigma_current[i-2]+sigma_current[i]+sigma_current[i+N-1])
                elif bottom:
                    #Bottom side, not corners
                    argument=2*beta*(sigma_current[i-2]+sigma_current[i]+sigma_current[i-N-1])
                else:
                    #Point on the interior of the grid
                    argument=2*beta*(sigma_current[i-2]+sigma_current[i]+sigma_current[i+N-1]+sigma_current[i-N-1])
                prob_1=np.exp(argument)/(2*np.cosh(argument))
                jump_1=np.random.rand()<prob_1
                if jump_1:
                    sigma_next[i-1]=1
                else:
                    sigma_next[i-1]=-1
            #Repeat for even numbered points
            for i in range(2,N**2,2): 
                # Identify if point on edge (note there is no corners is N odd)
                left=(i-1)%N==0 #if true, on left edge of grid
                right=(i+1)%N==1 #if true, on right edge of grid
                top=(i-N)<=0 #if true, on top edge of grid
                bottom=(i+N)>N**2 #if true, on bottom edge of grid
                if left:
                   argument=2*beta*(sigma_next[i]+sigma_next[i+N-1]+sigma_next[i-N-1])
                elif right:
                   argument=2*beta*(sigma_next[i-2]+sigma_next[i+N-1]+sigma_next[i-N-1])
                elif top:
                   argument=2*beta*(sigma_next[i-2]+sigma_next[i]+sigma_next[i+N-1])
                elif bottom:
                   argument=2*beta*(sigma_next[i-2]+sigma_next[i]+sigma_next[i-N-1])
                else:
                   argument=2*beta*(sigma_next[i-2]+sigma_next[i]+sigma_next[i+N-1]+sigma_next[i-N-1])
                prob_1=np.exp(argument)/(2*np.cosh(argument))
                jump_1=np.random.rand()<prob_1
                if jump_1:
                    sigma_next[i-1]=1
                else:
                    sigma_next[i-1]=-1
        else: 
            #In this case, N is an even number.
            #Indices of blue/red points more complicated, make a vector of indices that are blue (or red)
            indices_holder_blue=np.zeros(int(N**2/2))
            row_indicator=1
            for i in range(1,N+1):
                if row_indicator%2==1: #row number is odd
                    indices_holder_blue[int((N/2)*(row_indicator-1)):int(N/2*(row_indicator))]=(row_indicator-1)*N*np.ones(int(N/2))+np.arange(1,N+1,2)
                else: #row number is even
                    indices_holder_blue[int((N/2)*(row_indicator-1)):int(N/2*(row_indicator))]=(row_indicator-1)*N*np.ones(int(N/2))+np.arange(2,N+2,2)
                row_indicator=row_indicator+1
            indices_holder_blue=indices_holder_blue.astype(int)

            indices_holder_red=np.zeros(int(N**2/2))
            row_indicator=1
            for i in range(1,N+1):
                if row_indicator%2==1: #row number is odd
                    indices_holder_red[int((N/2)*(row_indicator-1)):int(N/2*(row_indicator))]=(row_indicator-1)*N*np.ones(int(N/2))+np.arange(2,N+2,2)
                else: #row number is even
                    indices_holder_red[int((N/2)*(row_indicator-1)):int(N/2*(row_indicator))]=(row_indicator-1)*N*np.ones(int(N/2))+np.arange(1,N+1,2)
                row_indicator=row_indicator+1
            indices_holder_red=indices_holder_red.astype(int)
            #Continue with problem now
            for i in indices_holder_blue:
                #Now need to identify if point on edge or corner
                left=(i-1)%N==0 #if true, on left edge of grid
                right=(i+1)%N==1 #if true, on right edge of grid
                top=(i-N)<=0 #if true, on top edge of grid
                bottom=(i+N)>N**2 #if true, on bottom edge of grid

                if left:
                    if top:
                        #Point 1 here
                        argument=2*beta*(sigma_current[i]+sigma_current[i+N-1])
                    else:
                        #Point on left-hand side, not on corner
                        argument=2*beta*(sigma_current[i]+sigma_current[i+N-1]+sigma_current[i-N-1])
                elif right:
                    if bottom:
                        #Point N^2 here
                        argument=2*beta*(sigma_current[i-2]+sigma_current[i-N-1])
                    else:
                        #Point on right-hand side, not corner
                        argument=2*beta*(sigma_current[i-2]+sigma_current[i+N-1]+sigma_current[i-N-1])
                elif top:
                    #Top side, not corners
                    argument=2*beta*(sigma_current[i-2]+sigma_current[i]+sigma_current[i+N-1])
                elif bottom:
                    #Bottom side, not corners
                    argument=2*beta*(sigma_current[i-2]+sigma_current[i]+sigma_current[i-N-1])
                else:
                    #Point on the interior of the grid
                    argument=2*beta*(sigma_current[i-2]+sigma_current[i]+sigma_current[i+N-1]+sigma_current[i-N-1])
                prob_1=np.exp(argument)/(2*np.cosh(argument))
                jump_1=np.random.rand()<prob_1
                if jump_1:
                    sigma_next[i-1]=1
                else:
                    sigma_next[i-1]=-1
            for i in indices_holder_red:
                # Identify if point on edge or at corner
                left=(i-1)%N==0 #if true, on left edge of grid
                right=(i+1)%N==1 #if true, on right edge of grid
                top=(i-N)<=0 #if true, on top edge of grid
                bottom=(i+N)>N**2 #if true, on bottom edge of grid
                if left:
                    if bottom:
                        argument=2*beta*(sigma_next[i]+sigma_next[i-N-1])
                    else:
                        argument=2*beta*(sigma_next[i]+sigma_next[i+N-1]+sigma_next[i-N-1])
                elif right:
                    if top:
                        argument=2*beta*(sigma_next[i-2]+sigma_next[i+N-1])
                    else:
                        argument=2*beta*(sigma_next[i-2]+sigma_next[i+N-1]+sigma_next[i-N-1])
                elif top:
                   argument=2*beta*(sigma_next[i-2]+sigma_next[i]+sigma_next[i+N-1])
                elif bottom:
                   argument=2*beta*(sigma_next[i-2]+sigma_next[i]+sigma_next[i-N-1])
                else:
                   argument=2*beta*(sigma_next[i-2]+sigma_next[i]+sigma_next[i+N-1]+sigma_next[i-N-1])
                prob_1=np.exp(argument)/(2*np.cosh(argument))
                jump_1=np.random.rand()<prob_1
                if jump_1:
                    sigma_next[i-1]=1
                else:
                    sigma_next[i-1]=-1 
        #Now we have the full iteration at the next time point, need to store it 
        sample_holder[t,:]=sigma_next
        sigma_current=sigma_next #update the iteration
    
    expectation=np.mean(sample_holder,axis=0)
    m=(1/N**2)*np.sum(expectation) #This is the magnetisation
    return sample_holder,m

#Now we want to test this function by using it for problem 11b
N=3 #asks us to use 100x100 grid
beta=np.linspace(0.1,10.1,num=4) #vector of beta values to test
sigma_initial=np.ones(N**2)
iteration_numbers=np.array([100,100,100,100])
m_holder2=np.zeros(len(beta))
for i in range(0,len(beta)):
    sample,magnetisation=GibbsSamplerIsing(N,beta[i],sigma_initial,iteration_numbers[i]) #perform sampling for each beta
    print(sample)
    del(sample)
    m_holder2[i]=magnetisation       
    del(magnetisation)       

beta_m_holder2=np.zeros((len(beta),2))
beta_m_holder2[:,0]=beta
beta_m_holder2[:,1]=m_holder2
print(beta_m_holder2)
                

#Code is correct, but the magnetisation was defined incorrectly in the notes which led to an incorrect output here. 
#Analysis of the behaviour of the system using the above is as expected.


             

                        
                




















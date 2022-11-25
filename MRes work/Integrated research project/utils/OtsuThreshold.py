#This python file uses Otsu's algorithm to calculate a threshold
import numpy as np

def OtsuThreshold(counts,threshold):
    #Find total number of pixels in image
    total_pixels=np.sum(counts)
    #Calculate probabilities from histogram information
    probs=counts/total_pixels
    #Calculate class probabilities at given threshold
    w_0=np.sum(probs[0:threshold]) #remember max element not included in vector
    w_1=1-w_0
    #Calculate class means 
    mu_0=np.dot(np.arange(threshold),probs[0:threshold])/w_0
    mu_1=np.dot(np.arange(threshold,256),probs[threshold:256])/w_1
    #Return inter-class variance
    inter_cv=w_0*w_1*(mu_0-mu_1)**2
    return inter_cv

    


    
    

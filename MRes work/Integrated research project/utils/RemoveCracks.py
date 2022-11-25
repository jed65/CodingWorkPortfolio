#This python file uses the signal to remove cracks from an image
import numpy as np

def RemoveCracks(isolated_row,thresh_width):
    black_pixels=np.where(isolated_row==0) #gives the indices in isolated_row that represent black in image
    black_pixels2=black_pixels[0] #this is so we can index on the values
    #Now need to find widths of the above
    width_vector=np.array([0]) #first value is zero (can be ignored)
    i=0 #initialise variables
    counter=0
    while counter<len(black_pixels2): #want it to check all of the values in the vector
         width_counter=1
         i=counter #ensures that we start at the first black pixel
         while black_pixels2[i+1]==black_pixels2[i]+1: #i.e. chain of black pixels
            width_counter=width_counter+1 #add one to the width
            i=i+1 #to check the next element
            if i>=len(black_pixels2)-1:
                break #leave the while loop before it spits an error
         width_vector=np.append(width_vector,width_counter) #append the width of this signal amplitude
         counter=counter+width_counter #add the next width on so we keep track of where we get up to
    indices=np.arange(1,len(width_vector))
    tracker=0
    for i in indices:
        if width_vector[i]<thresh_width:  #in this case, want to set these parts to white
            set_white=black_pixels2[tracker:tracker+width_vector[i]] #indices that need changing 
            isolated_row[set_white]=1 #change them to white
        tracker=tracker+width_vector[i]
    
    return isolated_row


#import sys
#sys.path.append(r'C:\Users\yc2634\AppData\Local\Programs\Python\Python39\Lib\site-packages')

import numpy as np
from PIL import Image as im
import matplotlib.pyplot as plt
import matplotlib
import time



# read the raw image from vol file
fname = r'C:\Users\jaked\OneDrive\Documents\MRes Course\Integrated Research Project\AM1-Lam16,001,007,008,0013,8.3-1.16_subvol_1530x1020x161_8bit.raw'
nrows = 1020
ncols = 1530 
nsls = 161
dtype = np.uint8
vxsiz = 0.0332592 #mm

#
fd = open(fname)
print('reading image data:',fname)
t = time.time()
A = np.fromfile(fd, dtype=dtype, sep="", count=nrows*ncols*nsls)
fd.close()
A = A.reshape([nsls,nrows,ncols])
print('---> reading completed, time consumed: ', time.time()-t,'s')

# Take a look at the images
from utils.IndexTracker import *
fig, ax = plt.subplots(1, 1)
tracker = IndexTracker(ax,np.moveaxis(A,0,-1))
fig.canvas.mpl_connect('scroll_event', tracker.on_scroll)
plt.show()

#Take a slice and look at its histogram and an image of it
A_slice=A[77,:,:] #take the 78th slice

#Histogram of this slice
counts,bins=np.histogram(A_slice,bins=256) #get counts and bins
plt.stairs(counts, bins) #produce plot
plt.show() #show plot

#Image of this slice
slice=im.fromarray(A_slice)
slice.save('slice_78.png')

#Perform a thresholding based upon histogram
A_slice_threshold=A_slice > 125 #Makes all points above (below) this value black (white)
slice_threshold=im.fromarray(A_slice_threshold)
slice_threshold.save('slice_78_threshold_125.png') #Take a look at the image (saves to folder)


#Now wish to perform Otsu thresholding algorithm, need counts and total pixels in image
from utils.OtsuThreshold import *
maximum=OtsuThreshold(counts,1) #initialise maximum and threshold value
otsu_thr=1
#Now loop over different thresholds, update if new inter-class variance exceeds previous maximum
for i in np.arange(2,257):
    inter_class_variance=OtsuThreshold(counts,i)
    if inter_class_variance>maximum:
        maximum=inter_class_variance
        otsu_thr=i


#Check answer using code from wikipedia
def compute_otsu_criteria(im, th):
    # create the thresholded image
    thresholded_im = np.zeros(im.shape)
    thresholded_im[im >= th] = 1

    # compute weights
    nb_pixels = im.size
    nb_pixels1 = np.count_nonzero(thresholded_im)
    weight1 = nb_pixels1 / nb_pixels
    weight0 = 1 - weight1

    # if one the classes is empty, eg all pixels are below or above the threshold, that threshold will not be considered
    # in the search for the best threshold
    if weight1 == 0 or weight0 == 0:
        return np.inf

    # find all pixels belonging to each class
    val_pixels1 = im[thresholded_im == 1]
    val_pixels0 = im[thresholded_im == 0]

    # compute variance of these classes
    var0 = np.var(val_pixels0) if len(val_pixels0) > 0 else 0
    var1 = np.var(val_pixels1) if len(val_pixels1) > 0 else 0

    return weight0 * var0 + weight1 * var1
im =A_slice # load your image as a numpy array.
threshold_range = range(np.max(im)+1)
criterias = [compute_otsu_criteria(im, th) for th in threshold_range] #seems to be slower than mine
best_threshold = threshold_range[np.argmin(criterias)] 
print(best_threshold==otsu_thr) #they are the same!

#Look at Otsu threshold image
A_slice_otsu_threshold=A_slice > otsu_thr #Makes all points above (below) this value black (white)
slice_otsu_threshold=im.fromarray(A_slice_otsu_threshold)
slice_otsu_threshold.save('slice_78_otsu_threshold.png') #Take a look at the image (saves to folder)

#Want to do this for every slice and scroll through the images
A_otsu=A #make a copy of A 
for slice in np.arange(0,nsls):
    A_slice=A[slice,:,:] #take slice
    counts,bins=np.histogram(A_slice,bins=256) #find counts
    maximum=OtsuThreshold(counts,1) #initialise maximum and threshold value
    otsu_thr=1
    for i in np.arange(2,257): #perform otsu thresholding algorithm
        inter_class_variance=OtsuThreshold(counts,i)
        if inter_class_variance>maximum:
            maximum=inter_class_variance
            otsu_thr=i
    A_otsu[slice,:,:]=A_otsu[slice,:,:] > otsu_thr #apply the threshold to the slice

# Take a look at the thresholded images
from utils.IndexTracker import *
fig, ax = plt.subplots(1, 1)
tracker = IndexTracker(ax,np.moveaxis(A_otsu,0,-1))
fig.canvas.mpl_connect('scroll_event', tracker.on_scroll)
plt.show()

#Aim to get rid of the cracks
isolated_row=A_otsu[77,480,:]
column_number=np.arange(1,1531)

#Plot showing the distribution of black and white values in this row
fig, ax = plt.subplots()
ax.plot(column_number, isolated_row)
ax.set(xlabel='Column', ylabel='Row 481 in image',
       title='Distribution of black (value 0) vs white (value 1) pixels')
ax.grid()
plt.show()

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

print(width_vector) #only two significant regions with black pixels, can probably set the rest to white (value 1)

indices=np.arange(1,len(width_vector))
tracker=0
for i in indices:
    if width_vector[i]<20:  #in this case, want to set these parts to white
        set_white=black_pixels2[tracker:tracker+width_vector[i]] #indices that need changing 
        isolated_row[set_white]=1 #change them to white
    tracker=tracker+width_vector[i]


#I've written this into a function
from utils.RemoveCracks import *
row_number=np.arange(1,1021)
im_number=77
new_image=np.zeros((1020,1530))
for i in row_number:
    isolated_row=A_otsu[im_number,i-1,:]
    crack_removed_row=RemoveCracks(isolated_row,20) #set the threshold width=20
    new_image[i-1,:]=crack_removed_row

#Look at the new image, some issue going on here, sort out another day, for some reason giving black images
new_image=(new_image==1)
slice_crack_removed=im.fromarray(new_image)
slice_crack_removed.save('slice_78_crack_threshold20.png') #Take a look at the image (saves to folder)















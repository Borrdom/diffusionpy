import cv2 
import numpy as np
from matplotlib import pyplot as plt
from scipy import ndimage
from skimage import io,color,measure
from tkinter.filedialog import askopenfilename
import pandas as pd
import os


def BinaryImage(pathtoimage="particles.jpg",cropx=50,cropy=False,threshold=70,sharpen=False):
    if cropx and not cropy:
        img=cv2.imread(pathtoimage,0)[:-cropx,:]
    elif cropx and cropy:
        img=cv2.imread(pathtoimage,0)[:-cropx,:-cropy]
    elif not cropx and cropy:
        img=cv2.imread(pathtoimage,0)[:,:-cropy]
    else:
        img=cv2.imread(pathtoimage,0)
    kernel= np.array([[-1,-1,-1], [-1,9,-1], [-1,-1,-1]])
    img=cv2.filter2D(img, -1, kernel) if sharpen else img
    cv2.imshow("Original_BlackandWhite",img)
    cv2.imwrite("Original_BlackandWhite.jpg",img)
    ret, thresh=cv2.threshold(img,threshold,255,cv2.THRESH_BINARY_INV)
    return thresh,img

def DilateOrErode(thresh,it=1,kernelsize=2,keyword=None):
    kernel=np.ones((kernelsize,kernelsize),np.uint8)
    eroded=cv2.erode(thresh,kernel,iterations=it)
    dilated=cv2.dilate(thresh,kernel,iterations=it)
    cv2.imshow("Binary",thresh)
    cv2.imwrite("Binary.jpg",thresh)
    cv2.imshow("Binary_eroded",eroded)
    cv2.imshow("Binary_dilated",dilated)
    if "eroded"==keyword:
        return eroded
    elif "dilated"==keyword:
        return dilated
    else:
        return thresh
    
def DetectLength(thresh,img,pix2mum,name="default"):
    mask=thresh==255
    s=[[1,1,1],[1,1,1],[1,1,1]]
    labeled_mask,num_labels=ndimage.label(mask,structure=s)
    img2=color.label2rgb(labeled_mask,bg_label=0)
    cv2.imshow("colored labels",img2)
    clusters=measure.regionprops(labeled_mask,img)
    Lengths=[]
    for prop in clusters:
        print('Label:{} Length:{}'.format(prop.label,prop.major_axis_length))
        Lengths.append(prop.major_axis_length) 
    Lengths=np.asarray(Lengths)
    MeanLength=np.mean(Lengths)/pix2mum
    NumberOfCrystals=len(Lengths)
    print("Overall "+str(NumberOfCrystals)+" were detected with an average size of "+str(MeanLength)+" \u03BCm")
    fig,ax=plt.subplots()
    n,bins, patches = ax.hist(x=Lengths,bins='auto',color='#0504aa',alpha=0.7,rwidth=0.85)
    ax.grid(axis='y', alpha=0.75)
    ax.set_ylabel('Anzahl')
    ax.set_xlabel('$L_{cr} $/$ \mu m$')
    ax.plot([MeanLength,MeanLength],[0, NumberOfCrystals/2],'k--')
    path=os.path.join(os.getcwd(),"CrystalDistribution",name)
    os.makedirs(path) if not os.path.isdir(path) else None
    filex=os.path.join(path,name+".xlsx")
    fileimg1=os.path.join(path,name+"_dist.jpg")
    fileimg2=os.path.join(path,name+"_detect.jpg")
    datawrite=pd.concat((pd.DataFrame(Lengths),pd.DataFrame(np.asarray([MeanLength])),pd.DataFrame(np.asarray([NumberOfCrystals]))),ignore_index=True, axis=1)
    datawrite.columns=["Length[micrometer]","Length_bar[micrometer]","Number[-]"]
    datawrite.to_excel(filex)
    fig.savefig(fileimg1)
    cv2.imwrite(fileimg2,img2*255)
    return NumberOfCrystals,MeanLength,Lengths

if __name__=="__main__":
    pathtoimage=askopenfilename()
    plt.close('all')
    #pathtoimage="image0621.jpg"
    pix2mum=145/50
    name=os.path.splitext(os.path.basename(pathtoimage))[0]

    thresh,img=BinaryImage(pathtoimage,cropx=50,cropy=False,threshold=70,sharpen=True)
    #thresh=DilateOrErode(thresh,it=1,kernelsize=2,keyword="eroded")
    NumberOfCrystals,MeanLength,Lengths=DetectLength(thresh,img,pix2mum,name)

    

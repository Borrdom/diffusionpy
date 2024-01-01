import tensorflow_hub as hub
import tensorflow as tf
from matplotlib import pyplot as plt
import numpy as np
import cv2

model=hub.load("https://tfhub.dev/google/magenta/arbitrary-image-stylization-v1-256/2")

def load_image(img_path): 
    img=tf.io.read_file(img_path)
    img=tf.image.decode_image(img,channels=3)
    img=tf.image.convert_image_dtype(img,tf.float32)
    img=img[tf.newaxis,:]
    return img

#style_path="VanGogh.jpg"

#style_path="Gustav-klimt-the-tree-of-life.jpg"
style_path="VanGogh2.jpg"
#content_path="Tabletten-1024x683.jpg"
#content_path="Splash7droplessmod7.jpg"

content_path="gud.jpg"
#style_path="Paper2Fig.jpg"
#content_path="photo-1550165703-d2ed566a8c4f.jpg"
#style_path="Splash4.jpg"

#style_path="Gustav-klimt-the-tree-of-life.jpg"
style_image=load_image(style_path)
content_image=load_image(content_path)



stylized_image=model(tf.constant(content_image),tf.constant(style_image))[0]

stylized_path=style_path.split(".")[0]+"_"+content_path
cv2.imwrite(stylized_path,cv2.cvtColor(np.squeeze(stylized_image)*255,cv2.COLOR_BGR2RGB))
plt.imshow(np.squeeze(stylized_image))
plt.show()

import numpy as np
import matplotlib.pyplot as plt
from skimage.color import rgb2gray
from skimage.io import imread
from skimage import measure
import glob


def FindContour(fname):
    im = imread(fname)
    img = rgb2gray(im)
    contours = measure.find_contours(img, 0.20)
    maxlen = 0
    for n, contour in enumerate(contours):
        if len(contour)>maxlen:
            maxlen = len(contour)
            shape = contour
    fig, ax = plt.subplots()
    ax.plot(shape[:,1],shape[:,0],'b-')
    print shape
    ax.imshow(im)
    ax.axis('image')
    ax.set_xticks([])
    ax.set_yticks([])
    fname = fname[:-4] + '-contour.png'
    plt.savefig(fname, bbox_inches='tight',pad_inches=0)
    plt.close(fig)
    return


def main():
    FindContour('shape001.png')
    # figs = glob.iglob('shape[0-9]*[0-9].png')
    # for fig in figs:
    #     FindContour(fig)
    return

if __name__ == "__main__":
    main()



#gPb-junctions

This is the code for my undergraduate thesis: Detecting Junctions in Photographs of Objects (April 2011), available at https://raw.github.com/hrichardlee/gPb-junctions/master/Lee2011-Junctions.pdf.

##Abstract
Line drawings are inherently simple, yet retain most of the information of a full image. This suggest that line drawings, also called contour images, could be useful as an intermediate representation for computer understanding of images. However, extracting physically meaningful contours from photographs has proved to be diffcult. Recently, there has been significant progress in detecting contours. As one example, the gPb algorithm (Arbeláez et al., 2010) uses local cues such as brightness, color, and texture, along with global cues that interpret the image as a whole in order to achieve high-performance contour detection.

However, contour detectors perform poorly near junctions. Junctions are points where two or more contours meet, creating two or more regions near a single point. Contour detectors, including gPb, assume that contours occur at the meeting of only two regions. We find that this assumption breaks down near junctions, resulting in poor performance. This is unfortunate because the location of junctions and the orientation of the incident edges at a junction are useful in reasoning about the geometric relationships between contours, and thus reasoning about the shape of the object.

These two observations, that junctions carry significant information and that con-tour detectors perform poorly near junctions, suggest that finding a complete and useful contour image requires complementing a contour detector with a junction detector. However, while there has been substantial recent progress in contour detection, there has not been comparable work and progress on the problem of junction detection.

We thus build and implement a junction detector. We are informed by the insight that junction points are a generalization of contour points, and we adapt ideas from the gPb contour detector in developing our algorithm. We also capture images of objects, and run our implementation on our dataset. Although our results are qualitative, they suggest that our junction detector could be useful as a complement to existing contour detection schemes.

#Notes
The paper and the code are an extension of the work in Contour Detection and Hierarchical Image Segmentation by Pablo Arbeláez, Michael Maire, Charless Fowlkes, and Jitendra Malik, published in CVPR 2011 (http://vision.caltech.edu/~mmaire/). The code here starts from the code used in that group's research (http://vision.caltech.edu/~mmaire/software/gpb_src.tar.gz and http://www.eecs.berkeley.edu/Research/Projects/CS/vision/grouping/resources.html)

Additionally, this work was done in April 2011. I have reconstructed the history of this code in June 2013 in order to preserve this for posterity.
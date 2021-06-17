Abdul Rehman (3440146) M.Sc. INFOTECH
An Trieu (3523966) M.Sc. INFOTECH
Khanh Quynh Nguyen (3517671) M.Sc. INFOTECH
Valentina Roldan (3519666) M.Sc. INFOTECH

Problem 3
b)
- The global shows the total frequency components of the image. We can see some diagonal lines, correspondings to the masts of the boats.
= The 8x8 DCT shows the localized frequency component of a 8x8 block. Hence we can still see the general shape of the boats. The grid-structure dots are the 0-frequency component of the 8x8 blocks.

c) The high frequencies have been removed. 8x8 DCT shows better result, since it does not have the ringing effect comparing to global DCT. This could be because the global DCT has large high-frequency coefficient, hence eliminating them affect the image badly.

d) Using the weighting matrix yields better result. This can be seen in the difference between the clouds. In the first strategy, the cloud has blocky edges.

Problem 4
input image:                      machine.pgm
available point transformations:
(0) affine rescaling
(1) gamma correction
your choice:                      0
smallest grey value:              50
largest  grey value:              255
output image:                     machine_rescaled.pgm

Problem 5
b) A low gamma (0.5) yields a better result than a high gamma (2). This is because the image has generally high intensity. To equalize the intensity, we need to increase contrast in the high intensity range. Thus a low gamma value.
Abdul Rehman (3440146) M.Sc. INFOTECH
An Trieu (3523966) M.Sc. INFOTECH
Khanh Quynh Nguyen (3517671) M.Sc. INFOTECH
Valentina Roldan (3519666) M.Sc. INFOTECH

Problem 3
b) filter parameter: 0.0005
c) filter parameter: 0.005
d) filter parameter: 0.01
Deblurring takes so long for hogblur.pgm because this image's size is 240x320 which is not a power of 2.
Hence FFT cannot be computed efficiently.
bus1.pgm and bus2.pgm both have size 512x512 which if optimal for FFT computation.
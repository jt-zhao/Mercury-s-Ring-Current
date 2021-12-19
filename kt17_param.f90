!*******************************************
! Save the following code as kt17_param.f90
!*******************************************
!
! *** Parameters v1.0 ***
!
! Model parameters from Korth et al., Modular model for Mercury’s magnetospheric
! magnetic field confined within the average observed magnetopause, J. Geophys.
! Res. Space Physics, 120, doi: 10.1002/2015JA021022, 2015.
       mu= 190.0d0 ! dipole moment [nT Rp^3]
       tilt= 0.0d0*3.14159265d0/180.0d0 ! dipole tilt angle [radians]
       pdip= 0.0d0*3.14159265d0/180.0d0 ! dipole longitude [radians]
       rdx= 0.0d0 ! dipole x offset
       rdy= 0.0d0 ! dipole y offset
       rdz= 0.196d0 ! dipole z offset
       rss= 1.41d0 ! distance from Mercury center to sub-
! solar magnetopause [Rp]
       r0= 1.42d0 ! distance from Mercury center to fitted
! sub-solar magnetopause [Rp]
       alfa= 0.5d0 ! magnetopause flaring factor
       tamp1= 7.64d0 ! tail disk current magnitude
       tamp2= 2.06d0 ! harris sheet current magntidue
       d0= 0.09d0 ! half-width of current sheet in Z at
! inner edge of tail current [Rp]
       deltadx= 1.0d0 ! expansion magnitudes of tail current
! sheet in x direction
       deltady= 0.1d0 ! expansion magnitudes of tail current
! sheet in y direction
       scalex= 1.5d0 ! e-folding distance for the sunward
! expansion of the harris sheet
       scaley= 9.0d0 ! scale distance for the flankward
! expansion of the harris sheet
       zshift= 3.5d0 ! location of image sheets for harris
! sheet
       mptol= 1.0d-3 ! Tolerance for magnetopause encounter
       r_taildisk=(/ &
59048.35734, &
-135664.4246, &
-913.4507339, &
209989.1008, &
-213142.9370, &
19.69235037, &
-18.16704312, &
12.69175932, &
-14.13692134, &
14.13449724, &
7.682813704, &
9.663177797, &
0.6465427021, &
1.274059603, &
1.280231032 &
/)

       r_dipshld=(/ &
7.792407683, &
74.37654983, &
4.119647072, &
-131.3308600, &
546.6006311, &
-1077.694401, &
52.46268495, &
1057.273707, &
-74.91550119, &
-141.8047123, &
3.876004886, &
156.2250932, &
-506.6470185, &
1439.804381, &
-64.55225925, &
-1443.754088, &
0.1412297078, &
0.7439847555, &
1.042798338, &
0.7057116022  &
/)

       r_diskshld=(/ &
-398.4670279, &
-1143.001682, &
-1836.300383, &
-73.92180417, &
-326.3986853, &
-29.96868107, &
-1157.035602, &
-604.1846034, &
-52.04876183, &
-2030.691236, &
-1529.120337, &
-6.382209946, &
2587.666032, &
213.8979183, &
-28.30225993, &
630.1309859, &
2968.552238, &
888.6328623, &
497.3863092, &
2304.254471, &
858.4176875, &
1226.958595, &
850.1684953, &
-20.90110940, &
-203.9184239, &
-792.6099018, &
1115.955690, &
527.3226825, &
22.47634041, &
-0.0704405637, &
-1405.093137, &
-97.20408343, &
5.656730182, &
-138.7129102, &
-1979.755673, &
5.407603749, &
1.091088905, &
0.6733299808, &
0.3266747827, &
0.9533161464, &
1.362763038, &
0.0014515208 &
/)

       r_slabshld=(/ &
-91.67686636, &
-87.31240824, &
251.8848107, &
95.65629983, &
-80.96810700, &
198.1447476, &
-283.1968987, &
-269.1514899, &
504.6322310, &
166.0272150, &
-214.9025413, &
623.7920115, &
-35.99544615, &
-322.8644690, &
345.7105790, &
928.8553184, &
810.1295090, &
19.62627762, &
-12.70326428, &
490.4662048, &
-814.0985363, &
-1781.184984, &
-1371.261326, &
60.31364790, &
116.6305510, &
-178.3347065, &
604.0308838, &
1155.151174, &
770.3896601, &
-202.8545948, &
298.6337705, &
304.7964641, &
33.70850254, &
393.6080147, &
308.1194271, &
-660.1691658, &
1.677629714, &
1.292226584, &
0.3116253398, &
-0.4392669057, &
0.7578074817, &
1.497779521 &
/)


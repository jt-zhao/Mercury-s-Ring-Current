!*******************************************
! Save the following code as kt17_common.f90
!*******************************************
!
!-----------------------------------------------------------------------
! KT17 Common Blocks
!-----------------------------------------------------------------------
!
       integer*4 n_taildisk,n_dipshld,n_diskshld,n_slabshld
       parameter (n_taildisk=15,n_dipshld=20,n_diskshld=42,n_slabshld=42)
       common /par/ mu,tilt,pdip,psun,rdx,rdy,rdz,rss,r0,alfa,tamp1,tamp2,d0, &
       deltadx,deltady,scalex,scaley,zshift,r_taildisk,r_dipshld,r_diskshld, &
       r_slabshld,mptol
       real*8 mu,tilt
       real*8 pdip,psun
       real*8 rdx,rdy,rdz
       real*8 rss,r0,alfa
       real*8 tamp1,tamp2
       real*8 d0,deltadx,deltady,scalex,scaley,zshift
       real*8 r_taildisk(n_taildisk),r_dipshld(n_dipshld)
       real*8 r_diskshld(n_diskshld),r_slabshld(n_slabshld)
       real*8 mptol
       real*8 rhel,act
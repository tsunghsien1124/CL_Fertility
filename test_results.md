
[80~66]
new (0807) - w/o Threads.@threads
19.259 ms (5765 allocations: 45.35 MiB)
19.295 ms (5765 allocations: 45.35 MiB)
new (0807) - w/ Threads.@threads for age 66~79
14.106 ms (9167 allocations: 46.04 MiB)
14.454 ms (9167 allocations: 46.04 MiB)
new (0807) - w/ Threads.@threads for age 66~79, ind_age
14.285 ms (9141 allocations: 45.83 MiB)
13.982 ms (9141 allocations: 45.83 MiB)
14.588 ms (9141 allocations: 45.83 MiB)
new (0807) - w/ Threads.@threads for age 66~80, ind_age
14.169 ms (9382 allocations: 45.86 MiB)
14.583 ms (9382 allocations: 45.86 MiB)
old (0806) - w/o Threads.@threads
25.503 ms (44974 allocations: 81.61 MiB)
25.388 ms (44974 allocations: 81.61 MiB)
old (0806) - w/ Threads.@threads
20.075 ms (38819 allocations: 77.71 MiB)
20.497 ms (38819 allocations: 77.71 MiB)
old (0806) - w/ Threads.@threads, ind_age
19.982 ms (38791 allocations: 77.48 MiB)
19.850 ms (38791 allocations: 77.48 MiB)

[80~65]
new (0807) - w/o Threads.@threads
8.099 s (5776 allocations: 45.35 MiB)
8.186 s (5775 allocations: 45.35 MiB)
new (0807) - w/o Threads.@threads, nested loop 
8.126 s (9393 allocations: 45.86 MiB)
new (0807) - w/ Threads.@threads (so nested loop) 
270.068 ms (9639 allocations: 46.37 MiB)
270.272 ms (9636 allocations: 46.14 MiB)
old (0806) - w/o Threads.@threads
578.308 ms (83028 allocations: 2.04 GiB)
583.879 ms (83028 allocations: 2.04 GiB)

[66]
new (0807) - w/o Threads.@threads
8.680 ms (5624 allocations: 45.34 MiB)
8.931 ms (5624 allocations: 45.34 MiB)
11.405 ms (5624 allocations: 45.34 MiB)
new (0807) - w/o Threads.@threads, nested loop
8.999 ms (5624 allocations: 45.34 MiB)
9.188 ms (5624 allocations: 45.34 MiB)
new (0807) - w Threads.@threads (so nested loop)
8.074 ms (5867 allocations: 45.39 MiB)
8.415 ms (5867 allocations: 45.39 MiB)
8.276 ms (5867 allocations: 45.39 MiB)
8.274 ms (5867 allocations: 45.39 MiB)

old (0806) - w/o Threads.@threads
14.124 ms (8433 allocations: 68.55 MiB)
13.920 ms (8433 allocations: 68.55 MiB)
old (0806) - w/ Threads.@threads
12.911 ms (8676 allocations: 68.60 MiB)
old (0806) - w/ Threads.@threads, w/ @views
12.945 ms (7976 allocations: 68.27 MiB)
13.417 ms (7976 allocations: 68.27 MiB)

[65]
new (0807)
8.086 s (5694 allocations: 45.35 MiB)
8.111 s (5624 allocations: 45.34 MiB)

new (0807) - put n-if-else within the innermost loop of a
7.986 s (5694 allocations: 45.35 MiB)
8.005 s (5694 allocations: 45.35 MiB)

new (0807) - put n-if-else within the innermost loop of a, shutting down the case of n > 0
8.987 ms (5624 allocations: 45.34 MiB)
8.970 ms (5624 allocations: 45.34 MiB)

new (0807) - put n-if-else within the innermost loop of a, shutting down the case of n > 0, nested loop for choice variables
8.106 s (5662 allocations: 45.35 MiB)
8.218 s (5673 allocations: 45.35 MiB)

old (0806) - w/ Threads.@threads
655.123 ms (43677 allocations: 2.03 GiB)
699.888 ms (43677 allocations: 2.03 GiB)
688.801 ms (43677 allocations: 2.03 GiB)
720.288 ms (43677 allocations: 2.03 GiB)
783.925 ms (43677 allocations: 2.03 GiB)
687.185 ms (43677 allocations: 2.03 GiB)

old (0806) - w/ Threads.@threads, Sol_all container outside loops
1.531 s (248283466 allocations: 4.16 GiB)
1.551 s (248283419 allocations: 4.16 GiB)

old (0806) - w/o Threads.@threads
8.585 s (43433 allocations: 2.03 GiB)

[80~46]
new (0807) - w/ Threads.@threads, old algo for all n
64.231 s (62745496 allocations: 65.95 GiB)
45.452 s (62737207 allocations: 65.94 GiB)
63.165 s (62737207 allocations: 65.94 GiB)

new (0807) - w/ Threads.@threads, old algo for all n, no n>0
349.457 ms (626206 allocations: 330.11 MiB)
347.709 ms (626206 allocations: 330.11 MiB)
344.385 ms (626206 allocations: 330.11 MiB)

new (0807) - w/ Threads.@threads, new algo for n=0
79.785 s (62125407 allocations: 65.67 GiB)
94.925 s (62125407 allocations: 65.67 GiB)

new (0807) - w/ Threads.@threads, new algo for n=0, no n>0
283.760 ms (14406 allocations: 46.80 MiB)
284.835 ms (14406 allocations: 46.80 MiB)
283.199 ms (14406 allocations: 46.80 MiB)
282.021 ms (14406 allocations: 46.80 MiB)
278.769 ms (14406 allocations: 46.80 MiB) [c break]
283.565 ms (14406 allocations: 46.80 MiB) [c break]
277.982 ms (14406 allocations: 46.80 MiB) [c break]

new (0807) - w/ Threads.@threads, new algo for n=0, no n>0, temp in func
315.504 ms (11087006 allocations: 215.76 MiB) [near operation]
278.139 ms (14406 allocations: 47.12 MiB) [outside solve function]

new (0807) - w/ Threads.@threads, new algo for all n
164.370 s (16624 allocations: 47.02 MiB)
127.070 s (14407 allocations: 46.87 MiB) [c & q break]
125.085 s (14407 allocations: 46.87 MiB) [c break, then q break]

new (0807) - w/ Threads.@threads, new algo for all n, temp in func
[near operation]
115.609 s (14407 allocations: 47.19 MiB) [outside solve function]

# ======= #
# macbook #
# ======= #

[80~64]
new (0807) - w/ Threads.@threads, new algo for all n
12.438 s (6776 allocations: 45.73 MiB)
24.673 s (3075608414 allocations: 45.87 GiB) [temp in func near operation]
12.504 s (6776 allocations: 45.73 MiB) [temp in func outside solve function]

[80~64] typo fixed!
new (0807) - w/ Threads.@threads, old algo for all n
1.594 s (3307976 allocations: 3.51 GiB)

new (0807) - w/ Threads.@threads, new algo for all n
2.687 s (6776 allocations: 45.73 MiB) [temp in func outside solve function]
2.723 s (6776 allocations: 45.73 MiB) [temp in func outside solve function]

new (0807) - w/ Threads.@threads, new algo for n=0, reshape
770.961 ms (6827 allocations: 45.73 MiB) [benchmark]
752.717 ms (142499 allocations: 57.85 MiB) [reshape, reduce+]
778.980 ms (142499 allocations: 57.85 MiB) [reshape, reduce+]
767.211 ms (142499 allocations: 57.85 MiB) [reshape, sum]

new (0807) - w/ Threads.@threads, new algo for all n, reshape
2.709 s (142499 allocations: 57.86 MiB) [benchmark]
5.058 s (98295140 allocations: 29.89 GiB) [reshape, sum]
2.988 s (29588291 allocations: 13.37 GiB) [reshape, sum, inbounds/views]
2.953 s (29510735 allocations: 13.36 GiB) [reshape, sum, inbounds/views]
2.768 s (9894207 allocations: 12.48 GiB) [outside reshape, sum, inbounds/views]
2.918 s (9856829 allocations: 12.48 GiB) [outside reshape, sum, inbounds/views]
2.918 s (9856829 allocations: 12.48 GiB) [outside reshape, sum, inbounds/views]
2.823 s (9841433 allocations: 12.48 GiB) [outside-most reshape, sum, inbounds/views]
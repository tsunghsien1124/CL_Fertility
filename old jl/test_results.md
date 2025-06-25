
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

446.508 ms (6892 allocations: 45.85 MiB)
456.918 ms (6712 allocations: 45.83 MiB)
458.209 ms (6712 allocations: 45.83 MiB)
449.036 ms (6841 allocations: 45.85 MiB)
451.049 ms (6841 allocations: 45.85 MiB)
907.349 ms (6841 allocations: 45.85 MiB)
898.166 ms (6841 allocations: 45.85 MiB)
897.428 ms (6841 allocations: 45.85 MiB) [80~64]
9.448 s (9165 allocations: 46.12 MiB) [80~46]
9.184 s (9272 allocations: 46.14 MiB)
9.081 s (9167 allocations: 46.14 MiB)

[c precomputation]
9.041 s (9272 allocations: 46.14 MiB)
9.129 s (12054 allocations: 46.33 MiB)
9.235 s (9167 allocations: 46.14 MiB)
9.147 s (9167 allocations: 46.14 MiB)
9.115 s (9167 allocations: 46.14 MiB)
9.233 s (9271 allocations: 46.14 MiB)
9.165 s (12144 allocations: 46.33 MiB)
9.161 s (9170 allocations: 46.72 MiB)

11.910 s (9880702 allocations: 10.24 GiB) [80~45]

[B202]
[80~46]
new (0809) - Benchmark
3.619 s (18998 allocations: 48.20 MiB)
3.576 s (18995 allocations: 40.80 MiB)
3.538 s (18995 allocations: 40.80 MiB)

[80~45] only n=0
3.551 s (19487 allocations: 40.86 MiB)
3.517 s (19487 allocations: 40.86 MiB)

[80~45] only n=0 and n=n_max
3.582 s (19487 allocations: 40.86 MiB)
3.612 s (19487 allocations: 40.86 MiB)

[80~45] only n=0 and n=n_max and other n with f_i = 2 (infertile)
3.789 s (19487 allocations: 40.86 MiB)
3.780 s (19487 allocations: 40.86 MiB)

[80~45] all n
3.902 s (19487 allocations: 40.86 MiB)
3.957 s (19487 allocations: 40.86 MiB)
4.205 s (19487 allocations: 40.86 MiB)
3.969 s (19487 allocations: 40.86 MiB)
3.917 s (19487 allocations: 40.86 MiB)

[80~44] 44 with only EV_inf
3.984 s (19982 allocations: 41.16 MiB)
3.951 s (19982 allocations: 41.16 MiB)

[80~18] 44~18 with only EV_inf
3.998 s (32776 allocations: 42.79 MiB)
3.948 s (32776 allocations: 42.79 MiB)

[80~18] 44~18 with EV_inf, n=0
4.054 s (32777 allocations: 42.89 MiB)
4.106 s (32776 allocations: 42.89 MiB)
4.022 s (32776 allocations: 42.89 MiB)
4.052 s (32776 allocations: 42.87 MiB)

[80~18] 44~18 with EV_inf, n=0, n=n_max
5.389 s (32777 allocations: 42.97 MiB)
5.383 s (32776 allocations: 42.97 MiB)
5.192 s (32776 allocations: 42.97 MiB)
5.190 s (32776 allocations: 42.97 MiB)

[80~18] 44~18 with EV_inf, n=0, n=n_max, other n with f_i = 2 (infertile)
9.339 s (32776 allocations: 42.97 MiB)
9.592 s (32776 allocations: 42.97 MiB)

[80~18] 44~18 with EV_inf, n=0, n=n_max, other n
13.661 s (32971 allocations: 42.98 MiB)
13.678 s (33165 allocations: 43.00 MiB)

[80~18] olddest algo [main_gridpoints.jl]
683.934 s (90982116783 allocations: 1623.86 GiB) [false]
866.745 s (514635325 allocations: 434.42 GiB) [true]

[macbook]
36.017 s (12787 allocations: 40.05 MiB)
35.024 s (12787 allocations: 40.05 MiB)
35.523 s (17559 allocations: 40.31 MiB)
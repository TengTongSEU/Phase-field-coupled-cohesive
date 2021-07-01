!**********************************************************************************************************
!                                                                                                         !      
!                                             UEL for PF-CZM                                              !
!                                        BFGS quasi-newton solver                                         !                                       
!                                                                                                         !
!**********************************************************************************************************
*  Copyright (C) 2019 South China University of Technology, China. All rights reserved.
*      
!*******************************************************************************************************F***
!
      module NumKind
!
!**********************************************************************************************************
        implicit none
        integer (kind(1)), parameter :: ikind = kind(1), 
     &                                  rkind = kind(0.D0), 
     &                                  lkind = kind(.true.)
!       
      end module Numkind
      

!**********************************************************************************************************
!
      module ComVarible
!
!**********************************************************************************************************
        use NumKind
        implicit none
        
        real(rkind) :: USRVAR(100000,18,4) 
        common USRVAR
!       
      end module ComVarible


        
        
!**********************************************************************************************************
!
      module ModelParam
!
!**********************************************************************************************************
        use NumKind
        use ComVarible
        implicit none

        ! Constants

        ! Flag of initilization
        logical (lkind) :: bInitialized = .false.

        ! Tolerance
        real    (rkind), parameter :: TOL = 1.0d-12 
        ! number of guass points
        integer (ikind), parameter :: ngp = 2
        
        ! geometric function parameter
        real(rkind), parameter :: c0 = 3.1415926535897932384626433832d0
        
        ! 
        real(rkind) :: thk, EA, nu, Gf, ft, lb
        real(rkind) :: De(3, 3)        
        real(rkind) :: p, a1, a2, a3
        real(rkind) :: gp(ngp), gw(ngp)
        real(rkind) :: QQ(12,12) 

        
        !
        contains

        !===================================
          subroutine Initialize(props, nprops, istype, jelem)
        !===================================

            integer (ikind), intent (in) :: nprops, istype, jelem
            real    (rkind), intent (in) :: props(nprops)

            !********************************************
            real(rkind) :: G0, K11, K12
            integer(ikind) :: indexq(12), i
            
            real(rkind) :: coeff
            
            integer (ikind), parameter :: N_ELEM=100000
            integer (ikind), parameter :: NSTV=18
        
        



            ! material properties
            EA   =  props(1) ! props(1) -- Young's modulus
            nu   =  props(2) ! props(2) -- Poisson's ratio
            ft   =  props(3) ! props(3) -- failure strength
            Gf   =  props(4) ! props(4) -- fracture energy
            lb   =  props(5) ! props(5) -- length scale            
            thk  =  props(6) ! props(6) -- thickness
            if (thk < TOL) thk = 1.0
            
            
                    coeff   = 1.0d0

        if (JELEM	.eq.	1	)		coeff = 	1.084736072
        if (JELEM	.eq.	2	)		coeff = 	0.828696135
        if (JELEM	.eq.	3	)		coeff = 	0.90269218
        if (JELEM	.eq.	4	)		coeff = 	0.889364444
        if (JELEM	.eq.	5	)		coeff = 	0.919059234
        if (JELEM	.eq.	6	)		coeff = 	0.959463177
        if (JELEM	.eq.	7	)		coeff = 	1.011367345
        if (JELEM	.eq.	8	)		coeff = 	0.915328319
        if (JELEM	.eq.	9	)		coeff = 	1.079087828
        if (JELEM	.eq.	10	)		coeff = 	1.066949179
        if (JELEM	.eq.	11	)		coeff = 	1.146492837
        if (JELEM	.eq.	12	)		coeff = 	0.717900186
        if (JELEM	.eq.	13	)		coeff = 	0.704401613
        if (JELEM	.eq.	14	)		coeff = 	1.076429305
        if (JELEM	.eq.	15	)		coeff = 	0.972999654
        if (JELEM	.eq.	16	)		coeff = 	0.917023433
        if (JELEM	.eq.	17	)		coeff = 	1.021516082
        if (JELEM	.eq.	18	)		coeff = 	0.880990444
        if (JELEM	.eq.	19	)		coeff = 	0.947503053
        if (JELEM	.eq.	20	)		coeff = 	0.983777177
        if (JELEM	.eq.	21	)		coeff = 	1.028135031
        if (JELEM	.eq.	22	)		coeff = 	0.881321934
        if (JELEM	.eq.	23	)		coeff = 	0.799844927
        if (JELEM	.eq.	24	)		coeff = 	0.891424599
        if (JELEM	.eq.	25	)		coeff = 	0.989448819
        if (JELEM	.eq.	26	)		coeff = 	0.759263596
        if (JELEM	.eq.	27	)		coeff = 	1.031590745
        if (JELEM	.eq.	28	)		coeff = 	0.954684931
        if (JELEM	.eq.	29	)		coeff = 	0.735459412
        if (JELEM	.eq.	30	)		coeff = 	1.027973714
        if (JELEM	.eq.	31	)		coeff = 	1.033177506
        if (JELEM	.eq.	32	)		coeff = 	0.771880865
        if (JELEM	.eq.	33	)		coeff = 	1.103593015
        if (JELEM	.eq.	34	)		coeff = 	1.018905196
        if (JELEM	.eq.	35	)		coeff = 	0.937587894
        if (JELEM	.eq.	36	)		coeff = 	1.047667979
        if (JELEM	.eq.	37	)		coeff = 	0.923814816
        if (JELEM	.eq.	38	)		coeff = 	0.861496511
        if (JELEM	.eq.	39	)		coeff = 	0.963538837
        if (JELEM	.eq.	40	)		coeff = 	0.918966713
        if (JELEM	.eq.	41	)		coeff = 	0.862648697
        if (JELEM	.eq.	42	)		coeff = 	1.038220741
        if (JELEM	.eq.	43	)		coeff = 	0.934777634
        if (JELEM	.eq.	44	)		coeff = 	1.08130412
        if (JELEM	.eq.	45	)		coeff = 	0.959605086
        if (JELEM	.eq.	46	)		coeff = 	0.840968029
        if (JELEM	.eq.	47	)		coeff = 	0.779215508
        if (JELEM	.eq.	48	)		coeff = 	0.964504887
        if (JELEM	.eq.	49	)		coeff = 	1.025116002
        if (JELEM	.eq.	50	)		coeff = 	0.918404481
        if (JELEM	.eq.	51	)		coeff = 	0.782730316
        if (JELEM	.eq.	52	)		coeff = 	0.961267009
        if (JELEM	.eq.	53	)		coeff = 	0.694612909
        if (JELEM	.eq.	54	)		coeff = 	1.049628241
        if (JELEM	.eq.	55	)		coeff = 	1.081868077
        if (JELEM	.eq.	56	)		coeff = 	1.019483091
        if (JELEM	.eq.	57	)		coeff = 	0.992264694
        if (JELEM	.eq.	58	)		coeff = 	0.985682144
        if (JELEM	.eq.	59	)		coeff = 	1.01551433
        if (JELEM	.eq.	60	)		coeff = 	0.904252766
        if (JELEM	.eq.	61	)		coeff = 	1.090917242
        if (JELEM	.eq.	62	)		coeff = 	0.990731071
        if (JELEM	.eq.	63	)		coeff = 	1.020098092
        if (JELEM	.eq.	64	)		coeff = 	1.016908504
        if (JELEM	.eq.	65	)		coeff = 	1.084397766
        if (JELEM	.eq.	66	)		coeff = 	0.936909466
        if (JELEM	.eq.	67	)		coeff = 	1.023645723
        if (JELEM	.eq.	68	)		coeff = 	1.064189541
        if (JELEM	.eq.	69	)		coeff = 	1.220344622
        if (JELEM	.eq.	70	)		coeff = 	1.023397087
        if (JELEM	.eq.	71	)		coeff = 	0.949622782
        if (JELEM	.eq.	72	)		coeff = 	0.820462714
        if (JELEM	.eq.	73	)		coeff = 	1.122082028
        if (JELEM	.eq.	74	)		coeff = 	0.794429521
        if (JELEM	.eq.	75	)		coeff = 	1.073515997
        if (JELEM	.eq.	76	)		coeff = 	0.843277665
        if (JELEM	.eq.	77	)		coeff = 	0.860487113
        if (JELEM	.eq.	78	)		coeff = 	0.782153984
        if (JELEM	.eq.	79	)		coeff = 	1.070997134
        if (JELEM	.eq.	80	)		coeff = 	0.962694025
        if (JELEM	.eq.	81	)		coeff = 	0.989954963
        if (JELEM	.eq.	82	)		coeff = 	1.057627242
        if (JELEM	.eq.	83	)		coeff = 	0.942492045
        if (JELEM	.eq.	84	)		coeff = 	0.933129186
        if (JELEM	.eq.	85	)		coeff = 	1.044103543
        if (JELEM	.eq.	86	)		coeff = 	0.958432525
        if (JELEM	.eq.	87	)		coeff = 	0.636260985
        if (JELEM	.eq.	88	)		coeff = 	0.966809257
        if (JELEM	.eq.	89	)		coeff = 	0.903885097
        if (JELEM	.eq.	90	)		coeff = 	0.988207076
        if (JELEM	.eq.	91	)		coeff = 	1.128803893
        if (JELEM	.eq.	92	)		coeff = 	1.020764163
        if (JELEM	.eq.	93	)		coeff = 	0.860016096
        if (JELEM	.eq.	94	)		coeff = 	1.005831252
        if (JELEM	.eq.	95	)		coeff = 	1.095303953
        if (JELEM	.eq.	96	)		coeff = 	0.960912277
        if (JELEM	.eq.	97	)		coeff = 	1.000284382
        if (JELEM	.eq.	98	)		coeff = 	0.887101221
        if (JELEM	.eq.	99	)		coeff = 	0.957075082
        if (JELEM	.eq.	100	)		coeff = 	0.858514423
        if (JELEM	.eq.	101	)		coeff = 	0.85228229
        if (JELEM	.eq.	102	)		coeff = 	1.05221071
        if (JELEM	.eq.	103	)		coeff = 	1.076498841
        if (JELEM	.eq.	104	)		coeff = 	0.850156187
        if (JELEM	.eq.	105	)		coeff = 	0.923167389
        if (JELEM	.eq.	106	)		coeff = 	1.152314122
        if (JELEM	.eq.	107	)		coeff = 	0.80183708
        if (JELEM	.eq.	108	)		coeff = 	0.95971483
        if (JELEM	.eq.	109	)		coeff = 	0.951435674
        if (JELEM	.eq.	110	)		coeff = 	0.933060203
        if (JELEM	.eq.	111	)		coeff = 	0.878523418
        if (JELEM	.eq.	112	)		coeff = 	0.83058006
        if (JELEM	.eq.	113	)		coeff = 	0.995932743
        if (JELEM	.eq.	114	)		coeff = 	1.094602224
        if (JELEM	.eq.	115	)		coeff = 	0.889313083
        if (JELEM	.eq.	116	)		coeff = 	1.009819759
        if (JELEM	.eq.	117	)		coeff = 	0.839880036
        if (JELEM	.eq.	118	)		coeff = 	0.998955655
        if (JELEM	.eq.	119	)		coeff = 	0.846305235
        if (JELEM	.eq.	120	)		coeff = 	1.056606177
        if (JELEM	.eq.	121	)		coeff = 	1.074104192
        if (JELEM	.eq.	122	)		coeff = 	0.814159226
        if (JELEM	.eq.	123	)		coeff = 	1.120572796
        if (JELEM	.eq.	124	)		coeff = 	0.906773843
        if (JELEM	.eq.	125	)		coeff = 	0.889352396
        if (JELEM	.eq.	126	)		coeff = 	0.981234927
        if (JELEM	.eq.	127	)		coeff = 	0.99675379
        if (JELEM	.eq.	128	)		coeff = 	0.67808804
        if (JELEM	.eq.	129	)		coeff = 	0.991568194
        if (JELEM	.eq.	130	)		coeff = 	0.980416957
        if (JELEM	.eq.	131	)		coeff = 	1.063607287
        if (JELEM	.eq.	132	)		coeff = 	1.011465283
        if (JELEM	.eq.	133	)		coeff = 	1.014792028
        if (JELEM	.eq.	134	)		coeff = 	0.803014614
        if (JELEM	.eq.	135	)		coeff = 	1.034091753
        if (JELEM	.eq.	136	)		coeff = 	1.015736983
        if (JELEM	.eq.	137	)		coeff = 	0.98889748
        if (JELEM	.eq.	138	)		coeff = 	0.899123046
        if (JELEM	.eq.	139	)		coeff = 	1.068539327
        if (JELEM	.eq.	140	)		coeff = 	0.820176165
        if (JELEM	.eq.	141	)		coeff = 	1.095388671
        if (JELEM	.eq.	142	)		coeff = 	0.974550505
        if (JELEM	.eq.	143	)		coeff = 	1.133256543
        if (JELEM	.eq.	144	)		coeff = 	0.881543257
        if (JELEM	.eq.	145	)		coeff = 	0.902027624
        if (JELEM	.eq.	146	)		coeff = 	1.044082514
        if (JELEM	.eq.	147	)		coeff = 	0.909152049
        if (JELEM	.eq.	148	)		coeff = 	0.947743225
        if (JELEM	.eq.	149	)		coeff = 	0.83344418
        if (JELEM	.eq.	150	)		coeff = 	0.947372773
        if (JELEM	.eq.	151	)		coeff = 	0.796982114
        if (JELEM	.eq.	152	)		coeff = 	0.98601713
        if (JELEM	.eq.	153	)		coeff = 	1.00265445
        if (JELEM	.eq.	154	)		coeff = 	0.967065757
        if (JELEM	.eq.	155	)		coeff = 	1.031432581
        if (JELEM	.eq.	156	)		coeff = 	0.770207361
        if (JELEM	.eq.	157	)		coeff = 	0.973178421
        if (JELEM	.eq.	158	)		coeff = 	1.032011226
        if (JELEM	.eq.	159	)		coeff = 	0.98284912
        if (JELEM	.eq.	160	)		coeff = 	0.901126439
        if (JELEM	.eq.	161	)		coeff = 	0.990665827
        if (JELEM	.eq.	162	)		coeff = 	1.054782591
        if (JELEM	.eq.	163	)		coeff = 	0.830016674
        if (JELEM	.eq.	164	)		coeff = 	0.939806919
        if (JELEM	.eq.	165	)		coeff = 	0.998451896
        if (JELEM	.eq.	166	)		coeff = 	1.041827002
        if (JELEM	.eq.	167	)		coeff = 	1.042671486
        if (JELEM	.eq.	168	)		coeff = 	0.957783702
        if (JELEM	.eq.	169	)		coeff = 	0.982251721
        if (JELEM	.eq.	170	)		coeff = 	0.886374244
        if (JELEM	.eq.	171	)		coeff = 	1.102481697
        if (JELEM	.eq.	172	)		coeff = 	0.835453007
        if (JELEM	.eq.	173	)		coeff = 	0.909160523
        if (JELEM	.eq.	174	)		coeff = 	1.071253657
        if (JELEM	.eq.	175	)		coeff = 	0.828664887
        if (JELEM	.eq.	176	)		coeff = 	1.048793064
        if (JELEM	.eq.	177	)		coeff = 	0.932783879
        if (JELEM	.eq.	178	)		coeff = 	0.951859956
        if (JELEM	.eq.	179	)		coeff = 	1.061605428
        if (JELEM	.eq.	180	)		coeff = 	1.178679441
        if (JELEM	.eq.	181	)		coeff = 	0.873782485
        if (JELEM	.eq.	182	)		coeff = 	0.87667565
        if (JELEM	.eq.	183	)		coeff = 	0.985597262
        if (JELEM	.eq.	184	)		coeff = 	1.11111169
        if (JELEM	.eq.	185	)		coeff = 	0.939343874
        if (JELEM	.eq.	186	)		coeff = 	1.057424519
        if (JELEM	.eq.	187	)		coeff = 	0.89136023
        if (JELEM	.eq.	188	)		coeff = 	0.954363688
        if (JELEM	.eq.	189	)		coeff = 	1.032291151
        if (JELEM	.eq.	190	)		coeff = 	0.782969292
        if (JELEM	.eq.	191	)		coeff = 	0.879465905
        if (JELEM	.eq.	192	)		coeff = 	0.808859513
        if (JELEM	.eq.	193	)		coeff = 	1.103463551
        if (JELEM	.eq.	194	)		coeff = 	1.054209399
        if (JELEM	.eq.	195	)		coeff = 	0.88805483
        if (JELEM	.eq.	196	)		coeff = 	0.903225725
        if (JELEM	.eq.	197	)		coeff = 	0.871356482
        if (JELEM	.eq.	198	)		coeff = 	0.963482487
        if (JELEM	.eq.	199	)		coeff = 	0.984398842
        if (JELEM	.eq.	200	)		coeff = 	0.931580443
        if (JELEM	.eq.	201	)		coeff = 	0.830315608
        if (JELEM	.eq.	202	)		coeff = 	0.912281119
        if (JELEM	.eq.	203	)		coeff = 	0.957399677
        if (JELEM	.eq.	204	)		coeff = 	1.019070312
        if (JELEM	.eq.	205	)		coeff = 	0.900602796
        if (JELEM	.eq.	206	)		coeff = 	0.996273511
        if (JELEM	.eq.	207	)		coeff = 	0.944704778
        if (JELEM	.eq.	208	)		coeff = 	0.808226044
        if (JELEM	.eq.	209	)		coeff = 	0.838032982
        if (JELEM	.eq.	210	)		coeff = 	0.799493532
        if (JELEM	.eq.	211	)		coeff = 	0.758385924
        if (JELEM	.eq.	212	)		coeff = 	0.853041424
        if (JELEM	.eq.	213	)		coeff = 	1.207721189
        if (JELEM	.eq.	214	)		coeff = 	1.191757691
        if (JELEM	.eq.	215	)		coeff = 	1.093140939
        if (JELEM	.eq.	216	)		coeff = 	1.030027869
        if (JELEM	.eq.	217	)		coeff = 	1.14225029
        if (JELEM	.eq.	218	)		coeff = 	0.984781107
        if (JELEM	.eq.	219	)		coeff = 	1.007322192
        if (JELEM	.eq.	220	)		coeff = 	0.95234507
        if (JELEM	.eq.	221	)		coeff = 	0.773532683
        if (JELEM	.eq.	222	)		coeff = 	1.019159531
        if (JELEM	.eq.	223	)		coeff = 	1.008141042
        if (JELEM	.eq.	224	)		coeff = 	0.827981253
        if (JELEM	.eq.	225	)		coeff = 	1.007483292
        if (JELEM	.eq.	226	)		coeff = 	1.070676798
        if (JELEM	.eq.	227	)		coeff = 	0.96183772
        if (JELEM	.eq.	228	)		coeff = 	0.829763189
        if (JELEM	.eq.	229	)		coeff = 	0.995541227
        if (JELEM	.eq.	230	)		coeff = 	0.903592645
        if (JELEM	.eq.	231	)		coeff = 	0.926359638
        if (JELEM	.eq.	232	)		coeff = 	0.977644569
        if (JELEM	.eq.	233	)		coeff = 	0.971298358
        if (JELEM	.eq.	234	)		coeff = 	0.743476125
        if (JELEM	.eq.	235	)		coeff = 	1.095207943
        if (JELEM	.eq.	236	)		coeff = 	1.024477716
        if (JELEM	.eq.	237	)		coeff = 	0.978564183
        if (JELEM	.eq.	238	)		coeff = 	0.938796327
        if (JELEM	.eq.	239	)		coeff = 	0.815773375
        if (JELEM	.eq.	240	)		coeff = 	0.972536798
        if (JELEM	.eq.	241	)		coeff = 	0.981168185
        if (JELEM	.eq.	242	)		coeff = 	0.884411554
        if (JELEM	.eq.	243	)		coeff = 	0.972862667
        if (JELEM	.eq.	244	)		coeff = 	0.827120823
        if (JELEM	.eq.	245	)		coeff = 	0.973245405
        if (JELEM	.eq.	246	)		coeff = 	0.964539186
        if (JELEM	.eq.	247	)		coeff = 	0.967495958
        if (JELEM	.eq.	248	)		coeff = 	1.039418486
        if (JELEM	.eq.	249	)		coeff = 	1.094130906
        if (JELEM	.eq.	250	)		coeff = 	1.104317252
        if (JELEM	.eq.	251	)		coeff = 	0.807819627
        if (JELEM	.eq.	252	)		coeff = 	1.038283757
        if (JELEM	.eq.	253	)		coeff = 	0.826627048
        if (JELEM	.eq.	254	)		coeff = 	0.897747519
        if (JELEM	.eq.	255	)		coeff = 	0.819152448
        if (JELEM	.eq.	256	)		coeff = 	0.75966158
        if (JELEM	.eq.	257	)		coeff = 	1.070064352
        if (JELEM	.eq.	258	)		coeff = 	0.992946031
        if (JELEM	.eq.	259	)		coeff = 	0.674976784
        if (JELEM	.eq.	260	)		coeff = 	0.920935198
        if (JELEM	.eq.	261	)		coeff = 	0.801465639
        if (JELEM	.eq.	262	)		coeff = 	0.968930444
        if (JELEM	.eq.	263	)		coeff = 	1.15601129
        if (JELEM	.eq.	264	)		coeff = 	0.927946424
        if (JELEM	.eq.	265	)		coeff = 	1.038918629
        if (JELEM	.eq.	266	)		coeff = 	0.956312387
        if (JELEM	.eq.	267	)		coeff = 	0.892752119
        if (JELEM	.eq.	268	)		coeff = 	0.932760511
        if (JELEM	.eq.	269	)		coeff = 	0.938557289
        if (JELEM	.eq.	270	)		coeff = 	0.982248434
        if (JELEM	.eq.	271	)		coeff = 	1.034947157
        if (JELEM	.eq.	272	)		coeff = 	0.983460887
        if (JELEM	.eq.	273	)		coeff = 	1.164550611
        if (JELEM	.eq.	274	)		coeff = 	0.932327132
        if (JELEM	.eq.	275	)		coeff = 	0.729927169
        if (JELEM	.eq.	276	)		coeff = 	1.089160289
        if (JELEM	.eq.	277	)		coeff = 	1.128024078
        if (JELEM	.eq.	278	)		coeff = 	0.809462969
        if (JELEM	.eq.	279	)		coeff = 	1.034116548
        if (JELEM	.eq.	280	)		coeff = 	1.167869526
        if (JELEM	.eq.	281	)		coeff = 	0.853305131
        if (JELEM	.eq.	282	)		coeff = 	1.069749551
        if (JELEM	.eq.	283	)		coeff = 	0.81417162
        if (JELEM	.eq.	284	)		coeff = 	1.089194172
        if (JELEM	.eq.	285	)		coeff = 	1.004174133
        if (JELEM	.eq.	286	)		coeff = 	0.937032659
        if (JELEM	.eq.	287	)		coeff = 	0.939513326
        if (JELEM	.eq.	288	)		coeff = 	0.913339153
        if (JELEM	.eq.	289	)		coeff = 	0.919880116
        if (JELEM	.eq.	290	)		coeff = 	0.982266263
        if (JELEM	.eq.	291	)		coeff = 	1.070037315
        if (JELEM	.eq.	292	)		coeff = 	0.882067166
        if (JELEM	.eq.	293	)		coeff = 	1.035665779
        if (JELEM	.eq.	294	)		coeff = 	0.919082661
        if (JELEM	.eq.	295	)		coeff = 	0.829313039
        if (JELEM	.eq.	296	)		coeff = 	1.094748337
        if (JELEM	.eq.	297	)		coeff = 	0.700149097
        if (JELEM	.eq.	298	)		coeff = 	1.132128084
        if (JELEM	.eq.	299	)		coeff = 	0.842341953
        if (JELEM	.eq.	300	)		coeff = 	0.842168691
        if (JELEM	.eq.	301	)		coeff = 	1.116070568
        if (JELEM	.eq.	302	)		coeff = 	0.95104327
        if (JELEM	.eq.	303	)		coeff = 	0.752873825
        if (JELEM	.eq.	304	)		coeff = 	1.01272964
        if (JELEM	.eq.	305	)		coeff = 	0.85755609
        if (JELEM	.eq.	306	)		coeff = 	0.934606606
        if (JELEM	.eq.	307	)		coeff = 	0.865622462
        if (JELEM	.eq.	308	)		coeff = 	0.861105473
        if (JELEM	.eq.	309	)		coeff = 	1.116289553
        if (JELEM	.eq.	310	)		coeff = 	1.02351073
        if (JELEM	.eq.	311	)		coeff = 	0.918095364
        if (JELEM	.eq.	312	)		coeff = 	0.966881111
        if (JELEM	.eq.	313	)		coeff = 	0.698183675
        if (JELEM	.eq.	314	)		coeff = 	0.883475444
        if (JELEM	.eq.	315	)		coeff = 	0.944658878
        if (JELEM	.eq.	316	)		coeff = 	1.019028317
        if (JELEM	.eq.	317	)		coeff = 	1.031388927
        if (JELEM	.eq.	318	)		coeff = 	0.809214097
        if (JELEM	.eq.	319	)		coeff = 	0.978620311
        if (JELEM	.eq.	320	)		coeff = 	0.852758026
        if (JELEM	.eq.	321	)		coeff = 	1.087768965
        if (JELEM	.eq.	322	)		coeff = 	0.827909121
        if (JELEM	.eq.	323	)		coeff = 	1.13630784
        if (JELEM	.eq.	324	)		coeff = 	0.79919473
        if (JELEM	.eq.	325	)		coeff = 	0.798541186
        if (JELEM	.eq.	326	)		coeff = 	0.957254404
        if (JELEM	.eq.	327	)		coeff = 	1.077965705
        if (JELEM	.eq.	328	)		coeff = 	1.056166729
        if (JELEM	.eq.	329	)		coeff = 	0.899821537
        if (JELEM	.eq.	330	)		coeff = 	0.84458811
        if (JELEM	.eq.	331	)		coeff = 	1.128748765
        if (JELEM	.eq.	332	)		coeff = 	0.879614927
        if (JELEM	.eq.	333	)		coeff = 	0.731445176
        if (JELEM	.eq.	334	)		coeff = 	1.006826662
        if (JELEM	.eq.	335	)		coeff = 	0.923055952
        if (JELEM	.eq.	336	)		coeff = 	1.006789547
        if (JELEM	.eq.	337	)		coeff = 	1.043464236
        if (JELEM	.eq.	338	)		coeff = 	0.867200103
        if (JELEM	.eq.	339	)		coeff = 	0.893478692
        if (JELEM	.eq.	340	)		coeff = 	1.024762458
        if (JELEM	.eq.	341	)		coeff = 	0.94033606
        if (JELEM	.eq.	342	)		coeff = 	0.985614465
        if (JELEM	.eq.	343	)		coeff = 	1.090819385
        if (JELEM	.eq.	344	)		coeff = 	1.140654597
        if (JELEM	.eq.	345	)		coeff = 	0.966468781
        if (JELEM	.eq.	346	)		coeff = 	1.024927077
        if (JELEM	.eq.	347	)		coeff = 	1.007681385
        if (JELEM	.eq.	348	)		coeff = 	1.022324064
        if (JELEM	.eq.	349	)		coeff = 	1.058559423
        if (JELEM	.eq.	350	)		coeff = 	0.991495195
        if (JELEM	.eq.	351	)		coeff = 	0.90289007
        if (JELEM	.eq.	352	)		coeff = 	1.047545981
        if (JELEM	.eq.	353	)		coeff = 	0.913797537
        if (JELEM	.eq.	354	)		coeff = 	0.979635194
        if (JELEM	.eq.	355	)		coeff = 	0.982286268
        if (JELEM	.eq.	356	)		coeff = 	1.057048444
        if (JELEM	.eq.	357	)		coeff = 	1.050970778
        if (JELEM	.eq.	358	)		coeff = 	0.929969685
        if (JELEM	.eq.	359	)		coeff = 	1.027605974
        if (JELEM	.eq.	360	)		coeff = 	0.947048619
        if (JELEM	.eq.	361	)		coeff = 	0.750644033
        if (JELEM	.eq.	362	)		coeff = 	0.896729214
        if (JELEM	.eq.	363	)		coeff = 	0.909389318
        if (JELEM	.eq.	364	)		coeff = 	0.727395672
        if (JELEM	.eq.	365	)		coeff = 	0.87209163
        if (JELEM	.eq.	366	)		coeff = 	0.932664642
        if (JELEM	.eq.	367	)		coeff = 	0.746023887
        if (JELEM	.eq.	368	)		coeff = 	1.109214897
        if (JELEM	.eq.	369	)		coeff = 	1.027692763
        if (JELEM	.eq.	370	)		coeff = 	0.649742293
        if (JELEM	.eq.	371	)		coeff = 	0.873467107
        if (JELEM	.eq.	372	)		coeff = 	0.970821706
        if (JELEM	.eq.	373	)		coeff = 	0.908806894
        if (JELEM	.eq.	374	)		coeff = 	0.986715942
        if (JELEM	.eq.	375	)		coeff = 	0.996669677
        if (JELEM	.eq.	376	)		coeff = 	1.044475639
        if (JELEM	.eq.	377	)		coeff = 	0.995913901
        if (JELEM	.eq.	378	)		coeff = 	1.134033316
        if (JELEM	.eq.	379	)		coeff = 	0.971654133
        if (JELEM	.eq.	380	)		coeff = 	1.009438299
        if (JELEM	.eq.	381	)		coeff = 	0.689971158
        if (JELEM	.eq.	382	)		coeff = 	0.948280608
        if (JELEM	.eq.	383	)		coeff = 	0.836064277
        if (JELEM	.eq.	384	)		coeff = 	0.989115287
        if (JELEM	.eq.	385	)		coeff = 	0.974473937
        if (JELEM	.eq.	386	)		coeff = 	0.847335473
        if (JELEM	.eq.	387	)		coeff = 	0.623195339
        if (JELEM	.eq.	388	)		coeff = 	0.957299171
        if (JELEM	.eq.	389	)		coeff = 	0.774326631
        if (JELEM	.eq.	390	)		coeff = 	0.887283014
        if (JELEM	.eq.	391	)		coeff = 	0.9447797
        if (JELEM	.eq.	392	)		coeff = 	0.70816346
        if (JELEM	.eq.	393	)		coeff = 	0.848302934
        if (JELEM	.eq.	394	)		coeff = 	0.726976484
        if (JELEM	.eq.	395	)		coeff = 	0.920429407
        if (JELEM	.eq.	396	)		coeff = 	0.996825747
        if (JELEM	.eq.	397	)		coeff = 	0.970486478
        if (JELEM	.eq.	398	)		coeff = 	0.787934471
        if (JELEM	.eq.	399	)		coeff = 	1.154577282
        if (JELEM	.eq.	400	)		coeff = 	1.063649045
        if (JELEM	.eq.	401	)		coeff = 	0.971860668
        if (JELEM	.eq.	402	)		coeff = 	0.951875254
        if (JELEM	.eq.	403	)		coeff = 	1.109182314
        if (JELEM	.eq.	404	)		coeff = 	0.916579187
        if (JELEM	.eq.	405	)		coeff = 	0.806859046
        if (JELEM	.eq.	406	)		coeff = 	1.08255861
        if (JELEM	.eq.	407	)		coeff = 	0.981071976
        if (JELEM	.eq.	408	)		coeff = 	1.024362418
        if (JELEM	.eq.	409	)		coeff = 	0.656466272
        if (JELEM	.eq.	410	)		coeff = 	0.932346493
        if (JELEM	.eq.	411	)		coeff = 	1.032087921
        if (JELEM	.eq.	412	)		coeff = 	1.072856755
        if (JELEM	.eq.	413	)		coeff = 	0.951296817
        if (JELEM	.eq.	414	)		coeff = 	0.846528031
        if (JELEM	.eq.	415	)		coeff = 	0.841440447
        if (JELEM	.eq.	416	)		coeff = 	0.843489002
        if (JELEM	.eq.	417	)		coeff = 	1.047533376
        if (JELEM	.eq.	418	)		coeff = 	0.951465636
        if (JELEM	.eq.	419	)		coeff = 	0.817672937
        if (JELEM	.eq.	420	)		coeff = 	1.077627465
        if (JELEM	.eq.	421	)		coeff = 	0.829954412
        if (JELEM	.eq.	422	)		coeff = 	0.798675318
        if (JELEM	.eq.	423	)		coeff = 	1.043025801
        if (JELEM	.eq.	424	)		coeff = 	1.098741565
        if (JELEM	.eq.	425	)		coeff = 	0.971134455
        if (JELEM	.eq.	426	)		coeff = 	0.842553182
        if (JELEM	.eq.	427	)		coeff = 	0.972456547
        if (JELEM	.eq.	428	)		coeff = 	0.98757201
        if (JELEM	.eq.	429	)		coeff = 	0.963247716
        if (JELEM	.eq.	430	)		coeff = 	1.07577446
        if (JELEM	.eq.	431	)		coeff = 	1.07298875
        if (JELEM	.eq.	432	)		coeff = 	0.820753228
        if (JELEM	.eq.	433	)		coeff = 	0.93413199
        if (JELEM	.eq.	434	)		coeff = 	1.028687377
        if (JELEM	.eq.	435	)		coeff = 	0.824542069
        if (JELEM	.eq.	436	)		coeff = 	1.11023456
        if (JELEM	.eq.	437	)		coeff = 	0.975637878
        if (JELEM	.eq.	438	)		coeff = 	0.893813827
        if (JELEM	.eq.	439	)		coeff = 	1.00788968
        if (JELEM	.eq.	440	)		coeff = 	0.990966025
        if (JELEM	.eq.	441	)		coeff = 	0.956441098
        if (JELEM	.eq.	442	)		coeff = 	0.80322773
        if (JELEM	.eq.	443	)		coeff = 	0.870747209
        if (JELEM	.eq.	444	)		coeff = 	1.103123257
        if (JELEM	.eq.	445	)		coeff = 	1.024778106
        if (JELEM	.eq.	446	)		coeff = 	0.996880957
        if (JELEM	.eq.	447	)		coeff = 	0.82462739
        if (JELEM	.eq.	448	)		coeff = 	0.985897353
        if (JELEM	.eq.	449	)		coeff = 	1.03624268
        if (JELEM	.eq.	450	)		coeff = 	0.935748191
        if (JELEM	.eq.	451	)		coeff = 	0.969708464
        if (JELEM	.eq.	452	)		coeff = 	0.799702615
        if (JELEM	.eq.	453	)		coeff = 	0.763735711
        if (JELEM	.eq.	454	)		coeff = 	0.851777899
        if (JELEM	.eq.	455	)		coeff = 	0.898792411
        if (JELEM	.eq.	456	)		coeff = 	0.885606692
        if (JELEM	.eq.	457	)		coeff = 	0.798730974
        if (JELEM	.eq.	458	)		coeff = 	1.105633684
        if (JELEM	.eq.	459	)		coeff = 	1.00873905
        if (JELEM	.eq.	460	)		coeff = 	1.184589851
        if (JELEM	.eq.	461	)		coeff = 	0.846370094
        if (JELEM	.eq.	462	)		coeff = 	0.961939331
        if (JELEM	.eq.	463	)		coeff = 	1.000466841
        if (JELEM	.eq.	464	)		coeff = 	1.040291163
        if (JELEM	.eq.	465	)		coeff = 	0.954209583
        if (JELEM	.eq.	466	)		coeff = 	1.021715513
        if (JELEM	.eq.	467	)		coeff = 	1.10372012
        if (JELEM	.eq.	468	)		coeff = 	1.094434032
        if (JELEM	.eq.	469	)		coeff = 	1.103739035
        if (JELEM	.eq.	470	)		coeff = 	0.988641221
        if (JELEM	.eq.	471	)		coeff = 	1.076628889
        if (JELEM	.eq.	472	)		coeff = 	0.979648638
        if (JELEM	.eq.	473	)		coeff = 	0.799375762
        if (JELEM	.eq.	474	)		coeff = 	1.003879435
        if (JELEM	.eq.	475	)		coeff = 	1.077974761
        if (JELEM	.eq.	476	)		coeff = 	0.94428533
        if (JELEM	.eq.	477	)		coeff = 	0.817609957
        if (JELEM	.eq.	478	)		coeff = 	1.005263745
        if (JELEM	.eq.	479	)		coeff = 	1.122361817
        if (JELEM	.eq.	480	)		coeff = 	1.069038106
        if (JELEM	.eq.	481	)		coeff = 	1.098956915
        if (JELEM	.eq.	482	)		coeff = 	0.886684441
        if (JELEM	.eq.	483	)		coeff = 	0.975972958
        if (JELEM	.eq.	484	)		coeff = 	0.913146855
        if (JELEM	.eq.	485	)		coeff = 	0.902315203
        if (JELEM	.eq.	486	)		coeff = 	0.943623367
        if (JELEM	.eq.	487	)		coeff = 	0.926109308
        if (JELEM	.eq.	488	)		coeff = 	0.81568091
        if (JELEM	.eq.	489	)		coeff = 	0.91513539
        if (JELEM	.eq.	490	)		coeff = 	0.817343027
        if (JELEM	.eq.	491	)		coeff = 	0.972970369
        if (JELEM	.eq.	492	)		coeff = 	1.069425508
        if (JELEM	.eq.	493	)		coeff = 	1.103862556
        if (JELEM	.eq.	494	)		coeff = 	0.896814763
        if (JELEM	.eq.	495	)		coeff = 	1.01649091
        if (JELEM	.eq.	496	)		coeff = 	0.912156121
        if (JELEM	.eq.	497	)		coeff = 	0.918437158
        if (JELEM	.eq.	498	)		coeff = 	0.955286407
        if (JELEM	.eq.	499	)		coeff = 	0.89649184
        if (JELEM	.eq.	500	)		coeff = 	0.962672181
        if (JELEM	.eq.	501	)		coeff = 	0.967338851
        if (JELEM	.eq.	502	)		coeff = 	0.964613136
        if (JELEM	.eq.	503	)		coeff = 	0.762190333
        if (JELEM	.eq.	504	)		coeff = 	0.994195938
        if (JELEM	.eq.	505	)		coeff = 	1.079268229
        if (JELEM	.eq.	506	)		coeff = 	1.03607751
        if (JELEM	.eq.	507	)		coeff = 	0.907409141
        if (JELEM	.eq.	508	)		coeff = 	0.84016368
        if (JELEM	.eq.	509	)		coeff = 	0.704950714
        if (JELEM	.eq.	510	)		coeff = 	1.0438748
        if (JELEM	.eq.	511	)		coeff = 	0.87856259
        if (JELEM	.eq.	512	)		coeff = 	0.93983636
        if (JELEM	.eq.	513	)		coeff = 	0.990497692
        if (JELEM	.eq.	514	)		coeff = 	0.961210782
        if (JELEM	.eq.	515	)		coeff = 	0.965222269
        if (JELEM	.eq.	516	)		coeff = 	0.918783782
        if (JELEM	.eq.	517	)		coeff = 	0.885411591
        if (JELEM	.eq.	518	)		coeff = 	1.01818551
        if (JELEM	.eq.	519	)		coeff = 	1.092050025
        if (JELEM	.eq.	520	)		coeff = 	0.84751874
        if (JELEM	.eq.	521	)		coeff = 	0.994113288
        if (JELEM	.eq.	522	)		coeff = 	0.872104651
        if (JELEM	.eq.	523	)		coeff = 	1.055626155
        if (JELEM	.eq.	524	)		coeff = 	1.082679926
        if (JELEM	.eq.	525	)		coeff = 	0.794046109
        if (JELEM	.eq.	526	)		coeff = 	0.816683225
        if (JELEM	.eq.	527	)		coeff = 	0.428245516
        if (JELEM	.eq.	528	)		coeff = 	0.824901116
        if (JELEM	.eq.	529	)		coeff = 	1.126817046
        if (JELEM	.eq.	530	)		coeff = 	0.951389533
        if (JELEM	.eq.	531	)		coeff = 	0.546737826
        if (JELEM	.eq.	532	)		coeff = 	0.960939788
        if (JELEM	.eq.	533	)		coeff = 	0.818668982
        if (JELEM	.eq.	534	)		coeff = 	1.102616433
        if (JELEM	.eq.	535	)		coeff = 	0.645398892
        if (JELEM	.eq.	536	)		coeff = 	0.777222481
        if (JELEM	.eq.	537	)		coeff = 	0.945709405
        if (JELEM	.eq.	538	)		coeff = 	0.982783109
        if (JELEM	.eq.	539	)		coeff = 	1.008207101
        if (JELEM	.eq.	540	)		coeff = 	0.894377995
        if (JELEM	.eq.	541	)		coeff = 	1.156863821
        if (JELEM	.eq.	542	)		coeff = 	0.99832123
        if (JELEM	.eq.	543	)		coeff = 	0.777252017
        if (JELEM	.eq.	544	)		coeff = 	0.950872432
        if (JELEM	.eq.	545	)		coeff = 	0.971224557
        if (JELEM	.eq.	546	)		coeff = 	0.964974416
        if (JELEM	.eq.	547	)		coeff = 	1.016214958
        if (JELEM	.eq.	548	)		coeff = 	0.741742953
        if (JELEM	.eq.	549	)		coeff = 	0.669788113
        if (JELEM	.eq.	550	)		coeff = 	0.960223543
        if (JELEM	.eq.	551	)		coeff = 	0.612553878
        if (JELEM	.eq.	552	)		coeff = 	0.976158724
        if (JELEM	.eq.	553	)		coeff = 	0.984251536
        if (JELEM	.eq.	554	)		coeff = 	1.044493265
        if (JELEM	.eq.	555	)		coeff = 	1.050955172
        if (JELEM	.eq.	556	)		coeff = 	0.843818558
        if (JELEM	.eq.	557	)		coeff = 	0.8921209
        if (JELEM	.eq.	558	)		coeff = 	0.955656603
        if (JELEM	.eq.	559	)		coeff = 	0.8458351
        if (JELEM	.eq.	560	)		coeff = 	0.960697077
        if (JELEM	.eq.	561	)		coeff = 	0.949289867
        if (JELEM	.eq.	562	)		coeff = 	1.044471041
        if (JELEM	.eq.	563	)		coeff = 	0.938721391
        if (JELEM	.eq.	564	)		coeff = 	1.068875083
        if (JELEM	.eq.	565	)		coeff = 	1.114332773
        if (JELEM	.eq.	566	)		coeff = 	0.907968537
        if (JELEM	.eq.	567	)		coeff = 	0.932407961
        if (JELEM	.eq.	568	)		coeff = 	1.042464132
        if (JELEM	.eq.	569	)		coeff = 	0.989597605
        if (JELEM	.eq.	570	)		coeff = 	0.925723448
        if (JELEM	.eq.	571	)		coeff = 	0.948320069
        if (JELEM	.eq.	572	)		coeff = 	1.074894983
        if (JELEM	.eq.	573	)		coeff = 	1.059158636
        if (JELEM	.eq.	574	)		coeff = 	1.213213833
        if (JELEM	.eq.	575	)		coeff = 	0.986383803
        if (JELEM	.eq.	576	)		coeff = 	0.967204466
        if (JELEM	.eq.	577	)		coeff = 	1.062492221
        if (JELEM	.eq.	578	)		coeff = 	0.913624471
        if (JELEM	.eq.	579	)		coeff = 	1.149323359
        if (JELEM	.eq.	580	)		coeff = 	1.078187862
        if (JELEM	.eq.	581	)		coeff = 	0.739727115
        if (JELEM	.eq.	582	)		coeff = 	0.689944646
        if (JELEM	.eq.	583	)		coeff = 	1.132699169
        if (JELEM	.eq.	584	)		coeff = 	0.965710995
        if (JELEM	.eq.	585	)		coeff = 	0.825899224
        if (JELEM	.eq.	586	)		coeff = 	1.035337374
        if (JELEM	.eq.	587	)		coeff = 	0.842976248
        if (JELEM	.eq.	588	)		coeff = 	0.853983237
        if (JELEM	.eq.	589	)		coeff = 	0.926020494
        if (JELEM	.eq.	590	)		coeff = 	1.198251049
        if (JELEM	.eq.	591	)		coeff = 	0.996784959
        if (JELEM	.eq.	592	)		coeff = 	0.794706508
        if (JELEM	.eq.	593	)		coeff = 	0.908997059
        if (JELEM	.eq.	594	)		coeff = 	0.997029183
        if (JELEM	.eq.	595	)		coeff = 	0.925071623
        if (JELEM	.eq.	596	)		coeff = 	1.035218923
        if (JELEM	.eq.	597	)		coeff = 	0.943606885
        if (JELEM	.eq.	598	)		coeff = 	0.670774552
        if (JELEM	.eq.	599	)		coeff = 	0.834046506
        if (JELEM	.eq.	600	)		coeff = 	1.023440396
        if (JELEM	.eq.	601	)		coeff = 	0.908264969
        if (JELEM	.eq.	602	)		coeff = 	1.002655492
        if (JELEM	.eq.	603	)		coeff = 	0.648520532
        if (JELEM	.eq.	604	)		coeff = 	1.094950239
        if (JELEM	.eq.	605	)		coeff = 	1.03311156
        if (JELEM	.eq.	606	)		coeff = 	0.855120029
        if (JELEM	.eq.	607	)		coeff = 	1.094720722
        if (JELEM	.eq.	608	)		coeff = 	0.9552281
        if (JELEM	.eq.	609	)		coeff = 	0.860412416
        if (JELEM	.eq.	610	)		coeff = 	0.887361844
        if (JELEM	.eq.	611	)		coeff = 	1.06929915
        if (JELEM	.eq.	612	)		coeff = 	0.98103911
        if (JELEM	.eq.	613	)		coeff = 	1.004769897
        if (JELEM	.eq.	614	)		coeff = 	0.969957922
        if (JELEM	.eq.	615	)		coeff = 	0.93884676
        if (JELEM	.eq.	616	)		coeff = 	1.067711475
        if (JELEM	.eq.	617	)		coeff = 	0.793894332
        if (JELEM	.eq.	618	)		coeff = 	0.922426757
        if (JELEM	.eq.	619	)		coeff = 	1.061384464
        if (JELEM	.eq.	620	)		coeff = 	0.945225831
        if (JELEM	.eq.	621	)		coeff = 	0.767404663
        if (JELEM	.eq.	622	)		coeff = 	0.868613448
        if (JELEM	.eq.	623	)		coeff = 	0.907136358
        if (JELEM	.eq.	624	)		coeff = 	0.973325494
        if (JELEM	.eq.	625	)		coeff = 	1.030148173
        if (JELEM	.eq.	626	)		coeff = 	0.944239069
        if (JELEM	.eq.	627	)		coeff = 	1.033569538
        if (JELEM	.eq.	628	)		coeff = 	1.01333286
        if (JELEM	.eq.	629	)		coeff = 	0.788948037
        if (JELEM	.eq.	630	)		coeff = 	0.810227818
        if (JELEM	.eq.	631	)		coeff = 	0.863295577
        if (JELEM	.eq.	632	)		coeff = 	0.773923098
        if (JELEM	.eq.	633	)		coeff = 	1.055807514
        if (JELEM	.eq.	634	)		coeff = 	0.959105852
        if (JELEM	.eq.	635	)		coeff = 	0.926644407
        if (JELEM	.eq.	636	)		coeff = 	0.78671604
        if (JELEM	.eq.	637	)		coeff = 	0.914592449
        if (JELEM	.eq.	638	)		coeff = 	0.994219993
        if (JELEM	.eq.	639	)		coeff = 	0.886891173
        if (JELEM	.eq.	640	)		coeff = 	0.851908364
        if (JELEM	.eq.	641	)		coeff = 	0.934928498
        if (JELEM	.eq.	642	)		coeff = 	1.094420284
        if (JELEM	.eq.	643	)		coeff = 	0.777593954
        if (JELEM	.eq.	644	)		coeff = 	1.113348421
        if (JELEM	.eq.	645	)		coeff = 	0.956428291
        if (JELEM	.eq.	646	)		coeff = 	1.078537422
        if (JELEM	.eq.	647	)		coeff = 	0.996671361
        if (JELEM	.eq.	648	)		coeff = 	0.854368642
        if (JELEM	.eq.	649	)		coeff = 	1.0349703
        if (JELEM	.eq.	650	)		coeff = 	0.810824725
        if (JELEM	.eq.	651	)		coeff = 	0.897408983
        if (JELEM	.eq.	652	)		coeff = 	0.99721213
        if (JELEM	.eq.	653	)		coeff = 	1.033525061
        if (JELEM	.eq.	654	)		coeff = 	1.032353661
        if (JELEM	.eq.	655	)		coeff = 	0.87562176
        if (JELEM	.eq.	656	)		coeff = 	1.116067296
        if (JELEM	.eq.	657	)		coeff = 	0.907276043
        if (JELEM	.eq.	658	)		coeff = 	0.928763797
        if (JELEM	.eq.	659	)		coeff = 	0.884209128
        if (JELEM	.eq.	660	)		coeff = 	0.685781853
        if (JELEM	.eq.	661	)		coeff = 	0.995649612
        if (JELEM	.eq.	662	)		coeff = 	1.030181252
        if (JELEM	.eq.	663	)		coeff = 	0.815890529
        if (JELEM	.eq.	664	)		coeff = 	0.857739278
        if (JELEM	.eq.	665	)		coeff = 	0.97471894
        if (JELEM	.eq.	666	)		coeff = 	1.091371339
        if (JELEM	.eq.	667	)		coeff = 	0.94570709
        if (JELEM	.eq.	668	)		coeff = 	1.052909733
        if (JELEM	.eq.	669	)		coeff = 	0.955104243
        if (JELEM	.eq.	670	)		coeff = 	1.00349876
        if (JELEM	.eq.	671	)		coeff = 	1.014590766
        if (JELEM	.eq.	672	)		coeff = 	0.892080767
        if (JELEM	.eq.	673	)		coeff = 	0.959603061
        if (JELEM	.eq.	674	)		coeff = 	0.865146139
        if (JELEM	.eq.	675	)		coeff = 	1.047282414
        if (JELEM	.eq.	676	)		coeff = 	0.909774344
        if (JELEM	.eq.	677	)		coeff = 	1.114143726
        if (JELEM	.eq.	678	)		coeff = 	0.860146624
        if (JELEM	.eq.	679	)		coeff = 	0.909614512
        if (JELEM	.eq.	680	)		coeff = 	0.74891382
        if (JELEM	.eq.	681	)		coeff = 	1.091073432
        if (JELEM	.eq.	682	)		coeff = 	0.791115498
        if (JELEM	.eq.	683	)		coeff = 	0.961235218
        if (JELEM	.eq.	684	)		coeff = 	0.930442455
        if (JELEM	.eq.	685	)		coeff = 	1.014231854
        if (JELEM	.eq.	686	)		coeff = 	1.09845918
        if (JELEM	.eq.	687	)		coeff = 	0.833483452
        if (JELEM	.eq.	688	)		coeff = 	1.068198956
        if (JELEM	.eq.	689	)		coeff = 	0.999291765
        if (JELEM	.eq.	690	)		coeff = 	0.928100228
        if (JELEM	.eq.	691	)		coeff = 	0.548196283
        if (JELEM	.eq.	692	)		coeff = 	0.959161279
        if (JELEM	.eq.	693	)		coeff = 	0.627954137
        if (JELEM	.eq.	694	)		coeff = 	1.040323713
        if (JELEM	.eq.	695	)		coeff = 	0.991835497
        if (JELEM	.eq.	696	)		coeff = 	0.903278151
        if (JELEM	.eq.	697	)		coeff = 	1.106006751
        if (JELEM	.eq.	698	)		coeff = 	0.883812546
        if (JELEM	.eq.	699	)		coeff = 	0.985778276
        if (JELEM	.eq.	700	)		coeff = 	0.855134746
        if (JELEM	.eq.	701	)		coeff = 	0.996817125
        if (JELEM	.eq.	702	)		coeff = 	1.013397852
        if (JELEM	.eq.	703	)		coeff = 	0.652767874
        if (JELEM	.eq.	704	)		coeff = 	0.895338443
        if (JELEM	.eq.	705	)		coeff = 	0.987730662
        if (JELEM	.eq.	706	)		coeff = 	1.087629383
        if (JELEM	.eq.	707	)		coeff = 	0.889044246
        if (JELEM	.eq.	708	)		coeff = 	0.923357652
        if (JELEM	.eq.	709	)		coeff = 	1.100511527
        if (JELEM	.eq.	710	)		coeff = 	1.077835147
        if (JELEM	.eq.	711	)		coeff = 	0.671264162
        if (JELEM	.eq.	712	)		coeff = 	0.964901653
        if (JELEM	.eq.	713	)		coeff = 	1.142763781
        if (JELEM	.eq.	714	)		coeff = 	1.113182753
        if (JELEM	.eq.	715	)		coeff = 	1.069604419
        if (JELEM	.eq.	716	)		coeff = 	0.803837709
        if (JELEM	.eq.	717	)		coeff = 	0.973434703
        if (JELEM	.eq.	718	)		coeff = 	0.946704739
        if (JELEM	.eq.	719	)		coeff = 	0.965551732
        if (JELEM	.eq.	720	)		coeff = 	1.104073383
        if (JELEM	.eq.	721	)		coeff = 	0.800449407
        if (JELEM	.eq.	722	)		coeff = 	1.021977436
        if (JELEM	.eq.	723	)		coeff = 	1.02759544
        if (JELEM	.eq.	724	)		coeff = 	0.936798467
        if (JELEM	.eq.	725	)		coeff = 	0.970677591
        if (JELEM	.eq.	726	)		coeff = 	0.999882687
        if (JELEM	.eq.	727	)		coeff = 	0.917382685
        if (JELEM	.eq.	728	)		coeff = 	0.759410186
        if (JELEM	.eq.	729	)		coeff = 	0.928717685
        if (JELEM	.eq.	730	)		coeff = 	1.02361346
        if (JELEM	.eq.	731	)		coeff = 	1.047060367
        if (JELEM	.eq.	732	)		coeff = 	0.980702699
        if (JELEM	.eq.	733	)		coeff = 	1.136744166
        if (JELEM	.eq.	734	)		coeff = 	0.816800837
        if (JELEM	.eq.	735	)		coeff = 	0.93193463
        if (JELEM	.eq.	736	)		coeff = 	1.047573152
        if (JELEM	.eq.	737	)		coeff = 	0.958436874
        if (JELEM	.eq.	738	)		coeff = 	1.113187882
        if (JELEM	.eq.	739	)		coeff = 	0.82624533
        if (JELEM	.eq.	740	)		coeff = 	0.97967089
        if (JELEM	.eq.	741	)		coeff = 	0.950431698
        if (JELEM	.eq.	742	)		coeff = 	0.944946912
        if (JELEM	.eq.	743	)		coeff = 	0.908982029
        if (JELEM	.eq.	744	)		coeff = 	0.999049274
        if (JELEM	.eq.	745	)		coeff = 	1.098049613
        if (JELEM	.eq.	746	)		coeff = 	0.976020506
        if (JELEM	.eq.	747	)		coeff = 	1.117595636
        if (JELEM	.eq.	748	)		coeff = 	0.887586352
        if (JELEM	.eq.	749	)		coeff = 	1.125787662
        if (JELEM	.eq.	750	)		coeff = 	0.736304659
        if (JELEM	.eq.	751	)		coeff = 	0.885946887
        if (JELEM	.eq.	752	)		coeff = 	0.760357943
        if (JELEM	.eq.	753	)		coeff = 	0.960279407
        if (JELEM	.eq.	754	)		coeff = 	1.035932065
        if (JELEM	.eq.	755	)		coeff = 	1.030251937
        if (JELEM	.eq.	756	)		coeff = 	0.879138835
        if (JELEM	.eq.	757	)		coeff = 	0.606002421
        if (JELEM	.eq.	758	)		coeff = 	1.003042414
        if (JELEM	.eq.	759	)		coeff = 	0.881684701
        if (JELEM	.eq.	760	)		coeff = 	1.082373047
        if (JELEM	.eq.	761	)		coeff = 	0.935933448
        if (JELEM	.eq.	762	)		coeff = 	0.983017993
        if (JELEM	.eq.	763	)		coeff = 	0.890543788
        if (JELEM	.eq.	764	)		coeff = 	1.029896006
        if (JELEM	.eq.	765	)		coeff = 	1.089471251
        if (JELEM	.eq.	766	)		coeff = 	0.977487988
        if (JELEM	.eq.	767	)		coeff = 	0.922465438
        if (JELEM	.eq.	768	)		coeff = 	1.073087648
        if (JELEM	.eq.	769	)		coeff = 	0.976981884
        if (JELEM	.eq.	770	)		coeff = 	0.918510071
        if (JELEM	.eq.	771	)		coeff = 	0.846963168
        if (JELEM	.eq.	772	)		coeff = 	1.016465442
        if (JELEM	.eq.	773	)		coeff = 	0.99065671
        if (JELEM	.eq.	774	)		coeff = 	0.810967454
        if (JELEM	.eq.	775	)		coeff = 	0.901833537
        if (JELEM	.eq.	776	)		coeff = 	1.035639984
        if (JELEM	.eq.	777	)		coeff = 	0.878778153
        if (JELEM	.eq.	778	)		coeff = 	1.02130556
        if (JELEM	.eq.	779	)		coeff = 	1.025165683
        if (JELEM	.eq.	780	)		coeff = 	1.176905482
        if (JELEM	.eq.	781	)		coeff = 	0.998144569
        if (JELEM	.eq.	782	)		coeff = 	0.981299936
        if (JELEM	.eq.	783	)		coeff = 	1.017525811
        if (JELEM	.eq.	784	)		coeff = 	1.021324157
        if (JELEM	.eq.	785	)		coeff = 	1.035445798
        if (JELEM	.eq.	786	)		coeff = 	0.76131472
        if (JELEM	.eq.	787	)		coeff = 	0.827528555
        if (JELEM	.eq.	788	)		coeff = 	0.992045654
        if (JELEM	.eq.	789	)		coeff = 	0.969704419
        if (JELEM	.eq.	790	)		coeff = 	0.945494002
        if (JELEM	.eq.	791	)		coeff = 	0.966891054
        if (JELEM	.eq.	792	)		coeff = 	1.027373039
        if (JELEM	.eq.	793	)		coeff = 	0.632904062
        if (JELEM	.eq.	794	)		coeff = 	1.054159434
        if (JELEM	.eq.	795	)		coeff = 	0.826587778
        if (JELEM	.eq.	796	)		coeff = 	1.13092457
        if (JELEM	.eq.	797	)		coeff = 	1.009828397
        if (JELEM	.eq.	798	)		coeff = 	0.883370659
        if (JELEM	.eq.	799	)		coeff = 	0.921074312
        if (JELEM	.eq.	800	)		coeff = 	1.059144673
        if (JELEM	.eq.	801	)		coeff = 	0.739600659
        if (JELEM	.eq.	802	)		coeff = 	0.951795284
        if (JELEM	.eq.	803	)		coeff = 	1.032781763
        if (JELEM	.eq.	804	)		coeff = 	0.941486189
        if (JELEM	.eq.	805	)		coeff = 	0.784514947
        if (JELEM	.eq.	806	)		coeff = 	0.802128368
        if (JELEM	.eq.	807	)		coeff = 	0.968854345
        if (JELEM	.eq.	808	)		coeff = 	0.979723767
        if (JELEM	.eq.	809	)		coeff = 	1.015436188
        if (JELEM	.eq.	810	)		coeff = 	1.112143573
        if (JELEM	.eq.	811	)		coeff = 	0.881298896
        if (JELEM	.eq.	812	)		coeff = 	1.073124318
        if (JELEM	.eq.	813	)		coeff = 	1.003256204
        if (JELEM	.eq.	814	)		coeff = 	0.992412826
        if (JELEM	.eq.	815	)		coeff = 	0.810001462
        if (JELEM	.eq.	816	)		coeff = 	1.144372443
        if (JELEM	.eq.	817	)		coeff = 	0.837372708
        if (JELEM	.eq.	818	)		coeff = 	1.022118429
        if (JELEM	.eq.	819	)		coeff = 	1.033100536
        if (JELEM	.eq.	820	)		coeff = 	0.967233772
        if (JELEM	.eq.	821	)		coeff = 	0.89119575
        if (JELEM	.eq.	822	)		coeff = 	1.047889096
        if (JELEM	.eq.	823	)		coeff = 	1.043525732
        if (JELEM	.eq.	824	)		coeff = 	0.688542213
        if (JELEM	.eq.	825	)		coeff = 	0.937088562
        if (JELEM	.eq.	826	)		coeff = 	1.017496271
        if (JELEM	.eq.	827	)		coeff = 	0.710603469
        if (JELEM	.eq.	828	)		coeff = 	0.801824695
        if (JELEM	.eq.	829	)		coeff = 	1.052025891
        if (JELEM	.eq.	830	)		coeff = 	1.202469512
        if (JELEM	.eq.	831	)		coeff = 	0.897736884
        if (JELEM	.eq.	832	)		coeff = 	0.822605521
        if (JELEM	.eq.	833	)		coeff = 	1.07877009
        if (JELEM	.eq.	834	)		coeff = 	1.124871184
        if (JELEM	.eq.	835	)		coeff = 	0.93558388
        if (JELEM	.eq.	836	)		coeff = 	0.933714134
        if (JELEM	.eq.	837	)		coeff = 	0.95941768
        if (JELEM	.eq.	838	)		coeff = 	1.17205354
        if (JELEM	.eq.	839	)		coeff = 	0.905994446
        if (JELEM	.eq.	840	)		coeff = 	0.748897515
        if (JELEM	.eq.	841	)		coeff = 	0.818652153
        if (JELEM	.eq.	842	)		coeff = 	1.080946679
        if (JELEM	.eq.	843	)		coeff = 	1.003625413
        if (JELEM	.eq.	844	)		coeff = 	1.035620517
        if (JELEM	.eq.	845	)		coeff = 	0.946857483
        if (JELEM	.eq.	846	)		coeff = 	0.931118673
        if (JELEM	.eq.	847	)		coeff = 	1.018507908
        if (JELEM	.eq.	848	)		coeff = 	0.861605981
        if (JELEM	.eq.	849	)		coeff = 	0.862797426
        if (JELEM	.eq.	850	)		coeff = 	0.869522041
        if (JELEM	.eq.	851	)		coeff = 	1.004573579
        if (JELEM	.eq.	852	)		coeff = 	1.112854073
        if (JELEM	.eq.	853	)		coeff = 	0.898867707
        if (JELEM	.eq.	854	)		coeff = 	0.609685023
        if (JELEM	.eq.	855	)		coeff = 	1.06155084
        if (JELEM	.eq.	856	)		coeff = 	1.08081976
        if (JELEM	.eq.	857	)		coeff = 	0.787006747
        if (JELEM	.eq.	858	)		coeff = 	0.969088398
        if (JELEM	.eq.	859	)		coeff = 	0.832759312
        if (JELEM	.eq.	860	)		coeff = 	0.855840267
        if (JELEM	.eq.	861	)		coeff = 	1.053119054
        if (JELEM	.eq.	862	)		coeff = 	1.034039165
        if (JELEM	.eq.	863	)		coeff = 	1.112930466
        if (JELEM	.eq.	864	)		coeff = 	0.932283921
        if (JELEM	.eq.	865	)		coeff = 	0.871249984
        if (JELEM	.eq.	866	)		coeff = 	0.96092412
        if (JELEM	.eq.	867	)		coeff = 	1.136169314
        if (JELEM	.eq.	868	)		coeff = 	0.628786602
        if (JELEM	.eq.	869	)		coeff = 	0.963750717
        if (JELEM	.eq.	870	)		coeff = 	1.00981751
        if (JELEM	.eq.	871	)		coeff = 	1.057519082
        if (JELEM	.eq.	872	)		coeff = 	0.92707737
        if (JELEM	.eq.	873	)		coeff = 	0.942506497
        if (JELEM	.eq.	874	)		coeff = 	0.882455093
        if (JELEM	.eq.	875	)		coeff = 	1.064818011
        if (JELEM	.eq.	876	)		coeff = 	1.003020197
        if (JELEM	.eq.	877	)		coeff = 	1.068419324
        if (JELEM	.eq.	878	)		coeff = 	0.833487335
        if (JELEM	.eq.	879	)		coeff = 	1.008200178
        if (JELEM	.eq.	880	)		coeff = 	1.025812666
        if (JELEM	.eq.	881	)		coeff = 	1.177273869
        if (JELEM	.eq.	882	)		coeff = 	0.859787919
        if (JELEM	.eq.	883	)		coeff = 	0.964732582
        if (JELEM	.eq.	884	)		coeff = 	0.953350058
        if (JELEM	.eq.	885	)		coeff = 	0.820455671
        if (JELEM	.eq.	886	)		coeff = 	0.893572565
        if (JELEM	.eq.	887	)		coeff = 	0.913202234
        if (JELEM	.eq.	888	)		coeff = 	1.055811685
        if (JELEM	.eq.	889	)		coeff = 	0.949713983
        if (JELEM	.eq.	890	)		coeff = 	0.726486668
        if (JELEM	.eq.	891	)		coeff = 	0.936244351
        if (JELEM	.eq.	892	)		coeff = 	0.856511634
        if (JELEM	.eq.	893	)		coeff = 	0.659599747
        if (JELEM	.eq.	894	)		coeff = 	0.80969765
        if (JELEM	.eq.	895	)		coeff = 	1.044296152
        if (JELEM	.eq.	896	)		coeff = 	1.128945859
        if (JELEM	.eq.	897	)		coeff = 	0.977444356
        if (JELEM	.eq.	898	)		coeff = 	1.156590675
        if (JELEM	.eq.	899	)		coeff = 	0.97127289
        if (JELEM	.eq.	900	)		coeff = 	0.741184819
        if (JELEM	.eq.	901	)		coeff = 	1.033515308
        if (JELEM	.eq.	902	)		coeff = 	0.994969976
        if (JELEM	.eq.	903	)		coeff = 	0.98279088
        if (JELEM	.eq.	904	)		coeff = 	0.844848149
        if (JELEM	.eq.	905	)		coeff = 	0.848224322
        if (JELEM	.eq.	906	)		coeff = 	0.976933903
        if (JELEM	.eq.	907	)		coeff = 	0.99655828
        if (JELEM	.eq.	908	)		coeff = 	0.773859383
        if (JELEM	.eq.	909	)		coeff = 	0.886558015
        if (JELEM	.eq.	910	)		coeff = 	0.887833334
        if (JELEM	.eq.	911	)		coeff = 	0.747610528
        if (JELEM	.eq.	912	)		coeff = 	0.961194164
        if (JELEM	.eq.	913	)		coeff = 	0.864569914
        if (JELEM	.eq.	914	)		coeff = 	0.977158589
        if (JELEM	.eq.	915	)		coeff = 	0.834337377
        if (JELEM	.eq.	916	)		coeff = 	0.993884885
        if (JELEM	.eq.	917	)		coeff = 	0.887538993
        if (JELEM	.eq.	918	)		coeff = 	0.688235251
        if (JELEM	.eq.	919	)		coeff = 	0.957482465
        if (JELEM	.eq.	920	)		coeff = 	0.983202197
        if (JELEM	.eq.	921	)		coeff = 	1.04642476
        if (JELEM	.eq.	922	)		coeff = 	1.012194678
        if (JELEM	.eq.	923	)		coeff = 	1.082007628
        if (JELEM	.eq.	924	)		coeff = 	0.998009186
        if (JELEM	.eq.	925	)		coeff = 	1.01039535
        if (JELEM	.eq.	926	)		coeff = 	1.007036658
        if (JELEM	.eq.	927	)		coeff = 	0.852175103
        if (JELEM	.eq.	928	)		coeff = 	0.955105458
        if (JELEM	.eq.	929	)		coeff = 	0.958097052
        if (JELEM	.eq.	930	)		coeff = 	0.872543778
        if (JELEM	.eq.	931	)		coeff = 	1.077938932
        if (JELEM	.eq.	932	)		coeff = 	0.927136207
        if (JELEM	.eq.	933	)		coeff = 	1.005791076
        if (JELEM	.eq.	934	)		coeff = 	1.009094766
        if (JELEM	.eq.	935	)		coeff = 	0.942659995
        if (JELEM	.eq.	936	)		coeff = 	0.825108345
        if (JELEM	.eq.	937	)		coeff = 	1.049207919
        if (JELEM	.eq.	938	)		coeff = 	0.911716934
        if (JELEM	.eq.	939	)		coeff = 	0.796932774
        if (JELEM	.eq.	940	)		coeff = 	1.049013864
        if (JELEM	.eq.	941	)		coeff = 	1.019220338
        if (JELEM	.eq.	942	)		coeff = 	0.964979148
        if (JELEM	.eq.	943	)		coeff = 	0.806650268
        if (JELEM	.eq.	944	)		coeff = 	0.963618593
        if (JELEM	.eq.	945	)		coeff = 	1.025293522
        if (JELEM	.eq.	946	)		coeff = 	0.954458533
        if (JELEM	.eq.	947	)		coeff = 	0.942768642
        if (JELEM	.eq.	948	)		coeff = 	0.987833155
        if (JELEM	.eq.	949	)		coeff = 	1.154744225
        if (JELEM	.eq.	950	)		coeff = 	0.901041118
        if (JELEM	.eq.	951	)		coeff = 	0.962131607
        if (JELEM	.eq.	952	)		coeff = 	0.996366949
        if (JELEM	.eq.	953	)		coeff = 	1.105825377
        if (JELEM	.eq.	954	)		coeff = 	1.00253133
        if (JELEM	.eq.	955	)		coeff = 	1.037952661
        if (JELEM	.eq.	956	)		coeff = 	1.047601709
        if (JELEM	.eq.	957	)		coeff = 	0.853887346
        if (JELEM	.eq.	958	)		coeff = 	0.993071999
        if (JELEM	.eq.	959	)		coeff = 	1.113364498
        if (JELEM	.eq.	960	)		coeff = 	0.998052188
        if (JELEM	.eq.	961	)		coeff = 	0.872241896
        if (JELEM	.eq.	962	)		coeff = 	1.060540933
        if (JELEM	.eq.	963	)		coeff = 	0.787625696
        if (JELEM	.eq.	964	)		coeff = 	1.013359495
        if (JELEM	.eq.	965	)		coeff = 	1.010429684
        if (JELEM	.eq.	966	)		coeff = 	1.047365332
        if (JELEM	.eq.	967	)		coeff = 	0.87563419
        if (JELEM	.eq.	968	)		coeff = 	1.102763406
        if (JELEM	.eq.	969	)		coeff = 	0.742977724
        if (JELEM	.eq.	970	)		coeff = 	1.063095224
        if (JELEM	.eq.	971	)		coeff = 	1.022582009
        if (JELEM	.eq.	972	)		coeff = 	0.906630051
        if (JELEM	.eq.	973	)		coeff = 	1.069497542
        if (JELEM	.eq.	974	)		coeff = 	0.960637597
        if (JELEM	.eq.	975	)		coeff = 	0.894149152
        if (JELEM	.eq.	976	)		coeff = 	0.770572374
        if (JELEM	.eq.	977	)		coeff = 	0.890004482
        if (JELEM	.eq.	978	)		coeff = 	0.882921034
        if (JELEM	.eq.	979	)		coeff = 	0.989316455
        if (JELEM	.eq.	980	)		coeff = 	1.036359448
        if (JELEM	.eq.	981	)		coeff = 	0.95817365
        if (JELEM	.eq.	982	)		coeff = 	1.042645391
        if (JELEM	.eq.	983	)		coeff = 	0.838359189
        if (JELEM	.eq.	984	)		coeff = 	0.914940373
        if (JELEM	.eq.	985	)		coeff = 	0.852630494
        if (JELEM	.eq.	986	)		coeff = 	0.863631491
        if (JELEM	.eq.	987	)		coeff = 	0.972535515
        if (JELEM	.eq.	988	)		coeff = 	1.016059879
        if (JELEM	.eq.	989	)		coeff = 	0.906471697
        if (JELEM	.eq.	990	)		coeff = 	0.648891834
        if (JELEM	.eq.	991	)		coeff = 	0.874458033
        if (JELEM	.eq.	992	)		coeff = 	0.845562582
        if (JELEM	.eq.	993	)		coeff = 	0.899829896
        if (JELEM	.eq.	994	)		coeff = 	0.936452273
        if (JELEM	.eq.	995	)		coeff = 	0.881677985
        if (JELEM	.eq.	996	)		coeff = 	0.96492309
        if (JELEM	.eq.	997	)		coeff = 	0.824325194
        if (JELEM	.eq.	998	)		coeff = 	1.10392698
        if (JELEM	.eq.	999	)		coeff = 	0.708695722
        if (JELEM	.eq.	1000	)		coeff = 	1.087569387
        if (JELEM	.eq.	1001	)		coeff = 	0.950729159
        if (JELEM	.eq.	1002	)		coeff = 	0.990492735
        if (JELEM	.eq.	1003	)		coeff = 	1.083724849
        if (JELEM	.eq.	1004	)		coeff = 	0.893067994
        if (JELEM	.eq.	1005	)		coeff = 	0.930822132
        if (JELEM	.eq.	1006	)		coeff = 	0.868674087
        if (JELEM	.eq.	1007	)		coeff = 	0.945017457
        if (JELEM	.eq.	1008	)		coeff = 	0.85513387
        if (JELEM	.eq.	1009	)		coeff = 	0.942019891
        if (JELEM	.eq.	1010	)		coeff = 	0.751693651
        if (JELEM	.eq.	1011	)		coeff = 	0.820087234
        if (JELEM	.eq.	1012	)		coeff = 	0.961892709
        if (JELEM	.eq.	1013	)		coeff = 	0.865991977
        if (JELEM	.eq.	1014	)		coeff = 	0.971459498
        if (JELEM	.eq.	1015	)		coeff = 	0.845986882
        if (JELEM	.eq.	1016	)		coeff = 	1.01244977
        if (JELEM	.eq.	1017	)		coeff = 	0.689093798
        if (JELEM	.eq.	1018	)		coeff = 	1.024943196
        if (JELEM	.eq.	1019	)		coeff = 	1.101092407
        if (JELEM	.eq.	1020	)		coeff = 	0.88235731
        if (JELEM	.eq.	1021	)		coeff = 	0.844682207
        if (JELEM	.eq.	1022	)		coeff = 	0.777614717
        if (JELEM	.eq.	1023	)		coeff = 	1.011191322
        if (JELEM	.eq.	1024	)		coeff = 	0.858735794
        if (JELEM	.eq.	1025	)		coeff = 	0.953233574
        if (JELEM	.eq.	1026	)		coeff = 	0.974125115
        if (JELEM	.eq.	1027	)		coeff = 	0.850285391
        if (JELEM	.eq.	1028	)		coeff = 	0.740072935
        if (JELEM	.eq.	1029	)		coeff = 	1.099135832
        if (JELEM	.eq.	1030	)		coeff = 	0.898880066
        if (JELEM	.eq.	1031	)		coeff = 	1.037746637
        if (JELEM	.eq.	1032	)		coeff = 	0.991594664
        if (JELEM	.eq.	1033	)		coeff = 	1.027863926
        if (JELEM	.eq.	1034	)		coeff = 	0.843951875
        if (JELEM	.eq.	1035	)		coeff = 	0.584282655
        if (JELEM	.eq.	1036	)		coeff = 	0.919315241
        if (JELEM	.eq.	1037	)		coeff = 	0.900608717
        if (JELEM	.eq.	1038	)		coeff = 	0.766602207
        if (JELEM	.eq.	1039	)		coeff = 	0.906445896
        if (JELEM	.eq.	1040	)		coeff = 	0.944508101
        if (JELEM	.eq.	1041	)		coeff = 	0.996480393
        if (JELEM	.eq.	1042	)		coeff = 	0.924232745
        if (JELEM	.eq.	1043	)		coeff = 	1.001265057
        if (JELEM	.eq.	1044	)		coeff = 	0.98923537
        if (JELEM	.eq.	1045	)		coeff = 	0.999777049
        if (JELEM	.eq.	1046	)		coeff = 	0.972728977
        if (JELEM	.eq.	1047	)		coeff = 	0.963061678
        if (JELEM	.eq.	1048	)		coeff = 	0.789194759
        if (JELEM	.eq.	1049	)		coeff = 	1.046658226
        if (JELEM	.eq.	1050	)		coeff = 	1.007998385
        if (JELEM	.eq.	1051	)		coeff = 	0.94280532
        if (JELEM	.eq.	1052	)		coeff = 	0.967633554
        if (JELEM	.eq.	1053	)		coeff = 	1.029589917
        if (JELEM	.eq.	1054	)		coeff = 	0.941182342
        if (JELEM	.eq.	1055	)		coeff = 	0.815278314
        if (JELEM	.eq.	1056	)		coeff = 	1.10835805
        if (JELEM	.eq.	1057	)		coeff = 	0.980229773
        if (JELEM	.eq.	1058	)		coeff = 	1.094807384
        if (JELEM	.eq.	1059	)		coeff = 	0.946008963
        if (JELEM	.eq.	1060	)		coeff = 	0.952930309
        if (JELEM	.eq.	1061	)		coeff = 	0.875270556
        if (JELEM	.eq.	1062	)		coeff = 	1.038307435
        if (JELEM	.eq.	1063	)		coeff = 	0.938859274
        if (JELEM	.eq.	1064	)		coeff = 	0.975305056
        if (JELEM	.eq.	1065	)		coeff = 	0.827021931
        if (JELEM	.eq.	1066	)		coeff = 	0.915643027
        if (JELEM	.eq.	1067	)		coeff = 	1.003813852
        if (JELEM	.eq.	1068	)		coeff = 	1.005643858
        if (JELEM	.eq.	1069	)		coeff = 	1.032097147
        if (JELEM	.eq.	1070	)		coeff = 	0.739086805
        if (JELEM	.eq.	1071	)		coeff = 	1.019243828
        if (JELEM	.eq.	1072	)		coeff = 	1.063023711
        if (JELEM	.eq.	1073	)		coeff = 	1.001790975
        if (JELEM	.eq.	1074	)		coeff = 	0.886244328
        if (JELEM	.eq.	1075	)		coeff = 	0.899897537
        if (JELEM	.eq.	1076	)		coeff = 	0.901720393
        if (JELEM	.eq.	1077	)		coeff = 	1.176463058
        if (JELEM	.eq.	1078	)		coeff = 	0.998243769
        if (JELEM	.eq.	1079	)		coeff = 	0.797219505
        if (JELEM	.eq.	1080	)		coeff = 	1.013598916
        if (JELEM	.eq.	1081	)		coeff = 	0.935921992
        if (JELEM	.eq.	1082	)		coeff = 	1.019358481
        if (JELEM	.eq.	1083	)		coeff = 	1.075950443
        if (JELEM	.eq.	1084	)		coeff = 	0.994446624
        if (JELEM	.eq.	1085	)		coeff = 	0.851880994
        if (JELEM	.eq.	1086	)		coeff = 	0.672796744
        if (JELEM	.eq.	1087	)		coeff = 	0.826372957
        if (JELEM	.eq.	1088	)		coeff = 	1.095037352
        if (JELEM	.eq.	1089	)		coeff = 	1.008243561
        if (JELEM	.eq.	1090	)		coeff = 	1.037380206
        if (JELEM	.eq.	1091	)		coeff = 	1.013748998
        if (JELEM	.eq.	1092	)		coeff = 	0.659959675
        if (JELEM	.eq.	1093	)		coeff = 	0.950361985
        if (JELEM	.eq.	1094	)		coeff = 	0.883164907
        if (JELEM	.eq.	1095	)		coeff = 	0.838669883
        if (JELEM	.eq.	1096	)		coeff = 	1.059974808
        if (JELEM	.eq.	1097	)		coeff = 	0.795843947
        if (JELEM	.eq.	1098	)		coeff = 	1.084597671
        if (JELEM	.eq.	1099	)		coeff = 	0.884852888
        if (JELEM	.eq.	1100	)		coeff = 	0.891065797
        if (JELEM	.eq.	1101	)		coeff = 	0.895606092
        if (JELEM	.eq.	1102	)		coeff = 	1.072530139
        if (JELEM	.eq.	1103	)		coeff = 	0.978895368
        if (JELEM	.eq.	1104	)		coeff = 	0.961561402
        if (JELEM	.eq.	1105	)		coeff = 	0.955445125
        if (JELEM	.eq.	1106	)		coeff = 	0.827830416
        if (JELEM	.eq.	1107	)		coeff = 	0.909905944
        if (JELEM	.eq.	1108	)		coeff = 	0.857866438
        if (JELEM	.eq.	1109	)		coeff = 	0.955231254
        if (JELEM	.eq.	1110	)		coeff = 	0.733539908
        if (JELEM	.eq.	1111	)		coeff = 	1.104747936
        if (JELEM	.eq.	1112	)		coeff = 	0.952298027
        if (JELEM	.eq.	1113	)		coeff = 	1.023951747
        if (JELEM	.eq.	1114	)		coeff = 	0.969296753
        if (JELEM	.eq.	1115	)		coeff = 	0.907424793
        if (JELEM	.eq.	1116	)		coeff = 	1.046072107
        if (JELEM	.eq.	1117	)		coeff = 	0.932530859
        if (JELEM	.eq.	1118	)		coeff = 	1.011426009
        if (JELEM	.eq.	1119	)		coeff = 	0.813460131
        if (JELEM	.eq.	1120	)		coeff = 	1.072544968
        if (JELEM	.eq.	1121	)		coeff = 	1.085851606
        if (JELEM	.eq.	1122	)		coeff = 	0.727876417
        if (JELEM	.eq.	1123	)		coeff = 	1.065046549
        if (JELEM	.eq.	1124	)		coeff = 	1.065181614
        if (JELEM	.eq.	1125	)		coeff = 	1.064067804
        if (JELEM	.eq.	1126	)		coeff = 	1.092071515
        if (JELEM	.eq.	1127	)		coeff = 	0.976545471
        if (JELEM	.eq.	1128	)		coeff = 	0.912927604
        if (JELEM	.eq.	1129	)		coeff = 	0.84462005
        if (JELEM	.eq.	1130	)		coeff = 	0.865336872
        if (JELEM	.eq.	1131	)		coeff = 	0.897385145
        if (JELEM	.eq.	1132	)		coeff = 	0.971578313
        if (JELEM	.eq.	1133	)		coeff = 	0.898910603
        if (JELEM	.eq.	1134	)		coeff = 	0.729777442
        if (JELEM	.eq.	1135	)		coeff = 	0.962402905
        if (JELEM	.eq.	1136	)		coeff = 	1.017313778
        if (JELEM	.eq.	1137	)		coeff = 	0.865533917
        if (JELEM	.eq.	1138	)		coeff = 	1.037301756
        if (JELEM	.eq.	1139	)		coeff = 	1.037936911
        if (JELEM	.eq.	1140	)		coeff = 	0.973741187
        if (JELEM	.eq.	1141	)		coeff = 	0.929043068
        if (JELEM	.eq.	1142	)		coeff = 	0.930309912
        if (JELEM	.eq.	1143	)		coeff = 	1.076945328
        if (JELEM	.eq.	1144	)		coeff = 	1.076457179
        if (JELEM	.eq.	1145	)		coeff = 	1.023149906
        if (JELEM	.eq.	1146	)		coeff = 	0.888583056
        if (JELEM	.eq.	1147	)		coeff = 	0.988237833
        if (JELEM	.eq.	1148	)		coeff = 	0.845888888
        if (JELEM	.eq.	1149	)		coeff = 	0.763244438
        if (JELEM	.eq.	1150	)		coeff = 	0.991548386
        if (JELEM	.eq.	1151	)		coeff = 	1.114341554
        if (JELEM	.eq.	1152	)		coeff = 	0.943673849
        if (JELEM	.eq.	1153	)		coeff = 	0.883809134
        if (JELEM	.eq.	1154	)		coeff = 	1.013071084
        if (JELEM	.eq.	1155	)		coeff = 	0.965973555
        if (JELEM	.eq.	1156	)		coeff = 	1.04184004
        if (JELEM	.eq.	1157	)		coeff = 	0.758037936
        if (JELEM	.eq.	1158	)		coeff = 	0.968909814
        if (JELEM	.eq.	1159	)		coeff = 	0.952734086
        if (JELEM	.eq.	1160	)		coeff = 	1.042025969
        if (JELEM	.eq.	1161	)		coeff = 	1.088918542
        if (JELEM	.eq.	1162	)		coeff = 	1.108870408
        if (JELEM	.eq.	1163	)		coeff = 	0.850934691
        if (JELEM	.eq.	1164	)		coeff = 	0.873785472
        if (JELEM	.eq.	1165	)		coeff = 	1.05014814
        if (JELEM	.eq.	1166	)		coeff = 	0.802516187
        if (JELEM	.eq.	1167	)		coeff = 	0.907621968
        if (JELEM	.eq.	1168	)		coeff = 	0.916972795
        if (JELEM	.eq.	1169	)		coeff = 	0.628813409
        if (JELEM	.eq.	1170	)		coeff = 	1.129864563
        if (JELEM	.eq.	1171	)		coeff = 	0.984735506
        if (JELEM	.eq.	1172	)		coeff = 	0.966790309
        if (JELEM	.eq.	1173	)		coeff = 	0.940015366
        if (JELEM	.eq.	1174	)		coeff = 	1.095328396
        if (JELEM	.eq.	1175	)		coeff = 	0.915870826
        if (JELEM	.eq.	1176	)		coeff = 	1.11427376
        if (JELEM	.eq.	1177	)		coeff = 	0.94787762
        if (JELEM	.eq.	1178	)		coeff = 	0.897639855
        if (JELEM	.eq.	1179	)		coeff = 	0.967363367
        if (JELEM	.eq.	1180	)		coeff = 	0.929601649
        if (JELEM	.eq.	1181	)		coeff = 	1.044314968
        if (JELEM	.eq.	1182	)		coeff = 	0.920653582
        if (JELEM	.eq.	1183	)		coeff = 	0.996536021
        if (JELEM	.eq.	1184	)		coeff = 	1.085246689
        if (JELEM	.eq.	1185	)		coeff = 	0.997384858
        if (JELEM	.eq.	1186	)		coeff = 	1.029401726
        if (JELEM	.eq.	1187	)		coeff = 	1.035817048
        if (JELEM	.eq.	1188	)		coeff = 	0.92793267
        if (JELEM	.eq.	1189	)		coeff = 	0.957587355
        if (JELEM	.eq.	1190	)		coeff = 	0.987715952
        if (JELEM	.eq.	1191	)		coeff = 	1.043048764
        if (JELEM	.eq.	1192	)		coeff = 	0.828569103
        if (JELEM	.eq.	1193	)		coeff = 	0.827004775
        if (JELEM	.eq.	1194	)		coeff = 	1.023298714
        if (JELEM	.eq.	1195	)		coeff = 	0.930290157
        if (JELEM	.eq.	1196	)		coeff = 	0.870242266
        if (JELEM	.eq.	1197	)		coeff = 	0.735305479
        if (JELEM	.eq.	1198	)		coeff = 	0.780424897
        if (JELEM	.eq.	1199	)		coeff = 	0.995405226
        if (JELEM	.eq.	1200	)		coeff = 	1.061490982
        if (JELEM	.eq.	1201	)		coeff = 	0.862269909
        if (JELEM	.eq.	1202	)		coeff = 	1.08071009
        if (JELEM	.eq.	1203	)		coeff = 	1.062871497
        if (JELEM	.eq.	1204	)		coeff = 	1.003281854
        if (JELEM	.eq.	1205	)		coeff = 	0.835199698
        if (JELEM	.eq.	1206	)		coeff = 	0.94023157
        if (JELEM	.eq.	1207	)		coeff = 	0.939214715
        if (JELEM	.eq.	1208	)		coeff = 	0.773892532
        if (JELEM	.eq.	1209	)		coeff = 	0.942523254
        if (JELEM	.eq.	1210	)		coeff = 	1.165053597
        if (JELEM	.eq.	1211	)		coeff = 	0.856108889
        if (JELEM	.eq.	1212	)		coeff = 	0.932331171
        if (JELEM	.eq.	1213	)		coeff = 	0.969568333
        if (JELEM	.eq.	1214	)		coeff = 	1.027771991
        if (JELEM	.eq.	1215	)		coeff = 	1.030802364
        if (JELEM	.eq.	1216	)		coeff = 	0.969266687
        if (JELEM	.eq.	1217	)		coeff = 	1.040074294
        if (JELEM	.eq.	1218	)		coeff = 	1.117015994
        if (JELEM	.eq.	1219	)		coeff = 	1.059143316
        if (JELEM	.eq.	1220	)		coeff = 	1.030699203
        if (JELEM	.eq.	1221	)		coeff = 	1.049421513
        if (JELEM	.eq.	1222	)		coeff = 	0.933290546
        if (JELEM	.eq.	1223	)		coeff = 	0.848726455
        if (JELEM	.eq.	1224	)		coeff = 	0.855488335
        if (JELEM	.eq.	1225	)		coeff = 	0.859629951
        if (JELEM	.eq.	1226	)		coeff = 	0.89908535
        if (JELEM	.eq.	1227	)		coeff = 	0.828046285
        if (JELEM	.eq.	1228	)		coeff = 	0.869525532
        if (JELEM	.eq.	1229	)		coeff = 	1.047507065
        if (JELEM	.eq.	1230	)		coeff = 	0.605997981
        if (JELEM	.eq.	1231	)		coeff = 	1.0900506
        if (JELEM	.eq.	1232	)		coeff = 	0.919030343
        if (JELEM	.eq.	1233	)		coeff = 	1.043872617
        if (JELEM	.eq.	1234	)		coeff = 	1.035037304
        if (JELEM	.eq.	1235	)		coeff = 	1.007700441
        if (JELEM	.eq.	1236	)		coeff = 	1.049437133
        if (JELEM	.eq.	1237	)		coeff = 	0.962107522
        if (JELEM	.eq.	1238	)		coeff = 	0.741866005
        if (JELEM	.eq.	1239	)		coeff = 	0.992754723
        if (JELEM	.eq.	1240	)		coeff = 	0.93972801
        if (JELEM	.eq.	1241	)		coeff = 	0.933030965
        if (JELEM	.eq.	1242	)		coeff = 	0.896665279
        if (JELEM	.eq.	1243	)		coeff = 	0.990880428
        if (JELEM	.eq.	1244	)		coeff = 	0.828482118
        if (JELEM	.eq.	1245	)		coeff = 	0.779527892
        if (JELEM	.eq.	1246	)		coeff = 	0.882515192
        if (JELEM	.eq.	1247	)		coeff = 	1.022826013
        if (JELEM	.eq.	1248	)		coeff = 	0.862239246
        if (JELEM	.eq.	1249	)		coeff = 	1.068874081
        if (JELEM	.eq.	1250	)		coeff = 	0.962744829
        if (JELEM	.eq.	1251	)		coeff = 	0.931750408
        if (JELEM	.eq.	1252	)		coeff = 	0.900665596
        if (JELEM	.eq.	1253	)		coeff = 	0.995805002
        if (JELEM	.eq.	1254	)		coeff = 	0.891332108
        if (JELEM	.eq.	1255	)		coeff = 	0.808666159
        if (JELEM	.eq.	1256	)		coeff = 	1.111773518
        if (JELEM	.eq.	1257	)		coeff = 	1.070639237
        if (JELEM	.eq.	1258	)		coeff = 	0.825676724
        if (JELEM	.eq.	1259	)		coeff = 	0.985414014
        if (JELEM	.eq.	1260	)		coeff = 	0.988236191
        if (JELEM	.eq.	1261	)		coeff = 	0.727831976
        if (JELEM	.eq.	1262	)		coeff = 	0.882758456
        if (JELEM	.eq.	1263	)		coeff = 	0.673438441
        if (JELEM	.eq.	1264	)		coeff = 	1.03817732
        if (JELEM	.eq.	1265	)		coeff = 	1.088782239
        if (JELEM	.eq.	1266	)		coeff = 	0.995468475
        if (JELEM	.eq.	1267	)		coeff = 	0.963936418
        if (JELEM	.eq.	1268	)		coeff = 	0.943948135
        if (JELEM	.eq.	1269	)		coeff = 	0.687662519
        if (JELEM	.eq.	1270	)		coeff = 	0.965995478
        if (JELEM	.eq.	1271	)		coeff = 	0.991057165
        if (JELEM	.eq.	1272	)		coeff = 	0.588966395
        if (JELEM	.eq.	1273	)		coeff = 	1.029948341
        if (JELEM	.eq.	1274	)		coeff = 	0.914135631
        if (JELEM	.eq.	1275	)		coeff = 	0.71797205
        if (JELEM	.eq.	1276	)		coeff = 	0.91216053
        if (JELEM	.eq.	1277	)		coeff = 	1.018968845
        if (JELEM	.eq.	1278	)		coeff = 	0.955264325
        if (JELEM	.eq.	1279	)		coeff = 	1.206347242
        if (JELEM	.eq.	1280	)		coeff = 	0.811296789
        if (JELEM	.eq.	1281	)		coeff = 	0.990110357
        if (JELEM	.eq.	1282	)		coeff = 	1.018396171
        if (JELEM	.eq.	1283	)		coeff = 	0.742138151
        if (JELEM	.eq.	1284	)		coeff = 	0.97484806
        if (JELEM	.eq.	1285	)		coeff = 	1.022239087
        if (JELEM	.eq.	1286	)		coeff = 	1.094614034
        if (JELEM	.eq.	1287	)		coeff = 	0.940410747
        if (JELEM	.eq.	1288	)		coeff = 	1.064984737
        if (JELEM	.eq.	1289	)		coeff = 	1.100940514
        if (JELEM	.eq.	1290	)		coeff = 	0.940893001
        if (JELEM	.eq.	1291	)		coeff = 	1.022419612
        if (JELEM	.eq.	1292	)		coeff = 	1.001621305
        if (JELEM	.eq.	1293	)		coeff = 	0.892816153
        if (JELEM	.eq.	1294	)		coeff = 	0.828722099
        if (JELEM	.eq.	1295	)		coeff = 	1.00544418
        if (JELEM	.eq.	1296	)		coeff = 	0.722956838
        if (JELEM	.eq.	1297	)		coeff = 	0.73741006
        if (JELEM	.eq.	1298	)		coeff = 	1.046784777
        if (JELEM	.eq.	1299	)		coeff = 	0.875189796
        if (JELEM	.eq.	1300	)		coeff = 	0.930246862
        if (JELEM	.eq.	1301	)		coeff = 	0.781111733
        if (JELEM	.eq.	1302	)		coeff = 	0.934258404
        if (JELEM	.eq.	1303	)		coeff = 	0.901268712
        if (JELEM	.eq.	1304	)		coeff = 	0.885424131
        if (JELEM	.eq.	1305	)		coeff = 	0.995327606
        if (JELEM	.eq.	1306	)		coeff = 	1.032754893
        if (JELEM	.eq.	1307	)		coeff = 	1.126926229
        if (JELEM	.eq.	1308	)		coeff = 	0.971708486
        if (JELEM	.eq.	1309	)		coeff = 	0.920844597
        if (JELEM	.eq.	1310	)		coeff = 	1.024724754
        if (JELEM	.eq.	1311	)		coeff = 	0.959015874
        if (JELEM	.eq.	1312	)		coeff = 	1.034498579
        if (JELEM	.eq.	1313	)		coeff = 	1.019439593
        if (JELEM	.eq.	1314	)		coeff = 	0.919077808
        if (JELEM	.eq.	1315	)		coeff = 	0.805486517
        if (JELEM	.eq.	1316	)		coeff = 	0.826941721
        if (JELEM	.eq.	1317	)		coeff = 	1.045542626
        if (JELEM	.eq.	1318	)		coeff = 	0.991541014
        if (JELEM	.eq.	1319	)		coeff = 	0.808209153
        if (JELEM	.eq.	1320	)		coeff = 	1.031265313
        if (JELEM	.eq.	1321	)		coeff = 	0.712595013
        if (JELEM	.eq.	1322	)		coeff = 	0.929114902
        if (JELEM	.eq.	1323	)		coeff = 	1.060523774
        if (JELEM	.eq.	1324	)		coeff = 	0.847393341
        if (JELEM	.eq.	1325	)		coeff = 	0.917355501
        if (JELEM	.eq.	1326	)		coeff = 	0.950880225
        if (JELEM	.eq.	1327	)		coeff = 	1.032806814
        if (JELEM	.eq.	1328	)		coeff = 	1.123874311
        if (JELEM	.eq.	1329	)		coeff = 	1.038220238
        if (JELEM	.eq.	1330	)		coeff = 	1.001841124
        if (JELEM	.eq.	1331	)		coeff = 	0.924592532
        if (JELEM	.eq.	1332	)		coeff = 	0.652545903
        if (JELEM	.eq.	1333	)		coeff = 	1.046425221
        if (JELEM	.eq.	1334	)		coeff = 	0.879930514
        if (JELEM	.eq.	1335	)		coeff = 	0.809392613
        if (JELEM	.eq.	1336	)		coeff = 	0.971679147
        if (JELEM	.eq.	1337	)		coeff = 	1.062839049
        if (JELEM	.eq.	1338	)		coeff = 	0.85533452
        if (JELEM	.eq.	1339	)		coeff = 	0.970504867
        if (JELEM	.eq.	1340	)		coeff = 	1.079637494
        if (JELEM	.eq.	1341	)		coeff = 	0.817127388
        if (JELEM	.eq.	1342	)		coeff = 	0.924038299
        if (JELEM	.eq.	1343	)		coeff = 	1.088274243
        if (JELEM	.eq.	1344	)		coeff = 	0.791114513
        if (JELEM	.eq.	1345	)		coeff = 	1.128573238
        if (JELEM	.eq.	1346	)		coeff = 	1.124229647
        if (JELEM	.eq.	1347	)		coeff = 	0.639826275
        if (JELEM	.eq.	1348	)		coeff = 	0.906962909
        if (JELEM	.eq.	1349	)		coeff = 	0.997607899
        if (JELEM	.eq.	1350	)		coeff = 	0.962807845
        if (JELEM	.eq.	1351	)		coeff = 	0.877226798
        if (JELEM	.eq.	1352	)		coeff = 	1.116809304
        if (JELEM	.eq.	1353	)		coeff = 	0.892393422
        if (JELEM	.eq.	1354	)		coeff = 	0.90156253
        if (JELEM	.eq.	1355	)		coeff = 	0.975327714
        if (JELEM	.eq.	1356	)		coeff = 	0.940376862
        if (JELEM	.eq.	1357	)		coeff = 	1.007864766
        if (JELEM	.eq.	1358	)		coeff = 	1.058656866
        if (JELEM	.eq.	1359	)		coeff = 	0.991514104
        if (JELEM	.eq.	1360	)		coeff = 	0.780248602
        if (JELEM	.eq.	1361	)		coeff = 	1.040475023
        if (JELEM	.eq.	1362	)		coeff = 	1.001870061
        if (JELEM	.eq.	1363	)		coeff = 	1.011873001
        if (JELEM	.eq.	1364	)		coeff = 	1.095163291
        if (JELEM	.eq.	1365	)		coeff = 	0.960475305
        if (JELEM	.eq.	1366	)		coeff = 	0.84375732
        if (JELEM	.eq.	1367	)		coeff = 	0.794524612
        if (JELEM	.eq.	1368	)		coeff = 	0.893285749
        if (JELEM	.eq.	1369	)		coeff = 	0.995898231
        if (JELEM	.eq.	1370	)		coeff = 	1.01929579
        if (JELEM	.eq.	1371	)		coeff = 	0.905011083
        if (JELEM	.eq.	1372	)		coeff = 	0.81374383
        if (JELEM	.eq.	1373	)		coeff = 	0.775281252
        if (JELEM	.eq.	1374	)		coeff = 	1.096404135
        if (JELEM	.eq.	1375	)		coeff = 	0.968808495
        if (JELEM	.eq.	1376	)		coeff = 	1.074614703
        if (JELEM	.eq.	1377	)		coeff = 	1.032336613
        if (JELEM	.eq.	1378	)		coeff = 	0.811168362
        if (JELEM	.eq.	1379	)		coeff = 	1.049957072
        if (JELEM	.eq.	1380	)		coeff = 	1.077478492
        if (JELEM	.eq.	1381	)		coeff = 	0.951673297
        if (JELEM	.eq.	1382	)		coeff = 	1.014636153
        if (JELEM	.eq.	1383	)		coeff = 	0.996156974
        if (JELEM	.eq.	1384	)		coeff = 	0.864728861
        if (JELEM	.eq.	1385	)		coeff = 	0.840206023
        if (JELEM	.eq.	1386	)		coeff = 	0.909037045
        if (JELEM	.eq.	1387	)		coeff = 	0.98671939
        if (JELEM	.eq.	1388	)		coeff = 	0.921554016
        if (JELEM	.eq.	1389	)		coeff = 	1.04421878
        if (JELEM	.eq.	1390	)		coeff = 	0.929704645
        if (JELEM	.eq.	1391	)		coeff = 	0.910778278
        if (JELEM	.eq.	1392	)		coeff = 	0.934721267
        if (JELEM	.eq.	1393	)		coeff = 	1.005883973
        if (JELEM	.eq.	1394	)		coeff = 	1.000945899
        if (JELEM	.eq.	1395	)		coeff = 	1.058356999
        if (JELEM	.eq.	1396	)		coeff = 	0.862930204
        if (JELEM	.eq.	1397	)		coeff = 	0.966047629
        if (JELEM	.eq.	1398	)		coeff = 	1.003610804
        if (JELEM	.eq.	1399	)		coeff = 	0.8722123
        if (JELEM	.eq.	1400	)		coeff = 	1.037174669
        if (JELEM	.eq.	1401	)		coeff = 	0.836931318
        if (JELEM	.eq.	1402	)		coeff = 	0.852480346
        if (JELEM	.eq.	1403	)		coeff = 	0.83610883
        if (JELEM	.eq.	1404	)		coeff = 	0.999372974
        if (JELEM	.eq.	1405	)		coeff = 	0.995834912
        if (JELEM	.eq.	1406	)		coeff = 	0.826795021
        if (JELEM	.eq.	1407	)		coeff = 	0.973957216
        if (JELEM	.eq.	1408	)		coeff = 	0.943861962
        if (JELEM	.eq.	1409	)		coeff = 	0.903730086
        if (JELEM	.eq.	1410	)		coeff = 	0.724538236
        if (JELEM	.eq.	1411	)		coeff = 	0.950920421
        if (JELEM	.eq.	1412	)		coeff = 	0.923592198
        if (JELEM	.eq.	1413	)		coeff = 	0.943760452
        if (JELEM	.eq.	1414	)		coeff = 	0.772500089
        if (JELEM	.eq.	1415	)		coeff = 	0.825220651
        if (JELEM	.eq.	1416	)		coeff = 	1.058933029
        if (JELEM	.eq.	1417	)		coeff = 	1.055855857
        if (JELEM	.eq.	1418	)		coeff = 	1.035148238
        if (JELEM	.eq.	1419	)		coeff = 	0.882128932
        if (JELEM	.eq.	1420	)		coeff = 	1.049021308
        if (JELEM	.eq.	1421	)		coeff = 	0.666146427
        if (JELEM	.eq.	1422	)		coeff = 	0.898523261
        if (JELEM	.eq.	1423	)		coeff = 	1.056980335
        if (JELEM	.eq.	1424	)		coeff = 	0.828731126
        if (JELEM	.eq.	1425	)		coeff = 	0.790229026
        if (JELEM	.eq.	1426	)		coeff = 	0.723115731
        if (JELEM	.eq.	1427	)		coeff = 	0.943846784
        if (JELEM	.eq.	1428	)		coeff = 	0.946113802
        if (JELEM	.eq.	1429	)		coeff = 	1.056556978
        if (JELEM	.eq.	1430	)		coeff = 	0.960191006
        if (JELEM	.eq.	1431	)		coeff = 	0.950298232
        if (JELEM	.eq.	1432	)		coeff = 	1.060548126
        if (JELEM	.eq.	1433	)		coeff = 	0.965707892
        if (JELEM	.eq.	1434	)		coeff = 	0.954128163
        if (JELEM	.eq.	1435	)		coeff = 	1.049128161
        if (JELEM	.eq.	1436	)		coeff = 	0.927855595
        if (JELEM	.eq.	1437	)		coeff = 	1.137840816
        if (JELEM	.eq.	1438	)		coeff = 	1.013474835
        if (JELEM	.eq.	1439	)		coeff = 	0.954731537
        if (JELEM	.eq.	1440	)		coeff = 	1.011260642
        if (JELEM	.eq.	1441	)		coeff = 	0.93436453
        if (JELEM	.eq.	1442	)		coeff = 	1.001617862
        if (JELEM	.eq.	1443	)		coeff = 	1.071937789
        if (JELEM	.eq.	1444	)		coeff = 	0.786115475
        if (JELEM	.eq.	1445	)		coeff = 	0.922308934
        if (JELEM	.eq.	1446	)		coeff = 	0.916331454
        if (JELEM	.eq.	1447	)		coeff = 	0.910730469
        if (JELEM	.eq.	1448	)		coeff = 	0.885068808
        if (JELEM	.eq.	1449	)		coeff = 	0.83848117
        if (JELEM	.eq.	1450	)		coeff = 	0.959354487
        if (JELEM	.eq.	1451	)		coeff = 	1.065430488
        if (JELEM	.eq.	1452	)		coeff = 	0.996530062
        if (JELEM	.eq.	1453	)		coeff = 	0.850144151
        if (JELEM	.eq.	1454	)		coeff = 	1.058397853
        if (JELEM	.eq.	1455	)		coeff = 	1.010375659
        if (JELEM	.eq.	1456	)		coeff = 	0.71331148
        if (JELEM	.eq.	1457	)		coeff = 	0.857642183
        if (JELEM	.eq.	1458	)		coeff = 	1.04167316
        if (JELEM	.eq.	1459	)		coeff = 	0.432099767
        if (JELEM	.eq.	1460	)		coeff = 	1.106572935
        if (JELEM	.eq.	1461	)		coeff = 	0.984402554
        if (JELEM	.eq.	1462	)		coeff = 	0.990122717
        if (JELEM	.eq.	1463	)		coeff = 	0.99121673
        if (JELEM	.eq.	1464	)		coeff = 	1.081541944
        if (JELEM	.eq.	1465	)		coeff = 	0.98471994
        if (JELEM	.eq.	1466	)		coeff = 	0.930864519
        if (JELEM	.eq.	1467	)		coeff = 	0.642623052
        if (JELEM	.eq.	1468	)		coeff = 	1.042387485
        if (JELEM	.eq.	1469	)		coeff = 	1.003758719
        if (JELEM	.eq.	1470	)		coeff = 	1.028412802
        if (JELEM	.eq.	1471	)		coeff = 	1.021142884
        if (JELEM	.eq.	1472	)		coeff = 	1.052572375
        if (JELEM	.eq.	1473	)		coeff = 	1.142169605
        if (JELEM	.eq.	1474	)		coeff = 	0.977912305
        if (JELEM	.eq.	1475	)		coeff = 	1.035107368
        if (JELEM	.eq.	1476	)		coeff = 	0.821935439
        if (JELEM	.eq.	1477	)		coeff = 	0.955978702
        if (JELEM	.eq.	1478	)		coeff = 	0.785807967
        if (JELEM	.eq.	1479	)		coeff = 	0.695322294
        if (JELEM	.eq.	1480	)		coeff = 	0.939440111
        if (JELEM	.eq.	1481	)		coeff = 	1.078485412
        if (JELEM	.eq.	1482	)		coeff = 	0.773136012
        if (JELEM	.eq.	1483	)		coeff = 	0.936990232
        if (JELEM	.eq.	1484	)		coeff = 	0.811425637
        if (JELEM	.eq.	1485	)		coeff = 	0.98467523
        if (JELEM	.eq.	1486	)		coeff = 	0.932809298
        if (JELEM	.eq.	1487	)		coeff = 	1.10229677
        if (JELEM	.eq.	1488	)		coeff = 	0.775041382
        if (JELEM	.eq.	1489	)		coeff = 	0.921816793
        if (JELEM	.eq.	1490	)		coeff = 	1.084884195
        if (JELEM	.eq.	1491	)		coeff = 	0.901961704
        if (JELEM	.eq.	1492	)		coeff = 	0.992430852
        if (JELEM	.eq.	1493	)		coeff = 	1.094468518
        if (JELEM	.eq.	1494	)		coeff = 	1.044093068
        if (JELEM	.eq.	1495	)		coeff = 	1.033561106
        if (JELEM	.eq.	1496	)		coeff = 	1.040287014
        if (JELEM	.eq.	1497	)		coeff = 	0.900954008
        if (JELEM	.eq.	1498	)		coeff = 	0.881149566
        if (JELEM	.eq.	1499	)		coeff = 	0.950639971
        if (JELEM	.eq.	1500	)		coeff = 	0.948849207
        if (JELEM	.eq.	1501	)		coeff = 	0.925511938
        if (JELEM	.eq.	1502	)		coeff = 	0.65551597
        if (JELEM	.eq.	1503	)		coeff = 	0.924329206
        if (JELEM	.eq.	1504	)		coeff = 	0.934894673
        if (JELEM	.eq.	1505	)		coeff = 	0.790434272
        if (JELEM	.eq.	1506	)		coeff = 	0.943776924
        if (JELEM	.eq.	1507	)		coeff = 	1.008874296
        if (JELEM	.eq.	1508	)		coeff = 	0.731398125
        if (JELEM	.eq.	1509	)		coeff = 	0.980488557
        if (JELEM	.eq.	1510	)		coeff = 	0.934562642
        if (JELEM	.eq.	1511	)		coeff = 	0.894552602
        if (JELEM	.eq.	1512	)		coeff = 	0.909542005
        if (JELEM	.eq.	1513	)		coeff = 	1.044637799
        if (JELEM	.eq.	1514	)		coeff = 	1.096205206
        if (JELEM	.eq.	1515	)		coeff = 	1.026019912
        if (JELEM	.eq.	1516	)		coeff = 	0.822744682
        if (JELEM	.eq.	1517	)		coeff = 	0.947142152
        if (JELEM	.eq.	1518	)		coeff = 	0.973760771
        if (JELEM	.eq.	1519	)		coeff = 	0.98309773
        if (JELEM	.eq.	1520	)		coeff = 	0.87268635
        if (JELEM	.eq.	1521	)		coeff = 	0.917942218
        if (JELEM	.eq.	1522	)		coeff = 	0.916679355
        if (JELEM	.eq.	1523	)		coeff = 	1.062075364
        if (JELEM	.eq.	1524	)		coeff = 	0.982535205
        if (JELEM	.eq.	1525	)		coeff = 	0.962595327
        if (JELEM	.eq.	1526	)		coeff = 	0.997976021
        if (JELEM	.eq.	1527	)		coeff = 	0.969442132
        if (JELEM	.eq.	1528	)		coeff = 	1.006949881
        if (JELEM	.eq.	1529	)		coeff = 	0.871288227
        if (JELEM	.eq.	1530	)		coeff = 	0.99564139
        if (JELEM	.eq.	1531	)		coeff = 	0.897813778
        if (JELEM	.eq.	1532	)		coeff = 	0.969287648
        if (JELEM	.eq.	1533	)		coeff = 	0.891139999
        if (JELEM	.eq.	1534	)		coeff = 	0.760220122
        if (JELEM	.eq.	1535	)		coeff = 	0.959186521
        if (JELEM	.eq.	1536	)		coeff = 	0.795868896
        if (JELEM	.eq.	1537	)		coeff = 	1.042922938
        if (JELEM	.eq.	1538	)		coeff = 	0.818871735
        if (JELEM	.eq.	1539	)		coeff = 	1.095634319
        if (JELEM	.eq.	1540	)		coeff = 	0.973548666
        if (JELEM	.eq.	1541	)		coeff = 	1.143418693
        if (JELEM	.eq.	1542	)		coeff = 	0.8566562
        if (JELEM	.eq.	1543	)		coeff = 	1.055680827
        if (JELEM	.eq.	1544	)		coeff = 	1.06050721
        if (JELEM	.eq.	1545	)		coeff = 	1.054863131
        if (JELEM	.eq.	1546	)		coeff = 	0.905108697
        if (JELEM	.eq.	1547	)		coeff = 	1.044319419
        if (JELEM	.eq.	1548	)		coeff = 	1.019270984
        if (JELEM	.eq.	1549	)		coeff = 	0.875150777
        if (JELEM	.eq.	1550	)		coeff = 	0.96369222
        if (JELEM	.eq.	1551	)		coeff = 	0.790180218
        if (JELEM	.eq.	1552	)		coeff = 	1.110407088
        if (JELEM	.eq.	1553	)		coeff = 	0.981348991
        if (JELEM	.eq.	1554	)		coeff = 	0.943358418
        if (JELEM	.eq.	1555	)		coeff = 	0.945473171
        if (JELEM	.eq.	1556	)		coeff = 	0.848666448
        if (JELEM	.eq.	1557	)		coeff = 	1.075502674
        if (JELEM	.eq.	1558	)		coeff = 	1.018703098
        if (JELEM	.eq.	1559	)		coeff = 	1.199293608
        if (JELEM	.eq.	1560	)		coeff = 	0.741324454
        if (JELEM	.eq.	1561	)		coeff = 	0.876027882
        if (JELEM	.eq.	1562	)		coeff = 	0.882324229
        if (JELEM	.eq.	1563	)		coeff = 	1.070385722
        if (JELEM	.eq.	1564	)		coeff = 	1.005060032
        if (JELEM	.eq.	1565	)		coeff = 	1.065627254
        if (JELEM	.eq.	1566	)		coeff = 	0.964923514
        if (JELEM	.eq.	1567	)		coeff = 	0.856470976
        if (JELEM	.eq.	1568	)		coeff = 	0.924779865
        if (JELEM	.eq.	1569	)		coeff = 	0.906182057
        if (JELEM	.eq.	1570	)		coeff = 	0.922628265
        if (JELEM	.eq.	1571	)		coeff = 	0.891085234
        if (JELEM	.eq.	1572	)		coeff = 	0.82774842
        if (JELEM	.eq.	1573	)		coeff = 	0.926660217
        if (JELEM	.eq.	1574	)		coeff = 	1.0552082
        if (JELEM	.eq.	1575	)		coeff = 	0.943047605
        if (JELEM	.eq.	1576	)		coeff = 	1.061159724
        if (JELEM	.eq.	1577	)		coeff = 	0.793256454
        if (JELEM	.eq.	1578	)		coeff = 	1.098539428
        if (JELEM	.eq.	1579	)		coeff = 	1.00801736
        if (JELEM	.eq.	1580	)		coeff = 	0.940877273
        if (JELEM	.eq.	1581	)		coeff = 	0.970854424
        if (JELEM	.eq.	1582	)		coeff = 	0.858121516
        if (JELEM	.eq.	1583	)		coeff = 	0.955363337
        if (JELEM	.eq.	1584	)		coeff = 	1.040083915
        if (JELEM	.eq.	1585	)		coeff = 	0.898579998
        if (JELEM	.eq.	1586	)		coeff = 	1.066641866
        if (JELEM	.eq.	1587	)		coeff = 	0.916550603
        if (JELEM	.eq.	1588	)		coeff = 	0.924423823
        if (JELEM	.eq.	1589	)		coeff = 	1.039468505
        if (JELEM	.eq.	1590	)		coeff = 	1.054650648
        if (JELEM	.eq.	1591	)		coeff = 	1.060165527
        if (JELEM	.eq.	1592	)		coeff = 	1.066275624
        if (JELEM	.eq.	1593	)		coeff = 	1.047846565
        if (JELEM	.eq.	1594	)		coeff = 	0.735119291
        if (JELEM	.eq.	1595	)		coeff = 	1.152682231
        if (JELEM	.eq.	1596	)		coeff = 	0.730741128
        if (JELEM	.eq.	1597	)		coeff = 	1.138588175
        if (JELEM	.eq.	1598	)		coeff = 	0.702601004
        if (JELEM	.eq.	1599	)		coeff = 	1.019414862
        if (JELEM	.eq.	1600	)		coeff = 	0.956980909
        if (JELEM	.eq.	1601	)		coeff = 	0.826147194
        if (JELEM	.eq.	1602	)		coeff = 	0.801469626
        if (JELEM	.eq.	1603	)		coeff = 	1.052364669
        if (JELEM	.eq.	1604	)		coeff = 	0.915682003
        if (JELEM	.eq.	1605	)		coeff = 	0.755480795
        if (JELEM	.eq.	1606	)		coeff = 	0.690368269
        if (JELEM	.eq.	1607	)		coeff = 	1.083320299
        if (JELEM	.eq.	1608	)		coeff = 	1.055787185
        if (JELEM	.eq.	1609	)		coeff = 	0.884263117
        if (JELEM	.eq.	1610	)		coeff = 	1.116361747
        if (JELEM	.eq.	1611	)		coeff = 	1.101990899
        if (JELEM	.eq.	1612	)		coeff = 	0.967027606
        if (JELEM	.eq.	1613	)		coeff = 	0.833918802
        if (JELEM	.eq.	1614	)		coeff = 	0.55869181
        if (JELEM	.eq.	1615	)		coeff = 	1.184302463
        if (JELEM	.eq.	1616	)		coeff = 	0.951985536
        if (JELEM	.eq.	1617	)		coeff = 	0.826786386
        if (JELEM	.eq.	1618	)		coeff = 	0.790478034
        if (JELEM	.eq.	1619	)		coeff = 	0.836627101
        if (JELEM	.eq.	1620	)		coeff = 	0.81488753
        if (JELEM	.eq.	1621	)		coeff = 	0.884412648
        if (JELEM	.eq.	1622	)		coeff = 	1.079120464
        if (JELEM	.eq.	1623	)		coeff = 	0.961495664
        if (JELEM	.eq.	1624	)		coeff = 	1.059287859
        if (JELEM	.eq.	1625	)		coeff = 	0.844724584
        if (JELEM	.eq.	1626	)		coeff = 	0.77150536
        if (JELEM	.eq.	1627	)		coeff = 	1.059057876
        if (JELEM	.eq.	1628	)		coeff = 	0.811335043
        if (JELEM	.eq.	1629	)		coeff = 	0.994580037
        if (JELEM	.eq.	1630	)		coeff = 	0.996013754
        if (JELEM	.eq.	1631	)		coeff = 	1.026896432
        if (JELEM	.eq.	1632	)		coeff = 	0.822500998
        if (JELEM	.eq.	1633	)		coeff = 	0.886294947
        if (JELEM	.eq.	1634	)		coeff = 	0.978328863
        if (JELEM	.eq.	1635	)		coeff = 	0.898523203
        if (JELEM	.eq.	1636	)		coeff = 	0.75127505
        if (JELEM	.eq.	1637	)		coeff = 	1.057437502
        if (JELEM	.eq.	1638	)		coeff = 	1.034819911
        if (JELEM	.eq.	1639	)		coeff = 	0.9221892
        if (JELEM	.eq.	1640	)		coeff = 	0.856490794
        if (JELEM	.eq.	1641	)		coeff = 	0.831801321
        if (JELEM	.eq.	1642	)		coeff = 	0.991805077
        if (JELEM	.eq.	1643	)		coeff = 	1.079980657
        if (JELEM	.eq.	1644	)		coeff = 	1.096928566
        if (JELEM	.eq.	1645	)		coeff = 	1.00201691
        if (JELEM	.eq.	1646	)		coeff = 	0.845930415
        if (JELEM	.eq.	1647	)		coeff = 	1.044051702
        if (JELEM	.eq.	1648	)		coeff = 	0.864961813
        if (JELEM	.eq.	1649	)		coeff = 	0.917688282
        if (JELEM	.eq.	1650	)		coeff = 	1.138042121
        if (JELEM	.eq.	1651	)		coeff = 	0.867394544
        if (JELEM	.eq.	1652	)		coeff = 	0.777379832
        if (JELEM	.eq.	1653	)		coeff = 	0.966145612
        if (JELEM	.eq.	1654	)		coeff = 	0.843120315
        if (JELEM	.eq.	1655	)		coeff = 	1.07336306
        if (JELEM	.eq.	1656	)		coeff = 	0.878798509
        if (JELEM	.eq.	1657	)		coeff = 	0.774002907
        if (JELEM	.eq.	1658	)		coeff = 	0.843844396
        if (JELEM	.eq.	1659	)		coeff = 	1.030417941
        if (JELEM	.eq.	1660	)		coeff = 	1.044554457
        if (JELEM	.eq.	1661	)		coeff = 	0.957760456
        if (JELEM	.eq.	1662	)		coeff = 	0.992010802
        if (JELEM	.eq.	1663	)		coeff = 	0.969789515
        if (JELEM	.eq.	1664	)		coeff = 	0.600675311
        if (JELEM	.eq.	1665	)		coeff = 	0.933663526
        if (JELEM	.eq.	1666	)		coeff = 	0.750469493
        if (JELEM	.eq.	1667	)		coeff = 	0.966663468
        if (JELEM	.eq.	1668	)		coeff = 	0.98102492
        if (JELEM	.eq.	1669	)		coeff = 	0.873270566
        if (JELEM	.eq.	1670	)		coeff = 	0.885266608
        if (JELEM	.eq.	1671	)		coeff = 	0.9796792
        if (JELEM	.eq.	1672	)		coeff = 	1.113774399
        if (JELEM	.eq.	1673	)		coeff = 	1.092960089
        if (JELEM	.eq.	1674	)		coeff = 	0.861680534
        if (JELEM	.eq.	1675	)		coeff = 	0.917392202
        if (JELEM	.eq.	1676	)		coeff = 	1.131226351
        if (JELEM	.eq.	1677	)		coeff = 	0.947808912
        if (JELEM	.eq.	1678	)		coeff = 	0.894726916
        if (JELEM	.eq.	1679	)		coeff = 	1.082212957
        if (JELEM	.eq.	1680	)		coeff = 	1.043408813
        if (JELEM	.eq.	1681	)		coeff = 	0.85528447
        if (JELEM	.eq.	1682	)		coeff = 	1.070464993
        if (JELEM	.eq.	1683	)		coeff = 	0.812691576
        if (JELEM	.eq.	1684	)		coeff = 	0.776335953
        if (JELEM	.eq.	1685	)		coeff = 	1.158687851
        if (JELEM	.eq.	1686	)		coeff = 	0.997480399
        if (JELEM	.eq.	1687	)		coeff = 	1.059648191
        if (JELEM	.eq.	1688	)		coeff = 	0.952669167
        if (JELEM	.eq.	1689	)		coeff = 	1.086199187
        if (JELEM	.eq.	1690	)		coeff = 	1.124654149
        if (JELEM	.eq.	1691	)		coeff = 	0.765509633
        if (JELEM	.eq.	1692	)		coeff = 	0.701404501
        if (JELEM	.eq.	1693	)		coeff = 	1.001891492
        if (JELEM	.eq.	1694	)		coeff = 	0.921126688
        if (JELEM	.eq.	1695	)		coeff = 	1.103975535
        if (JELEM	.eq.	1696	)		coeff = 	1.046183042
        if (JELEM	.eq.	1697	)		coeff = 	1.124357717
        if (JELEM	.eq.	1698	)		coeff = 	0.972465912
        if (JELEM	.eq.	1699	)		coeff = 	1.066092089
        if (JELEM	.eq.	1700	)		coeff = 	0.622453362
        if (JELEM	.eq.	1701	)		coeff = 	0.983974939
        if (JELEM	.eq.	1702	)		coeff = 	0.734426337
        if (JELEM	.eq.	1703	)		coeff = 	0.893037091
        if (JELEM	.eq.	1704	)		coeff = 	0.940795541
        if (JELEM	.eq.	1705	)		coeff = 	0.952659111
        if (JELEM	.eq.	1706	)		coeff = 	0.900065393
        if (JELEM	.eq.	1707	)		coeff = 	1.181317281
        if (JELEM	.eq.	1708	)		coeff = 	0.868879082
        if (JELEM	.eq.	1709	)		coeff = 	0.772777757
        if (JELEM	.eq.	1710	)		coeff = 	1.169639133
        if (JELEM	.eq.	1711	)		coeff = 	0.848234179
        if (JELEM	.eq.	1712	)		coeff = 	0.875582008
        if (JELEM	.eq.	1713	)		coeff = 	0.55685993
        if (JELEM	.eq.	1714	)		coeff = 	1.039978087
        if (JELEM	.eq.	1715	)		coeff = 	0.780486684
        if (JELEM	.eq.	1716	)		coeff = 	0.921842751
        if (JELEM	.eq.	1717	)		coeff = 	1.084507925
        if (JELEM	.eq.	1718	)		coeff = 	1.0278533
        if (JELEM	.eq.	1719	)		coeff = 	0.877078511
        if (JELEM	.eq.	1720	)		coeff = 	0.858028136
        if (JELEM	.eq.	1721	)		coeff = 	1.084997659
        if (JELEM	.eq.	1722	)		coeff = 	0.972356482
        if (JELEM	.eq.	1723	)		coeff = 	1.04265024
        if (JELEM	.eq.	1724	)		coeff = 	0.777228265
        if (JELEM	.eq.	1725	)		coeff = 	1.013049117
        if (JELEM	.eq.	1726	)		coeff = 	0.829206287
        if (JELEM	.eq.	1727	)		coeff = 	1.030286831
        if (JELEM	.eq.	1728	)		coeff = 	0.815467569
        if (JELEM	.eq.	1729	)		coeff = 	1.052613191
        if (JELEM	.eq.	1730	)		coeff = 	0.879046539
        if (JELEM	.eq.	1731	)		coeff = 	1.131890104
        if (JELEM	.eq.	1732	)		coeff = 	0.921732607
        if (JELEM	.eq.	1733	)		coeff = 	0.944943952
        if (JELEM	.eq.	1734	)		coeff = 	0.997683757
        if (JELEM	.eq.	1735	)		coeff = 	1.044704639
        if (JELEM	.eq.	1736	)		coeff = 	0.864439938
        if (JELEM	.eq.	1737	)		coeff = 	1.067849321
        if (JELEM	.eq.	1738	)		coeff = 	0.967022926
        if (JELEM	.eq.	1739	)		coeff = 	1.15850056
        if (JELEM	.eq.	1740	)		coeff = 	1.053168639
        if (JELEM	.eq.	1741	)		coeff = 	0.968103339
        if (JELEM	.eq.	1742	)		coeff = 	0.840749189
        if (JELEM	.eq.	1743	)		coeff = 	1.069533483
        if (JELEM	.eq.	1744	)		coeff = 	0.889960369
        if (JELEM	.eq.	1745	)		coeff = 	0.905239785
        if (JELEM	.eq.	1746	)		coeff = 	1.129079284
        if (JELEM	.eq.	1747	)		coeff = 	0.967101763
        if (JELEM	.eq.	1748	)		coeff = 	0.701907605
        if (JELEM	.eq.	1749	)		coeff = 	1.08130896
        if (JELEM	.eq.	1750	)		coeff = 	0.885609421
        if (JELEM	.eq.	1751	)		coeff = 	0.922960059
        if (JELEM	.eq.	1752	)		coeff = 	0.936801878
        if (JELEM	.eq.	1753	)		coeff = 	0.964395234
        if (JELEM	.eq.	1754	)		coeff = 	0.944655494
        if (JELEM	.eq.	1755	)		coeff = 	0.984126068
        if (JELEM	.eq.	1756	)		coeff = 	1.099154564
        if (JELEM	.eq.	1757	)		coeff = 	1.021401985
        if (JELEM	.eq.	1758	)		coeff = 	0.946565093
        if (JELEM	.eq.	1759	)		coeff = 	0.924631393
        if (JELEM	.eq.	1760	)		coeff = 	0.768376705
        if (JELEM	.eq.	1761	)		coeff = 	0.68419967
        if (JELEM	.eq.	1762	)		coeff = 	1.090064356
        if (JELEM	.eq.	1763	)		coeff = 	0.915342139
        if (JELEM	.eq.	1764	)		coeff = 	0.934186273
        if (JELEM	.eq.	1765	)		coeff = 	0.971243495
        if (JELEM	.eq.	1766	)		coeff = 	1.003165014
        if (JELEM	.eq.	1767	)		coeff = 	0.97076015
        if (JELEM	.eq.	1768	)		coeff = 	0.912204498
        if (JELEM	.eq.	1769	)		coeff = 	0.72691046
        if (JELEM	.eq.	1770	)		coeff = 	1.092317091
        if (JELEM	.eq.	1771	)		coeff = 	0.861796996
        if (JELEM	.eq.	1772	)		coeff = 	0.937831653
        if (JELEM	.eq.	1773	)		coeff = 	0.787646503
        if (JELEM	.eq.	1774	)		coeff = 	1.086448743
        if (JELEM	.eq.	1775	)		coeff = 	1.020632755
        if (JELEM	.eq.	1776	)		coeff = 	1.114793921
        if (JELEM	.eq.	1777	)		coeff = 	0.962862524
        if (JELEM	.eq.	1778	)		coeff = 	0.875133226
        if (JELEM	.eq.	1779	)		coeff = 	1.023571829
        if (JELEM	.eq.	1780	)		coeff = 	1.040687245
        if (JELEM	.eq.	1781	)		coeff = 	1.01001281
        if (JELEM	.eq.	1782	)		coeff = 	0.976865069
        if (JELEM	.eq.	1783	)		coeff = 	0.887931564
        if (JELEM	.eq.	1784	)		coeff = 	0.961254078
        if (JELEM	.eq.	1785	)		coeff = 	0.996028729
        if (JELEM	.eq.	1786	)		coeff = 	0.793759683
        if (JELEM	.eq.	1787	)		coeff = 	0.715898361
        if (JELEM	.eq.	1788	)		coeff = 	0.926244472
        if (JELEM	.eq.	1789	)		coeff = 	1.073090718
        if (JELEM	.eq.	1790	)		coeff = 	0.929382553
        if (JELEM	.eq.	1791	)		coeff = 	0.995891653
        if (JELEM	.eq.	1792	)		coeff = 	0.623262567
        if (JELEM	.eq.	1793	)		coeff = 	1.022473118
        if (JELEM	.eq.	1794	)		coeff = 	0.899790833
        if (JELEM	.eq.	1795	)		coeff = 	0.954102888
        if (JELEM	.eq.	1796	)		coeff = 	1.050967891
        if (JELEM	.eq.	1797	)		coeff = 	0.905817275
        if (JELEM	.eq.	1798	)		coeff = 	1.115626375
        if (JELEM	.eq.	1799	)		coeff = 	1.053902819
        if (JELEM	.eq.	1800	)		coeff = 	1.119302836
        if (JELEM	.eq.	1801	)		coeff = 	0.810361767
        if (JELEM	.eq.	1802	)		coeff = 	0.839854348
        if (JELEM	.eq.	1803	)		coeff = 	1.078835332
        if (JELEM	.eq.	1804	)		coeff = 	0.988479405
        if (JELEM	.eq.	1805	)		coeff = 	1.077953397
        if (JELEM	.eq.	1806	)		coeff = 	0.943406324
        if (JELEM	.eq.	1807	)		coeff = 	0.743954546
        if (JELEM	.eq.	1808	)		coeff = 	1.031307527
        if (JELEM	.eq.	1809	)		coeff = 	0.632122887
        if (JELEM	.eq.	1810	)		coeff = 	1.004926413
        if (JELEM	.eq.	1811	)		coeff = 	1.045987618
        if (JELEM	.eq.	1812	)		coeff = 	0.913966135
        if (JELEM	.eq.	1813	)		coeff = 	0.696889133
        if (JELEM	.eq.	1814	)		coeff = 	0.928003226
        if (JELEM	.eq.	1815	)		coeff = 	1.106699704
        if (JELEM	.eq.	1816	)		coeff = 	0.998470693
        if (JELEM	.eq.	1817	)		coeff = 	1.060201409
        if (JELEM	.eq.	1818	)		coeff = 	1.038862566
        if (JELEM	.eq.	1819	)		coeff = 	1.11434354
        if (JELEM	.eq.	1820	)		coeff = 	0.796996941
        if (JELEM	.eq.	1821	)		coeff = 	0.863907326
        if (JELEM	.eq.	1822	)		coeff = 	0.998605122
        if (JELEM	.eq.	1823	)		coeff = 	0.844205229
        if (JELEM	.eq.	1824	)		coeff = 	0.881281197
        if (JELEM	.eq.	1825	)		coeff = 	0.928266311
        if (JELEM	.eq.	1826	)		coeff = 	0.992893802
        if (JELEM	.eq.	1827	)		coeff = 	1.002341075
        if (JELEM	.eq.	1828	)		coeff = 	1.092434726
        if (JELEM	.eq.	1829	)		coeff = 	1.007154173
        if (JELEM	.eq.	1830	)		coeff = 	0.950240913
        if (JELEM	.eq.	1831	)		coeff = 	0.974875561
        if (JELEM	.eq.	1832	)		coeff = 	0.920720808
        if (JELEM	.eq.	1833	)		coeff = 	0.9602354
        if (JELEM	.eq.	1834	)		coeff = 	0.85355751
        if (JELEM	.eq.	1835	)		coeff = 	1.088321375
        if (JELEM	.eq.	1836	)		coeff = 	0.974010689
        if (JELEM	.eq.	1837	)		coeff = 	0.938120512
        if (JELEM	.eq.	1838	)		coeff = 	1.05298049
        if (JELEM	.eq.	1839	)		coeff = 	0.931551643
        if (JELEM	.eq.	1840	)		coeff = 	1.114536649
        if (JELEM	.eq.	1841	)		coeff = 	0.94233102
        if (JELEM	.eq.	1842	)		coeff = 	0.838384255
        if (JELEM	.eq.	1843	)		coeff = 	0.964088451
        if (JELEM	.eq.	1844	)		coeff = 	0.980732331
        if (JELEM	.eq.	1845	)		coeff = 	1.066482168
        if (JELEM	.eq.	1846	)		coeff = 	1.135568937
        if (JELEM	.eq.	1847	)		coeff = 	0.880103382
        if (JELEM	.eq.	1848	)		coeff = 	0.862576847
        if (JELEM	.eq.	1849	)		coeff = 	1.02055873
        if (JELEM	.eq.	1850	)		coeff = 	1.080105105
        if (JELEM	.eq.	1851	)		coeff = 	0.998041319
        if (JELEM	.eq.	1852	)		coeff = 	0.845936835
        if (JELEM	.eq.	1853	)		coeff = 	0.838713365
        if (JELEM	.eq.	1854	)		coeff = 	0.914164715
        if (JELEM	.eq.	1855	)		coeff = 	0.725995676
        if (JELEM	.eq.	1856	)		coeff = 	0.752944669
        if (JELEM	.eq.	1857	)		coeff = 	1.081199852
        if (JELEM	.eq.	1858	)		coeff = 	0.919794914
        if (JELEM	.eq.	1859	)		coeff = 	0.969323261
        if (JELEM	.eq.	1860	)		coeff = 	1.104843616
        if (JELEM	.eq.	1861	)		coeff = 	0.800346388
        if (JELEM	.eq.	1862	)		coeff = 	0.964782059
        if (JELEM	.eq.	1863	)		coeff = 	0.873861793
        if (JELEM	.eq.	1864	)		coeff = 	1.108741495
        if (JELEM	.eq.	1865	)		coeff = 	1.029519974
        if (JELEM	.eq.	1866	)		coeff = 	0.918882133
        if (JELEM	.eq.	1867	)		coeff = 	1.072461466
        if (JELEM	.eq.	1868	)		coeff = 	0.922958652
        if (JELEM	.eq.	1869	)		coeff = 	0.995370721
        if (JELEM	.eq.	1870	)		coeff = 	0.876285859
        if (JELEM	.eq.	1871	)		coeff = 	0.918273973
        if (JELEM	.eq.	1872	)		coeff = 	0.996306522
        if (JELEM	.eq.	1873	)		coeff = 	1.018730792
        if (JELEM	.eq.	1874	)		coeff = 	1.007576336
        if (JELEM	.eq.	1875	)		coeff = 	0.781107662
        if (JELEM	.eq.	1876	)		coeff = 	0.976043382
        if (JELEM	.eq.	1877	)		coeff = 	0.979789769
        if (JELEM	.eq.	1878	)		coeff = 	0.976610517
        if (JELEM	.eq.	1879	)		coeff = 	0.749944656
        if (JELEM	.eq.	1880	)		coeff = 	1.042632354
        if (JELEM	.eq.	1881	)		coeff = 	0.812321666
        if (JELEM	.eq.	1882	)		coeff = 	1.146329452
        if (JELEM	.eq.	1883	)		coeff = 	1.007130052
        if (JELEM	.eq.	1884	)		coeff = 	0.876144615
        if (JELEM	.eq.	1885	)		coeff = 	1.006844992
        if (JELEM	.eq.	1886	)		coeff = 	0.929224729
        if (JELEM	.eq.	1887	)		coeff = 	0.976927762
        if (JELEM	.eq.	1888	)		coeff = 	1.164587679
        if (JELEM	.eq.	1889	)		coeff = 	0.935313855
        if (JELEM	.eq.	1890	)		coeff = 	0.934554769
        if (JELEM	.eq.	1891	)		coeff = 	0.919424659
        if (JELEM	.eq.	1892	)		coeff = 	1.006867781
        if (JELEM	.eq.	1893	)		coeff = 	0.965872447
        if (JELEM	.eq.	1894	)		coeff = 	0.901400888
        if (JELEM	.eq.	1895	)		coeff = 	0.808270696
        if (JELEM	.eq.	1896	)		coeff = 	1.112321617
        if (JELEM	.eq.	1897	)		coeff = 	1.08775721
        if (JELEM	.eq.	1898	)		coeff = 	0.919304804
        if (JELEM	.eq.	1899	)		coeff = 	0.876981668
        if (JELEM	.eq.	1900	)		coeff = 	0.643175149
        if (JELEM	.eq.	1901	)		coeff = 	1.075822051
        if (JELEM	.eq.	1902	)		coeff = 	1.000925273
        if (JELEM	.eq.	1903	)		coeff = 	0.910421101
        if (JELEM	.eq.	1904	)		coeff = 	0.997860468
        if (JELEM	.eq.	1905	)		coeff = 	0.825420787
        if (JELEM	.eq.	1906	)		coeff = 	1.021006875
        if (JELEM	.eq.	1907	)		coeff = 	1.072512892
        if (JELEM	.eq.	1908	)		coeff = 	0.911648742
        if (JELEM	.eq.	1909	)		coeff = 	1.047899176
        if (JELEM	.eq.	1910	)		coeff = 	0.822077439
        if (JELEM	.eq.	1911	)		coeff = 	0.882384524
        if (JELEM	.eq.	1912	)		coeff = 	0.986054484
        if (JELEM	.eq.	1913	)		coeff = 	1.236764831
        if (JELEM	.eq.	1914	)		coeff = 	1.066329421
        if (JELEM	.eq.	1915	)		coeff = 	1.026206497
        if (JELEM	.eq.	1916	)		coeff = 	0.819419784
        if (JELEM	.eq.	1917	)		coeff = 	0.934651714
        if (JELEM	.eq.	1918	)		coeff = 	1.012808656
        if (JELEM	.eq.	1919	)		coeff = 	1.02319738
        if (JELEM	.eq.	1920	)		coeff = 	0.981738471
        if (JELEM	.eq.	1921	)		coeff = 	0.795270419
        if (JELEM	.eq.	1922	)		coeff = 	0.774683625
        if (JELEM	.eq.	1923	)		coeff = 	0.962537742
        if (JELEM	.eq.	1924	)		coeff = 	0.926461812
        if (JELEM	.eq.	1925	)		coeff = 	0.894930046
        if (JELEM	.eq.	1926	)		coeff = 	1.140800574
        if (JELEM	.eq.	1927	)		coeff = 	0.942566432
        if (JELEM	.eq.	1928	)		coeff = 	1.118611645
        if (JELEM	.eq.	1929	)		coeff = 	0.985201585
        if (JELEM	.eq.	1930	)		coeff = 	0.972910903
        if (JELEM	.eq.	1931	)		coeff = 	1.142477106
        if (JELEM	.eq.	1932	)		coeff = 	1.105736685
        if (JELEM	.eq.	1933	)		coeff = 	0.7759119
        if (JELEM	.eq.	1934	)		coeff = 	0.954405914
        if (JELEM	.eq.	1935	)		coeff = 	1.000294426
        if (JELEM	.eq.	1936	)		coeff = 	1.001069764
        if (JELEM	.eq.	1937	)		coeff = 	1.06561491
        if (JELEM	.eq.	1938	)		coeff = 	1.066275059
        if (JELEM	.eq.	1939	)		coeff = 	1.004654592
        if (JELEM	.eq.	1940	)		coeff = 	1.008723931
        if (JELEM	.eq.	1941	)		coeff = 	0.868192697
        if (JELEM	.eq.	1942	)		coeff = 	0.967686909
        if (JELEM	.eq.	1943	)		coeff = 	0.973714223
        if (JELEM	.eq.	1944	)		coeff = 	1.073403782
        if (JELEM	.eq.	1945	)		coeff = 	0.809344652
        if (JELEM	.eq.	1946	)		coeff = 	0.910995795
        if (JELEM	.eq.	1947	)		coeff = 	0.842479188
        if (JELEM	.eq.	1948	)		coeff = 	0.917091847
        if (JELEM	.eq.	1949	)		coeff = 	0.662198772
        if (JELEM	.eq.	1950	)		coeff = 	0.677640735
        if (JELEM	.eq.	1951	)		coeff = 	1.03315712
        if (JELEM	.eq.	1952	)		coeff = 	0.927413824
        if (JELEM	.eq.	1953	)		coeff = 	0.891502447
        if (JELEM	.eq.	1954	)		coeff = 	0.964520896
        if (JELEM	.eq.	1955	)		coeff = 	0.833958326
        if (JELEM	.eq.	1956	)		coeff = 	1.051728586
        if (JELEM	.eq.	1957	)		coeff = 	1.076311509
        if (JELEM	.eq.	1958	)		coeff = 	1.193848651
        if (JELEM	.eq.	1959	)		coeff = 	1.065027712
        if (JELEM	.eq.	1960	)		coeff = 	0.954400797
        if (JELEM	.eq.	1961	)		coeff = 	0.961044684
        if (JELEM	.eq.	1962	)		coeff = 	0.995296837
        if (JELEM	.eq.	1963	)		coeff = 	1.015760075
        if (JELEM	.eq.	1964	)		coeff = 	1.188837025
        if (JELEM	.eq.	1965	)		coeff = 	0.853140191
        if (JELEM	.eq.	1966	)		coeff = 	0.922999661
        if (JELEM	.eq.	1967	)		coeff = 	0.978201952
        if (JELEM	.eq.	1968	)		coeff = 	1.034973351
        if (JELEM	.eq.	1969	)		coeff = 	0.859070422
        if (JELEM	.eq.	1970	)		coeff = 	0.848584392
        if (JELEM	.eq.	1971	)		coeff = 	0.832524589
        if (JELEM	.eq.	1972	)		coeff = 	0.973042973
        if (JELEM	.eq.	1973	)		coeff = 	0.703610732
        if (JELEM	.eq.	1974	)		coeff = 	0.839021028
        if (JELEM	.eq.	1975	)		coeff = 	1.097873849
        if (JELEM	.eq.	1976	)		coeff = 	1.036933382
        if (JELEM	.eq.	1977	)		coeff = 	0.851941468
        if (JELEM	.eq.	1978	)		coeff = 	0.989719628
        if (JELEM	.eq.	1979	)		coeff = 	0.973300069
        if (JELEM	.eq.	1980	)		coeff = 	0.740655067
        if (JELEM	.eq.	1981	)		coeff = 	0.716426166
        if (JELEM	.eq.	1982	)		coeff = 	0.876462901
        if (JELEM	.eq.	1983	)		coeff = 	0.942684562
        if (JELEM	.eq.	1984	)		coeff = 	0.784078821
        if (JELEM	.eq.	1985	)		coeff = 	0.965281032
        if (JELEM	.eq.	1986	)		coeff = 	1.06028644
        if (JELEM	.eq.	1987	)		coeff = 	1.011475418
        if (JELEM	.eq.	1988	)		coeff = 	1.019742754
        if (JELEM	.eq.	1989	)		coeff = 	0.947450521
        if (JELEM	.eq.	1990	)		coeff = 	1.104260476
        if (JELEM	.eq.	1991	)		coeff = 	1.103355975
        if (JELEM	.eq.	1992	)		coeff = 	1.060012105
        if (JELEM	.eq.	1993	)		coeff = 	0.746851315
        if (JELEM	.eq.	1994	)		coeff = 	0.855250128
        if (JELEM	.eq.	1995	)		coeff = 	0.898221389
        if (JELEM	.eq.	1996	)		coeff = 	0.704708389
        if (JELEM	.eq.	1997	)		coeff = 	0.524460606
        if (JELEM	.eq.	1998	)		coeff = 	0.645835013
        if (JELEM	.eq.	1999	)		coeff = 	1.066095889
        if (JELEM	.eq.	2000	)		coeff = 	0.729028665
        if (JELEM	.eq.	2001	)		coeff = 	0.955454001
        if (JELEM	.eq.	2002	)		coeff = 	1.10037176
        if (JELEM	.eq.	2003	)		coeff = 	1.015417758
        if (JELEM	.eq.	2004	)		coeff = 	0.802473595
        if (JELEM	.eq.	2005	)		coeff = 	0.84269813
        if (JELEM	.eq.	2006	)		coeff = 	1.197300892
        if (JELEM	.eq.	2007	)		coeff = 	0.922421718
        if (JELEM	.eq.	2008	)		coeff = 	0.859169299
        if (JELEM	.eq.	2009	)		coeff = 	1.034659394
        if (JELEM	.eq.	2010	)		coeff = 	1.106331554
        if (JELEM	.eq.	2011	)		coeff = 	1.029317473
        if (JELEM	.eq.	2012	)		coeff = 	1.085706301
        if (JELEM	.eq.	2013	)		coeff = 	0.96852019
        if (JELEM	.eq.	2014	)		coeff = 	0.986188465
        if (JELEM	.eq.	2015	)		coeff = 	0.99636074
        if (JELEM	.eq.	2016	)		coeff = 	0.809061367
        if (JELEM	.eq.	2017	)		coeff = 	0.985735769
        if (JELEM	.eq.	2018	)		coeff = 	1.023324806
        if (JELEM	.eq.	2019	)		coeff = 	1.117337047
        if (JELEM	.eq.	2020	)		coeff = 	1.042617615
        if (JELEM	.eq.	2021	)		coeff = 	1.036455051
        if (JELEM	.eq.	2022	)		coeff = 	1.134480339
        if (JELEM	.eq.	2023	)		coeff = 	0.901206016
        if (JELEM	.eq.	2024	)		coeff = 	1.171647377
        if (JELEM	.eq.	2025	)		coeff = 	0.931678397
        if (JELEM	.eq.	2026	)		coeff = 	0.989108244
        if (JELEM	.eq.	2027	)		coeff = 	1.033516981
        if (JELEM	.eq.	2028	)		coeff = 	0.918424339
        if (JELEM	.eq.	2029	)		coeff = 	1.013061769
        if (JELEM	.eq.	2030	)		coeff = 	1.085265492
        if (JELEM	.eq.	2031	)		coeff = 	0.954000494
        if (JELEM	.eq.	2032	)		coeff = 	1.060693367
        if (JELEM	.eq.	2033	)		coeff = 	0.811556622
        if (JELEM	.eq.	2034	)		coeff = 	0.913753441
        if (JELEM	.eq.	2035	)		coeff = 	0.835209285
        if (JELEM	.eq.	2036	)		coeff = 	0.877581709
        if (JELEM	.eq.	2037	)		coeff = 	0.857285392
        if (JELEM	.eq.	2038	)		coeff = 	0.924753145
        if (JELEM	.eq.	2039	)		coeff = 	0.898231886
        if (JELEM	.eq.	2040	)		coeff = 	0.90608965
        if (JELEM	.eq.	2041	)		coeff = 	1.012875461
        if (JELEM	.eq.	2042	)		coeff = 	0.955115928
        if (JELEM	.eq.	2043	)		coeff = 	0.818889292
        if (JELEM	.eq.	2044	)		coeff = 	1.112683858
        if (JELEM	.eq.	2045	)		coeff = 	0.963900744
        if (JELEM	.eq.	2046	)		coeff = 	0.982430482
        if (JELEM	.eq.	2047	)		coeff = 	0.794814196
        if (JELEM	.eq.	2048	)		coeff = 	0.925635858
        if (JELEM	.eq.	2049	)		coeff = 	0.665768514
        if (JELEM	.eq.	2050	)		coeff = 	0.939507723
        if (JELEM	.eq.	2051	)		coeff = 	0.83937098
        if (JELEM	.eq.	2052	)		coeff = 	0.972615329
        if (JELEM	.eq.	2053	)		coeff = 	0.95123428
        if (JELEM	.eq.	2054	)		coeff = 	1.05571694
        if (JELEM	.eq.	2055	)		coeff = 	0.924270462
        if (JELEM	.eq.	2056	)		coeff = 	0.720585587
        if (JELEM	.eq.	2057	)		coeff = 	0.954442701
        if (JELEM	.eq.	2058	)		coeff = 	0.969650726
        if (JELEM	.eq.	2059	)		coeff = 	0.863731139
        if (JELEM	.eq.	2060	)		coeff = 	1.09049969
        if (JELEM	.eq.	2061	)		coeff = 	0.813494816
        if (JELEM	.eq.	2062	)		coeff = 	1.187018159
        if (JELEM	.eq.	2063	)		coeff = 	0.96079385
        if (JELEM	.eq.	2064	)		coeff = 	0.909650808
        if (JELEM	.eq.	2065	)		coeff = 	0.945273632
        if (JELEM	.eq.	2066	)		coeff = 	0.969968606
        if (JELEM	.eq.	2067	)		coeff = 	1.012996278
        if (JELEM	.eq.	2068	)		coeff = 	0.934551269
        if (JELEM	.eq.	2069	)		coeff = 	0.786744117
        if (JELEM	.eq.	2070	)		coeff = 	0.908244014
        if (JELEM	.eq.	2071	)		coeff = 	0.747874338
        if (JELEM	.eq.	2072	)		coeff = 	1.08741129
        if (JELEM	.eq.	2073	)		coeff = 	0.96093397
        if (JELEM	.eq.	2074	)		coeff = 	1.08233747
        if (JELEM	.eq.	2075	)		coeff = 	0.951222905
        if (JELEM	.eq.	2076	)		coeff = 	0.906046483
        if (JELEM	.eq.	2077	)		coeff = 	1.067100471
        if (JELEM	.eq.	2078	)		coeff = 	0.87110374
        if (JELEM	.eq.	2079	)		coeff = 	0.991552558
        if (JELEM	.eq.	2080	)		coeff = 	0.79990836
        if (JELEM	.eq.	2081	)		coeff = 	1.016755039
        if (JELEM	.eq.	2082	)		coeff = 	1.108292487
        if (JELEM	.eq.	2083	)		coeff = 	1.042522003
        if (JELEM	.eq.	2084	)		coeff = 	1.09555994
        if (JELEM	.eq.	2085	)		coeff = 	0.742424599
        if (JELEM	.eq.	2086	)		coeff = 	1.151882298
        if (JELEM	.eq.	2087	)		coeff = 	1.080338965
        if (JELEM	.eq.	2088	)		coeff = 	1.159415498
        if (JELEM	.eq.	2089	)		coeff = 	1.043541134
        if (JELEM	.eq.	2090	)		coeff = 	1.161575806
        if (JELEM	.eq.	2091	)		coeff = 	0.921696908
        if (JELEM	.eq.	2092	)		coeff = 	0.959260626
        if (JELEM	.eq.	2093	)		coeff = 	1.034535928
        if (JELEM	.eq.	2094	)		coeff = 	1.050790553
        if (JELEM	.eq.	2095	)		coeff = 	1.091416952
        if (JELEM	.eq.	2096	)		coeff = 	0.999847183
        if (JELEM	.eq.	2097	)		coeff = 	1.171156184
        if (JELEM	.eq.	2098	)		coeff = 	0.934208187
        if (JELEM	.eq.	2099	)		coeff = 	0.969859668
        if (JELEM	.eq.	2100	)		coeff = 	1.016454014
        if (JELEM	.eq.	2101	)		coeff = 	0.885114447
        if (JELEM	.eq.	2102	)		coeff = 	0.840109019
        if (JELEM	.eq.	2103	)		coeff = 	1.029527465
        if (JELEM	.eq.	2104	)		coeff = 	0.96003412
        if (JELEM	.eq.	2105	)		coeff = 	0.978627579
        if (JELEM	.eq.	2106	)		coeff = 	1.007282059
        if (JELEM	.eq.	2107	)		coeff = 	0.840227535
        if (JELEM	.eq.	2108	)		coeff = 	0.66789068
        if (JELEM	.eq.	2109	)		coeff = 	0.926815561
        if (JELEM	.eq.	2110	)		coeff = 	1.054974963
        if (JELEM	.eq.	2111	)		coeff = 	1.076780964
        if (JELEM	.eq.	2112	)		coeff = 	0.941070599
        if (JELEM	.eq.	2113	)		coeff = 	1.01077336
        if (JELEM	.eq.	2114	)		coeff = 	1.027845886
        if (JELEM	.eq.	2115	)		coeff = 	0.949788033
        if (JELEM	.eq.	2116	)		coeff = 	1.055233629
        if (JELEM	.eq.	2117	)		coeff = 	0.909641049
        if (JELEM	.eq.	2118	)		coeff = 	1.111883909
        if (JELEM	.eq.	2119	)		coeff = 	1.129504586
        if (JELEM	.eq.	2120	)		coeff = 	1.022561283
        if (JELEM	.eq.	2121	)		coeff = 	1.098513206
        if (JELEM	.eq.	2122	)		coeff = 	0.79800382
        if (JELEM	.eq.	2123	)		coeff = 	0.835885172
        if (JELEM	.eq.	2124	)		coeff = 	0.99245949
        if (JELEM	.eq.	2125	)		coeff = 	1.05915267
        if (JELEM	.eq.	2126	)		coeff = 	0.983056707
        if (JELEM	.eq.	2127	)		coeff = 	0.986909798
        if (JELEM	.eq.	2128	)		coeff = 	0.89130117
        if (JELEM	.eq.	2129	)		coeff = 	0.989545214
        if (JELEM	.eq.	2130	)		coeff = 	0.740229033
        if (JELEM	.eq.	2131	)		coeff = 	0.787845239
        if (JELEM	.eq.	2132	)		coeff = 	0.740845861
        if (JELEM	.eq.	2133	)		coeff = 	1.005968425
        if (JELEM	.eq.	2134	)		coeff = 	1.021499398
        if (JELEM	.eq.	2135	)		coeff = 	0.809110211
        if (JELEM	.eq.	2136	)		coeff = 	1.045505708
        if (JELEM	.eq.	2137	)		coeff = 	1.073555335
        if (JELEM	.eq.	2138	)		coeff = 	0.958268001
        if (JELEM	.eq.	2139	)		coeff = 	0.793777752
        if (JELEM	.eq.	2140	)		coeff = 	0.990611639
        if (JELEM	.eq.	2141	)		coeff = 	1.043688079
        if (JELEM	.eq.	2142	)		coeff = 	1.097768843
        if (JELEM	.eq.	2143	)		coeff = 	0.765710317
        if (JELEM	.eq.	2144	)		coeff = 	0.934156065
        if (JELEM	.eq.	2145	)		coeff = 	0.997390105
        if (JELEM	.eq.	2146	)		coeff = 	0.914268189
        if (JELEM	.eq.	2147	)		coeff = 	0.864424283
        if (JELEM	.eq.	2148	)		coeff = 	1.009405376
        if (JELEM	.eq.	2149	)		coeff = 	0.904674514
        if (JELEM	.eq.	2150	)		coeff = 	1.047500982
        if (JELEM	.eq.	2151	)		coeff = 	0.728606506
        if (JELEM	.eq.	2152	)		coeff = 	0.8977114
        if (JELEM	.eq.	2153	)		coeff = 	1.059968818
        if (JELEM	.eq.	2154	)		coeff = 	0.979713665
        if (JELEM	.eq.	2155	)		coeff = 	0.924740007
        if (JELEM	.eq.	2156	)		coeff = 	0.769302683
        if (JELEM	.eq.	2157	)		coeff = 	0.955774337
        if (JELEM	.eq.	2158	)		coeff = 	0.926812422
        if (JELEM	.eq.	2159	)		coeff = 	0.908834893
        if (JELEM	.eq.	2160	)		coeff = 	0.776714184
        if (JELEM	.eq.	2161	)		coeff = 	1.065071721
        if (JELEM	.eq.	2162	)		coeff = 	0.989748813
        if (JELEM	.eq.	2163	)		coeff = 	1.01523474
        if (JELEM	.eq.	2164	)		coeff = 	0.904232873
        if (JELEM	.eq.	2165	)		coeff = 	0.806039474
        if (JELEM	.eq.	2166	)		coeff = 	0.966600179
        if (JELEM	.eq.	2167	)		coeff = 	0.857873711
        if (JELEM	.eq.	2168	)		coeff = 	1.011353293
        if (JELEM	.eq.	2169	)		coeff = 	0.949892023
        if (JELEM	.eq.	2170	)		coeff = 	0.994330574
        if (JELEM	.eq.	2171	)		coeff = 	0.801122307
        if (JELEM	.eq.	2172	)		coeff = 	0.910458954
        if (JELEM	.eq.	2173	)		coeff = 	0.846206699
        if (JELEM	.eq.	2174	)		coeff = 	1.082355066
        if (JELEM	.eq.	2175	)		coeff = 	1.024650705
        if (JELEM	.eq.	2176	)		coeff = 	0.875452559
        if (JELEM	.eq.	2177	)		coeff = 	1.043594947
        if (JELEM	.eq.	2178	)		coeff = 	1.129499735
        if (JELEM	.eq.	2179	)		coeff = 	0.981403246
        if (JELEM	.eq.	2180	)		coeff = 	0.761090653
        if (JELEM	.eq.	2181	)		coeff = 	1.02962806
        if (JELEM	.eq.	2182	)		coeff = 	0.944098506
        if (JELEM	.eq.	2183	)		coeff = 	1.002266107
        if (JELEM	.eq.	2184	)		coeff = 	1.137223466
        if (JELEM	.eq.	2185	)		coeff = 	0.96389562
        if (JELEM	.eq.	2186	)		coeff = 	0.846954231
        if (JELEM	.eq.	2187	)		coeff = 	1.030541243
        if (JELEM	.eq.	2188	)		coeff = 	1.119123111
        if (JELEM	.eq.	2189	)		coeff = 	1.034257251
        if (JELEM	.eq.	2190	)		coeff = 	0.915677793
        if (JELEM	.eq.	2191	)		coeff = 	1.010531882
        if (JELEM	.eq.	2192	)		coeff = 	0.916090649
        if (JELEM	.eq.	2193	)		coeff = 	1.158174157
        if (JELEM	.eq.	2194	)		coeff = 	0.895381011
        if (JELEM	.eq.	2195	)		coeff = 	0.993698679
        if (JELEM	.eq.	2196	)		coeff = 	1.130049439
        if (JELEM	.eq.	2197	)		coeff = 	0.989669143
        if (JELEM	.eq.	2198	)		coeff = 	0.896042343
        if (JELEM	.eq.	2199	)		coeff = 	0.778653636
        if (JELEM	.eq.	2200	)		coeff = 	0.661787672
        if (JELEM	.eq.	2201	)		coeff = 	0.66423031
        if (JELEM	.eq.	2202	)		coeff = 	0.801550376
        if (JELEM	.eq.	2203	)		coeff = 	0.823949291
        if (JELEM	.eq.	2204	)		coeff = 	0.86024896
        if (JELEM	.eq.	2205	)		coeff = 	0.948415502
        if (JELEM	.eq.	2206	)		coeff = 	0.986191797
        if (JELEM	.eq.	2207	)		coeff = 	1.07508199
        if (JELEM	.eq.	2208	)		coeff = 	0.917710641
        if (JELEM	.eq.	2209	)		coeff = 	0.825090535
        if (JELEM	.eq.	2210	)		coeff = 	1.025985501
        if (JELEM	.eq.	2211	)		coeff = 	0.839642481
        if (JELEM	.eq.	2212	)		coeff = 	1.102301674
        if (JELEM	.eq.	2213	)		coeff = 	0.997037603
        if (JELEM	.eq.	2214	)		coeff = 	1.027850819
        if (JELEM	.eq.	2215	)		coeff = 	1.065039436
        if (JELEM	.eq.	2216	)		coeff = 	0.925376051
        if (JELEM	.eq.	2217	)		coeff = 	1.014147276
        if (JELEM	.eq.	2218	)		coeff = 	0.727885598
        if (JELEM	.eq.	2219	)		coeff = 	0.964380927
        if (JELEM	.eq.	2220	)		coeff = 	0.887447019
        if (JELEM	.eq.	2221	)		coeff = 	1.15868763
        if (JELEM	.eq.	2222	)		coeff = 	0.933394974
        if (JELEM	.eq.	2223	)		coeff = 	0.942116195
        if (JELEM	.eq.	2224	)		coeff = 	0.857104403
        if (JELEM	.eq.	2225	)		coeff = 	0.917595689
        if (JELEM	.eq.	2226	)		coeff = 	0.815349338
        if (JELEM	.eq.	2227	)		coeff = 	0.796468548
        if (JELEM	.eq.	2228	)		coeff = 	1.065295025
        if (JELEM	.eq.	2229	)		coeff = 	1.051176811
        if (JELEM	.eq.	2230	)		coeff = 	0.864991979
        if (JELEM	.eq.	2231	)		coeff = 	1.108517591
        if (JELEM	.eq.	2232	)		coeff = 	0.994048273
        if (JELEM	.eq.	2233	)		coeff = 	1.018745664
        if (JELEM	.eq.	2234	)		coeff = 	0.889192955
        if (JELEM	.eq.	2235	)		coeff = 	1.085017798
        if (JELEM	.eq.	2236	)		coeff = 	0.864243714
        if (JELEM	.eq.	2237	)		coeff = 	0.868782534
        if (JELEM	.eq.	2238	)		coeff = 	0.954902821
        if (JELEM	.eq.	2239	)		coeff = 	1.032205656
        if (JELEM	.eq.	2240	)		coeff = 	1.102184528
        if (JELEM	.eq.	2241	)		coeff = 	0.927024618
        if (JELEM	.eq.	2242	)		coeff = 	1.139830179
        if (JELEM	.eq.	2243	)		coeff = 	1.107652357
        if (JELEM	.eq.	2244	)		coeff = 	1.074066822
        if (JELEM	.eq.	2245	)		coeff = 	0.977583087
        if (JELEM	.eq.	2246	)		coeff = 	0.911756262
        if (JELEM	.eq.	2247	)		coeff = 	0.830104264
        if (JELEM	.eq.	2248	)		coeff = 	0.964444449
        if (JELEM	.eq.	2249	)		coeff = 	1.116877701
        if (JELEM	.eq.	2250	)		coeff = 	1.014856081
        if (JELEM	.eq.	2251	)		coeff = 	0.921961995
        if (JELEM	.eq.	2252	)		coeff = 	0.867114776
        if (JELEM	.eq.	2253	)		coeff = 	1.021810728
        if (JELEM	.eq.	2254	)		coeff = 	0.964604885
        if (JELEM	.eq.	2255	)		coeff = 	0.851493791
        if (JELEM	.eq.	2256	)		coeff = 	0.936515389
        if (JELEM	.eq.	2257	)		coeff = 	0.953755159
        if (JELEM	.eq.	2258	)		coeff = 	1.010127762
        if (JELEM	.eq.	2259	)		coeff = 	0.988134489
        if (JELEM	.eq.	2260	)		coeff = 	0.86357085
        if (JELEM	.eq.	2261	)		coeff = 	1.006734026
        if (JELEM	.eq.	2262	)		coeff = 	0.974313148
        if (JELEM	.eq.	2263	)		coeff = 	1.000015161
        if (JELEM	.eq.	2264	)		coeff = 	0.909267797
        if (JELEM	.eq.	2265	)		coeff = 	0.944677138
        if (JELEM	.eq.	2266	)		coeff = 	0.918649321
        if (JELEM	.eq.	2267	)		coeff = 	0.966477653
        if (JELEM	.eq.	2268	)		coeff = 	0.991712952
        if (JELEM	.eq.	2269	)		coeff = 	0.970236519
        if (JELEM	.eq.	2270	)		coeff = 	1.104802613
        if (JELEM	.eq.	2271	)		coeff = 	0.988313453
        if (JELEM	.eq.	2272	)		coeff = 	0.707435933
        if (JELEM	.eq.	2273	)		coeff = 	0.869689207
        if (JELEM	.eq.	2274	)		coeff = 	0.8912028
        if (JELEM	.eq.	2275	)		coeff = 	0.876306224
        if (JELEM	.eq.	2276	)		coeff = 	0.880140497
        if (JELEM	.eq.	2277	)		coeff = 	0.837846024
        if (JELEM	.eq.	2278	)		coeff = 	0.874361137
        if (JELEM	.eq.	2279	)		coeff = 	0.681373277
        if (JELEM	.eq.	2280	)		coeff = 	1.08179083
        if (JELEM	.eq.	2281	)		coeff = 	0.992357103
        if (JELEM	.eq.	2282	)		coeff = 	0.966215528
        if (JELEM	.eq.	2283	)		coeff = 	1.030802884
        if (JELEM	.eq.	2284	)		coeff = 	1.126734993
        if (JELEM	.eq.	2285	)		coeff = 	0.694093542
        if (JELEM	.eq.	2286	)		coeff = 	0.892202049
        if (JELEM	.eq.	2287	)		coeff = 	1.066893536
        if (JELEM	.eq.	2288	)		coeff = 	1.066920057
        if (JELEM	.eq.	2289	)		coeff = 	0.900290881
        if (JELEM	.eq.	2290	)		coeff = 	0.996440823
        if (JELEM	.eq.	2291	)		coeff = 	1.099058291
        if (JELEM	.eq.	2292	)		coeff = 	0.988364335
        if (JELEM	.eq.	2293	)		coeff = 	1.068787684
        if (JELEM	.eq.	2294	)		coeff = 	0.861232746
        if (JELEM	.eq.	2295	)		coeff = 	0.768981982
        if (JELEM	.eq.	2296	)		coeff = 	1.182727023
        if (JELEM	.eq.	2297	)		coeff = 	0.919220186
        if (JELEM	.eq.	2298	)		coeff = 	0.909627018
        if (JELEM	.eq.	2299	)		coeff = 	1.032125402
        if (JELEM	.eq.	2300	)		coeff = 	0.837898121
        
        
         EA = EA * coeff
         
        USRVAR(jelem,5,1)   = coeff
        USRVAR(jelem,5,2)   = coeff
        USRVAR(jelem,5,3)   = coeff
        USRVAR(jelem,5,4)   = coeff
         
            ! elastic stiffness matrix
            G0      = EA / (2.d0*(1.d0 + nu))
            K11     = EA / (1.D0 - nu * nu)
            K12     = nu * K11
            De(:,1) = (/ K11,  K12, 0.D0/)
            De(:,2) = (/ K12,  K11, 0.D0/)
            De(:,3) = (/0.D0, 0.D0,   G0/)
            
            ! softening parameters
            a1   =  4.d0/(c0*lb)*EA*Gf/(ft*ft)
            if      (istype == 1) then  ! linear softening
              p  =  2.d0
              a2 = -0.5d0
              a3 =  0.0d0
            else if (istype == 2) then  ! exponential softening
              p  =  2.5d0
              a2 =  2.0d0**(5.d0/3.d0) - 3.0d0
              a3 =  0.0d0
            else if (istype == 3) then  ! blinear softening
              p  =  2.0d0
              a2 =  0.03687d0
              a3 =  20.8343d0
            else if (istype == 4) then  ! concrete softening
              p  =  2.0d0
              a2 =  1.3868d0
              a3 =  0.6567d0
            else if (istype == 5) then  ! hyperbolic softening
              p  =  4.0d0 
              a2 =  2.0d0**(7.d0/3.d0) - 4.5d0
              a3 =  0.0d0
            else
              write (*,*) '**error: Softening law No. ', istype, 
     &                    'does not exist!'
            end if
            
            ! integration points
            gp = (/ -1.d0, 1.d0 /) / dsqrt(3.d0)
            gw = (/  1.d0, 1.d0 /)
            
            ! dof interchange
            indexq = (/ 1,2,9, 3,4,10, 5,6,11, 7,8,12 /)
            ! interchange the locations of dofs
            QQ = 0.d0
            do i = 1, 12
              QQ(indexq(i),i) = 1.d0
            end do             
            
            bInitialized = .true.
            
            return
          end subroutine Initialize
      !========================================================================= 
      end module ModelParam

!**********************************************************************************************************
!
      module FEM
!
!**********************************************************************************************************
        use NumKind
        implicit none

        contains      
          !==================shape function and its derivative with xi and eta======================    
          subroutine shapefuc(n, dn_xieta, xi, eta)
          
            implicit none      
            real(rkind) :: n(4), dn_xieta(2, 4), xi, eta

            n(1) = 0.25d0*(1.d0 - xi)*(1.d0 - eta)
            n(2) = 0.25d0*(1.d0 + xi)*(1.d0 - eta)
            n(3) = 0.25d0*(1.d0 + xi)*(1.d0 + eta)
            n(4) = 0.25d0*(1.d0 - xi)*(1.d0 + eta)
            
            dn_xieta(1, 1) = -0.25d0*(1.d0 - eta)
            dn_xieta(1, 2) =  0.25d0*(1.d0 - eta)
            dn_xieta(1, 3) =  0.25d0*(1.d0 + eta)
            dn_xieta(1, 4) = -0.25d0*(1.d0 + eta)
            
            dn_xieta(2, 1) = -0.25d0*(1.d0 - xi)
            dn_xieta(2, 2) = -0.25d0*(1.d0 + xi)
            dn_xieta(2, 3) =  0.25d0*(1.d0 + xi)
            dn_xieta(2, 4) =  0.25d0*(1.d0 - xi)
            
            return 
          end subroutine shapefuc

          !===============traditional b matrix==============================================      
          subroutine b_matrix(nd,bd,b,det_jacb, coords,xi,eta)
          
            implicit none
            real(rkind) :: nd(4), bd(2,4), b(3,8)
            real(rkind) :: jacb(2,2), inv_jacb(2,2), coords(2, 4)
            real(rkind) :: det_jacb, xi, eta
            
            !local varibles
            real(rkind) :: n(4), dn_xieta(2,4), dn_x(4), dn_y(4)
            integer(ikind) :: i, j
             
            ! shape functions 
            call shapefuc(n,dn_xieta,xi,eta)
            nd = n
            
            ! jacob matrix
            jacb = matmul(dn_xieta, transpose(coords))            
            det_jacb = jacb(1,1)*jacb(2,2) - jacb(1,2)*jacb(2,1)
            inv_jacb(1, 1) = jacb(2, 2)
            inv_jacb(1, 2) =-jacb(1, 2)
            inv_jacb(2, 1) =-jacb(2, 1)
            inv_jacb(2, 2) = jacb(1, 1)
            inv_jacb = 1.d0/det_jacb*inv_jacb            
            
            !initialize varibles
            do i = 1,4
              dn_x(i) = inv_jacb(1,1)*dn_xieta(1,i)
     &                + inv_jacb(1,2)*dn_xieta(2,i)
              dn_y(i) = inv_jacb(2,1)*dn_xieta(1,i)
     &                + inv_jacb(2,2)*dn_xieta(2,i)
            end do
            
            ! B matrix for displacement
            b = 0.d0
            do j = 1, 4
              b(1, 2*(j-1) + 1) = dn_x(j)
              b(2, 2*(j-1) + 2) = dn_y(j)
              b(3, 2*(j-1) + 1) = dn_y(j)
              b(3, 2*(j-1) + 2) = dn_x(j)
            end do
            
            ! B matrix for damage
            do j = 1,4
              bd(1,j) = dn_x(j)
              bd(2,j) = dn_y(j)
            end do
          
            return
          end subroutine b_matrix
      
        !********************************************************************
        ! define the dyadic function
          function dyadic(vector1,vector2, vlen)
        !********************************************************************
            integer (ikind) :: vlen, i, j
            real    (rkind) :: vector1(vlen),vector2(vlen)
            real    (rkind) :: dyadic(vlen,vlen)
          
            do i = 1, vlen
              do j = 1, vlen
                dyadic(i,j) = vector1(i) * vector2(j)
              end do
            end do

            return
          end function dyadic

      end module FEM
      

!**********************************************************************************************************
!
      subroutine pfczm(rhs,amatrix,coords,u,svars,jelem)
!
!**********************************************************************************************************
        use NumKind
        use ModelParam
        use FEM
        use ComVarible
        implicit none

        real(rkind):: rhs(12), amatrix(12,12), coords(2,4)
        real(rkind):: svars(4), u(12)
       
        ! local varibles
        real(rkind):: b(3,8), nd(4), bd(2,4)
        real(rkind):: uu(8), dd(4), rd(4), ru(8), kdd(4,4), kuu(8,8)
        real(rkind):: rr(12), kk(12,12)
        real(rkind):: strain(3), stressEff(3)
        real(rkind):: det_jacb, energy_crk, phi, omega, domega, ddomega
        real(rkind):: dalpha, ddalpha, phi_source, dphi_source, dvol
        integer(ikind):: i, j, k, jelem

        real(rkind):: savg, sdif, sdev, smax, smin

        
        
        integer (ikind), parameter :: N_ELEM=100000
        integer (ikind), parameter :: NSTV=18
        
        
        
        
        

                

        
        
        
        ! extrat nodal displacement and damage dofs
        do i = 1, 4
          uu(2*i - 1) = u(3*i - 2)
          uu(2*i    ) = u(3*i - 1)
          dd(i)       = u(3*i)
        end do
        
        ! initialize varibles
        rd  = 0.d0
        kdd = 0.d0
        kuu = 0.d0
        do i = 1, ngp
          do j = 1, ngp      
            call b_matrix(nd,bd,b,det_jacb, coords,gp(i),gp(j))
              
            strain = matmul(b, uu)  ! strain field
            stressEff = matmul(De, strain) ! effective stress

c           max/min pricipal stress
            savg = 0.5*(stressEff(1) + stressEff(2))
            sdif = 0.5*(stressEff(1) - stressEff(2))
            sdev = sqrt(sdif*sdif + stressEff(3)*stressEff(3))
            smax = savg + sdev
            smin = savg - sdev

            ! crack driving force
            k = (i - 1) * 2 + j
            energy_crk = 0.5d0*max(smax, ft)**2/EA
            energy_crk = max(energy_crk, svars(k))
            svars(k) = energy_crk
            
            phi  = dot_product(nd,dd) ! crack phase-field            
            call geometricFunc(dalpha,ddalpha,phi) ! geometric function
            call energeticFunc(omega,domega,ddomega,phi) ! energetic function
     
            phi_source  = domega *energy_crk + Gf/(c0*lb)*dalpha
            dphi_source = ddomega*energy_crk + Gf/(c0*lb)*ddalpha

            ! residual for damage
            dvol=  gw(i)*gw(j)*det_jacb*thk
            rd  =  rd  - dvol*(phi_source*nd + 2.d0*lb*Gf/c0
     &          *  matmul(transpose(bd), matmul(bd, dd)))

            ! element matrices
            kdd =  kdd + dvol*((dphi_source)*dyadic(nd, nd, 4)
     &          +  2.d0*lb*Gf/c0*matmul(transpose(bd),bd))
            
            kuu =  kuu + dvol*matmul(matmul(transpose(b), omega*De), b)

            
            
            
            USRVAR(jelem,1,i)   = omega * stressEff(1)
            USRVAR(jelem,1,i+2) = USRVAR(jelem,1,i)
            
            USRVAR(jelem,2,i)   = omega * stressEff(2)
            USRVAR(jelem,2,i+2) = USRVAR(jelem,2,i)
            
            USRVAR(jelem,3,i)   = omega * stressEff(3)
            USRVAR(jelem,3,i+2) = USRVAR(jelem,3,i)
            
            USRVAR(jelem,4,i)   = phi
            USRVAR(jelem,4,i+2) = USRVAR(jelem,4,i)
            
          end do
        end do        
       
        
        ru = -matmul(kuu,uu) ! applies to hybrid formulation
        
        rr(1:8 ) = ru
        rr(9:12) = rd
        
        kk = 0.d0
        kk(1:8 , 1:8 ) = kuu
        kk(9:12, 9:12) = kdd
        
        rhs     = matmul(transpose(QQ),rr)
        amatrix = matmul(matmul(transpose(QQ),kk),QQ)
  
        
        
        
                
        
    
        
        
        
        return 
      end subroutine pfczm
      
!**********************************************************************************************************
!
      subroutine energeticFunc(omega,domega,ddomega,phi)
!
!**********************************************************************************************************
        use NumKind
        use ModelParam
        implicit none
      
        real(rkind) :: omega, domega, ddomega, phi
        real(rkind) :: fac1, dfac1, ddfac1, fac2, dfac2, ddfac2
      
        fac1    =  (1.d0 - phi)**p
        dfac1   = -p*(1.d0 - phi)**(p - 1.d0); 
        ddfac1  =  p*(p - 1.d0)*(1.d0 - phi)**(p - 2.d0)
        
        fac2   =  fac1   + a1*phi + a1*a2*phi**2.d0 + a1*a2*a3*phi**3.d0
        dfac2  =  dfac1  + a1 + 2.d0*a1*a2*phi + 3.d0*a1*a2*a3*phi**2.d0
        ddfac2  =  ddfac1 + 2.d0*a1*a2 + 6.d0*a1*a2*a3*phi
        
        omega   =  fac1/fac2        
        domega  =  (dfac1*fac2  - fac1*dfac2)/(fac2**2.d0)
        ddomega = ((ddfac1*fac2 - fac1*ddfac2)*fac2 - 2.d0*
     &             (dfac1*fac2 - fac1*dfac2)*dfac2)/(fac2**3.d0)
     
        return
      end subroutine energeticFunc
      
!**********************************************************************************************************
!
      subroutine geometricFunc(dalpha,ddalpha,phi)
!
!**********************************************************************************************************
        use NumKind
        implicit none
        
        real(rkind) :: dalpha, phi, ddalpha
        
        dalpha  = 2.d0 - 2.d0*phi
        ddalpha =-2.d0
        
        return 
      end subroutine geometricFunc  
     
!**********************************************************************************************************
      subroutine UEL(rhs, amatrx, svars, energy, ndofel, nrhs, nsvars,
     &               props, nprops, coords, mcrd, nnode, 
     &               u, du, v, a, jtype, time, dtime, kstep, 
     &               kinc, jelem, params, ndload, jdltyp, adlmag,
     &               predef, npredf, lflags, mlvarx, ddlmag, mdload,
     &               pnewdt, jprops,njprop,period)
!**********************************************************************************************************

        use NumKind
        use ModelParam
        use ComVarible
        INCLUDE 'ABA_PARAM.INC' 

!**********************************************************************************************************
      ! interface of uel, DO NOT change !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! variables passed in
        integer (ikind), intent (in    ) :: ndofel, mlvarx, nrhs, 
     &    nsvars, nprops, mcrd, nnode, jtype, kstep, kinc, jelem, 
     &    ndload,npredf, mdload, njprop
     
        integer (ikind), intent (in    ) :: jdltyp(mdload,*), 
     &    lflags(*), jprops(njprop)
     
        real    (rkind), intent (in    ) :: props(nprops), 
     &    coords(mcrd,nnode), u(ndofel), du(mlvarx,*), v(ndofel), 
     &    a(ndofel), time(2), params(3), adlmag(mdload,*),
     &    ddlmag(mdload,*), predef(2,npredf,nnode), dtime, period
     
        ! variables can be updated
        real    (rkind), intent (in out) :: pnewdt
  
        ! variables to be updated (the update of energy(8) is optional)
        real    (rkind), intent (in out) :: rhs(mlvarx,nrhs), 
     &    amatrx(ndofel,ndofel), svars(nsvars), energy(8)
      ! interface of uel, DO NOT change !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        integer (ikind), parameter :: N_ELEM=100000
        integer (ikind), parameter :: NSTV=18
        
        

         
        
         
       
       
      !********************************************************************************************************
      !                   
      ! user coding to define rhs, amatrx, svars, energy and pnewdt (optional for the last two)
      !
        
!        write (*, *) 'TIME = ', time(2)

        ! initialize parameters, etc.                     
C        if (.not. bInitialized) then
         call Initialize(props, nprops, jtype, jelem)
C        end if
        
        
      IF (JTYPE .EQ. 3) THEN
          call pfczm(rhs(:,1),amatrx,coords,u,svars,jelem)
      ENDIF
        
      IF (JTYPE .EQ. 6) THEN
         CALL UEL6(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1 PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2 KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,NPREDF,
     3 LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPRO,PERIOD)
      ENDIF
        
      return
      end subroutine uel




! USER ELEMENT SUBROUTINE FOR QUADRATIC PLANE STRAIN COHESIVE ELEMENTS
! CODED BY TENG TONG (ASSISTANT PROF.SOUTHEAST UNIVERSITY) SIQI YUAN (PHD CANDIDATR.SOUTHEAST UNIVERSITY)   
C   
C--------------------------------------------------------------------------------
C              
      SUBROUTINE UEL6(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1 PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2 KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,NPREDF,
     3 LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPRO,PERIOD)
C     
C--------------------------------------------------------------------------------
C
      use NumKind
      use ComVarible
      INCLUDE 'ABA_PARAM.INC' !IMPLICIT REAL(A-H O-Z)TM0
C      
C--------------------------------------------------------------------------------
C
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),SVARS(*),
     1 ENERGY(*),COORDS(MCRD,NNODE),U(NDOFEL),DU(MLVARX,*),V(NDOFEL),
     2 A(NDOFEL),TIME(2),PARAMS(*),JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),
     3 DDLMAG(MDLOAD,*),PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)  
C        
C--------------------------------------------------------------------------------
C
       PARAMETER (ZERO = 0.D0, HALF=0.5D0, ONE= 1.0D0, TWO=2.0D0, 
     1     THREE= 3.0D0, TOL=-1E-5, NINTP = 3) 
       PARAMETER(N_ELEM=10000,NSDV=120,NGP=9)
C       
C--------------------------------------------------------------------------------
C
      DIMENSION SC(NDOFEL,NDOFEL),FC(NDOFEL,NRHS),T(MCRD,NRHS),
     1 T_D(MCRD,MCRD),R(MCRD,MCRD),BC(MCRD,NDOFEL),BCT(NDOFEL,MCRD),
     2 DEL(MCRD),TMP(NDOFEL,MCRD),RB(NDOFEL,NDOFEL),DELD(2)
C
C GAUSS INTEGRATION VARIABLES (3 INTEG POINT)
      DIMENSION GAUSS3(3), WEIGHT3(3)
C
C	ARRAYS FOR QUADRATIC LINE ELEMENT
      DIMENSION DNDXI(2), DU_CONT(MCRD), DU_LOC(MCRD)
      DIMENSION DV_CONT(MCRD),DV_LOC(MCRD)
      DIMENSION H(MCRD,4), C_COOR(MCRD,NNODE), PSI(4,NDOFEL) 
      DIMENSION B(MCRD, NDOFEL), BT(NDOFEL, MCRD)
      DIMENSION A1(NDOFEL, MCRD), A2(NDOFEL, NDOFEL) 
      DIMENSION AV_COOR(MCRD, 2)
C
C      
C--------------------------------------------------------------------------------
C
C	INITILIZATION 
C     IMPORTANT!! FORTRAN DOES NOT PUT ZEROS IN THERE AUTOMATICALLY 
C
      CALL KASET2(AMATRX, NDOFEL, NDOFEL)
      CALL KASET1(RHS, MLVARX) 
C
      CALL KASET2(PSI, 4, NDOFEL)
	CALL KASET2(H, MCRD, 4)
      CALL KASET2(AV_COOR, MCRD, 2) 
      CALL KASET2(R, MCRD, MCRD)
      CALL KASET2(RB, NDOFEL, NDOFEL)
      CALL KASET2(T, MCRD, NRHS)
      CALL KASET2(T_D, MCRD, MCRD)
C     PAPARMETER INPUT	 
      WIDTH = 250.0d0	! WIDTH OF ELEMENTS (SAME AS SOLID SECTION WIDTH FOR SOLID ELEMENTS)
C
C--------------------------------------------------------------------------------
C
C	RELATION MATRIX
      PSI(1, 1) = -ONE
      PSI(2, 2) = -ONE
      PSI(3, 3) = -ONE
      PSI(4, 4) = -ONE
      PSI(1, 7) = ONE
      PSI(2, 8) = ONE
      PSI(3, 5) = ONE
      PSI(4, 6) = ONE
C	COMPUTE NODAL COORDINATES IN DEFORMED STATE 
      DO  I=1,MCRD
          DO  J=1, NNODE
              NN=I+(J-1)*MCRD 
              C_COOR(I,J) = COORDS(I,J)  + U(NN)
	        END DO
      END DO
C	REFERENCE COORDINATE SYSTEM (MIDPOINT AVERAGES) 
      AV_COOR(1,1)=ONE/TWO*(C_COOR(1,1)+C_COOR(1,4))
      AV_COOR(2,1)=ONE/TWO*(C_COOR(2,1)+C_COOR(2,4))
      AV_COOR(1,2)=ONE/TWO*(C_COOR(1,2)+C_COOR(1,3))
      AV_COOR(2,2)=ONE/TWO*(C_COOR(2,2)+C_COOR(2,3))
C
C--------------------------------------------------------------------------------
C
C	GAUSSIAN POINT
      GAUSS3(1) = -SQRT(0.6)
      GAUSS3(2) = ZERO 
      GAUSS3(3) = SQRT(0.6)
C
      WEIGHT3(1) = 0.55555555555555
      WEIGHT3(2) = 0.88888888888888
      WEIGHT3(3) = 0.55555555555555
C
C--------------------------------------------------------------------------------
C
C     TRANSFORMATION MATRIX
      CALL KCOORDTRANS (R,COORDS,U,NDOFEL,NNODE,MCRD)
C      
      DO I=1,MCRD
          DO J=1,MCRD
              DO K=0,NNODE-1
                 RB(I+2*K,J+2*K) = R(I,J)
              ENDDO
          ENDDO
      ENDDO
C
C--------------------------------------------------------------------------------
C
      DO IINTP=1,NINTP    !BEGIN LOOP
C
C--------------------------------------------------------------------------------
C
      POINT = GAUSS3(IINTP) 
      WEIGHT = WEIGHT3(IINTP)
C      
C	SHAPE FUNCTION VALUE
      H1 = ONE/TWO*(-POINT + ONE) 
      H2 = ONE/TWO*( POINT + ONE) 
C      
C	DERIVATIVE OF SHAPE FUNCTION VALUE (3X1 MATRIX) 
      DNDXI(1) = -ONE/TWO
      DNDXI(2) = ONE/TWO
C
C	H MATRIX
      H(1,1) = H1
      H(2,2) = H1
      H(1,3) = H2
      H(2,4) = H2
C      
C	B MATRIX: FROM NODAL DISPLACEMENT (U) TO DISPLCAMENE JUMP (DEALT_U)
      CALL KASET2(B, MCRD, NDOFEL) 
      DO L=1, MCRD
         DO J=1, NDOFEL
            DO K=1, NDOFEL/2
               B(L,J) = B(L,J) + H(L,K)*PSI(K,J)
      	    END DO
         END DO
      END DO
C      
C	TRANSPOSED B MATRIX 
      DO L=1, MCRD
        DO J=1, NDOFEL 
           BT(J,L) = B(L,J)
        END DO
      END DO
C
C 	CALCULATE GLOBAL DISPLACEMENT AT INTEGRATION POINT 
C   FROM CONTINUOUS DISPLACEMENT
     	CALL KASET1(DU_CONT, MCRD)
      CALL KASET1(DV_CONT, MCRD)
	    DO L=1, MCRD
         DO J=1, NDOFEL
            DU_CONT(L) = DU_CONT(L) + B(L,J)*U(J)
            DV_CONT(L) = DV_CONT(L) + B(L,J)*V(J)
         END DO
      END DO
C
C     LOCAL COORDINATE SYSTEM		
C     (USE AVERAGE OF DEFORMED X-POSITIONS OF TOP AND BOTTOM)
	    X_XI = ZERO 
      Y_XI = ZERO
	    DO L=1,2
         X_XI = X_XI + DNDXI(L)*AV_COOR(1,L)
         Y_XI = Y_XI + DNDXI(L)*AV_COOR(2,L)
      END DO
C
C	JACOBIAN (VECTOR LENGTH IN XI-DIRECTION) 
      DETJ = SQRT(X_XI**TWO + Y_XI**TWO)
C
C	RELATIVE DISPLACEMENT IN LOCAL COORDINATE SYSTEM
      CALL KASET1(DU_LOC, MCRD) 
      CALL KASET1(DV_LOC, MCRD)
      DO I=1, MCRD
         DO J=1, MCRD
            DU_LOC(I) = DU_LOC(I) + R(I,J)*DU_CONT(J)
            DV_LOC(I) = DV_LOC(I) + R(I,J)*DV_CONT(J)
    	   END DO
      END DO
C      
      DEL  = DU_LOC
      DELD = DV_LOC  
C
C-----------------------------------------------------------------------------------------C
C
C     DETERMINE LOAD CASE
      CALL KSEPLAW(PROPS,DEL,DELD,T,T_D,DTIME,damage,coeff,JELEM)
C     MODE II FRACTURE
C
C-----------------------------------------------------------------------------------------C
C
C     ASSEMBLE AMATRX AND RHS
      BC  = MATMUL(B,RB)
      BCT = TRANSPOSE(BC)
      TMP = MATMUL(BCT,T_D)
      SC  = MATMUL(TMP,BC)
      FC  = MATMUL(BCT,T)
C       
      THICK = WIDTH*WEIGHT*DETJ
C
      CALL KMATRIXSCALAR(AMATRX,SC,THICK,NDOFEL,NDOFEL)
      CALL KMATRIXSCALAR(RHS,-FC,THICK,NDOFEL,NRHS)
      
      IF (IINTP .LE. 2.1D0) THEN
!          
      USRVAR(JELEM,6,IINTP)   = T(1,1)
      USRVAR(JELEM,6,IINTP+2) = USRVAR(JELEM,6,IINTP)
        
      USRVAR(JELEM,7,IINTP)   = T(2,1)
      USRVAR(JELEM,7,IINTP+2) = USRVAR(JELEM,7,IINTP)
      
      USRVAR(JELEM,8,IINTP)   = damage
      USRVAR(JELEM,8,IINTP+2) = USRVAR(JELEM,7,IINTP)
!      
      ENDIF


          
C
C
C--------------------------------------------------------------------------------
C
      ENDDO     !END LOOP
C
      RETURN
      END   
C
C--------------------------------------------------------------------------------
C
C     SUBROUTINE KCOORDSTRANS TO OBTRAIN R, THE TWO DIMENSIONAL
C     TRANSFORMATION MATRIX AND EL_LENGTH, THE ELEMENT SIZE.
      SUBROUTINE KCOORDTRANS(R,COORDS,U,NDOFEL,NNODE,MCRD)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION R(MCRD,MCRD),COORDS(MCRD,NNODE),U(NDOFEL)
      DIMENSION CO_DE(MCRD,NNODE),CO_DE_M(2,2)
C
      DO I=1,MCRD
          DO J=1,NNODE
              CO_DE(I,J)=COORDS(I,J)+U(2*(J-1)+I)
          END DO
      END DO
C
! CALCULATE OF THE DIRECTIONAL COSINE & THE TRANSFORMATION MATRIX
      D_X=CO_DE(1,2)-CO_DE(1,1)
      D_Y=CO_DE(2,2)-CO_DE(2,1)
      EL_LENGTH=(D_X**2+D_Y**2)**0.5D0
      COS_A=D_X/EL_LENGTH
      SIN_A=D_Y/EL_LENGTH
      R(1,1)=COS_A
      R(1,2)=SIN_A
      R(2,1)=-SIN_A
      R(2,2)=COS_A   
      RETURN
      END        
C
C--------------------------------------------------------------------------------
C 
C     SUBROUTINE KMATRIX_PLUSSCALAR TO MULTIPLY A MATRIX
C     WITH A SCALAR NUMBER
      SUBROUTINE KMATRIXSCALAR(A,B,C,N,M)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION A(N,M),B(N,M)
      DO I=1,N
         DO J=1,M
            A(I,J)=A(I,J)+C*B(I,J)
         ENDDO
      ENDDO
C      
      RETURN
      END
C
C--------------------------------------------------------------------------------
C 
      SUBROUTINE KSEPLAW(PROPS,DEL,DELD,T,T_D,DTIME,damage,coeff,JELEM)
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION PROPS(9),DEL(2),DELD(2),T(2,1),T_D(2,2)
C-----------------------------------------------------------------------------------------C      
      coeff = 1.0
      SIGMA_MAX = 2.0d0 
      TAU_MAX   = 2.0d0 * coeff
      DN = 0.05d0
      DT = 0.05d0    ! twice of DN
C      
      Q = 1.0 d0
      R = 0.0
      PHI_N = DEXP(1.D0)*SIGMA_MAX*DN
      
      DELT = ABS(DEL(1))
      DELN = DEL(2)
      
      if (del(1) .GE. 0) then
          sign_dt = 1.0d0
      else
          sign_dt = -1.0d0
      end if
      
      if ( del(2) .Ge. 0) then
          sign_dn = 1.0d0
      ELSE
          sign_dn = -1.0d0
      end if
C-----------------------------------------------------------------------------------------C      

C-----------------------------------------------------------------------------------------C      
      if (DELN .gt. 0.0d0) then
C-----------------------------------------------------------------------------------------C        
          T22_TERM1 = SIGMA_MAX * DEXP(1.0d0)
          T22_TERM2 = DELN / DN
          T22_TERM3 = DEXP( -1.0D0 * DELN / DN )
          T22_TERM4 = DEXP( -1.0D0 * DELT * DELT / DT / DT )
	        T(2,1) = T22_TERM1 * T22_TERM2 * T22_TERM3 * T22_TERM4

          TD22_TERM1 = SIGMA_MAX * DEXP(1.0d0)
          TD22_TERM2 = 1.0D0/DN - DELN/DN/DN
          TD22_TERM3 = DEXP( -1.0D0 * DELN/ DN)
	        T_D(2,2) = TD22_TERM1 * TD22_TERM2 * TD22_TERM3 
          
          TD21_TERM1 = SIGMA_MAX * DEXP(1.0d0)
          TD21_TERM2 = DELN / DN
          TD21_TERM3 = DEXP( -1.0D0 * DELN / DN )
          TD21_TERM4 = DEXP( -1.0D0 * DELT * DELT / DT / DT )
          TD21_TERM5 = ( -2.0D0 * DELT * sign_dt / DT / DT )
	        T_D(2,1) = TD21_TERM1 * TD21_TERM2 * TD21_TERM3 * TD21_TERM4 * 
     &        TD21_TERM5

          T11_TERM1 = TAU_MAX * SQRT(2.0D0 * DEXP(1.0d0) )
          T11_TERM2 = DELT / DT
          T11_TERM3 = DEXP( -1.0D0 * DELT*DELT / DT / DT)
          T11_TERM4 = DEXP( -1.0D0 * DELN / DN)
	    T(1,1) = T11_TERM1 * T11_TERM2 * T11_TERM3 * sign_dt * T11_TERM4
      
          TD11_TERM1 = TAU_MAX * SQRT(2.0D0 * DEXP(1.0d0) )
          TD11_TERM2 = 1.0D0/DT - 2.0D0*DELT*DELT/DT/DT/DT
          TD11_TERM3 = DEXP( -1.0D0 * DELT*DELT / DT / DT)
          TD11_TERM4 = DEXP( -1.0D0 * DELN / DN)
	        T_D(1,1) = TD11_TERM1 * TD11_TERM2 * TD11_TERM3 * TD11_TERM4
          
          TD12_TERM1 = TAU_MAX * SQRT(2.0D0 * DEXP(1.0d0) )
          TD12_TERM2 = DELT / DT
          TD12_TERM3 = DEXP( -1.0D0 * DELT*DELT / DT / DT)
          TD12_TERM4 = DEXP( -1.0D0 * DELN / DN)
          TD12_TERM5 =( -1.0D0  / DN)
          T_D(1,2) = TD12_TERM1 * TD12_TERM2 * TD12_TERM3 * TD12_TERM4 *
     &    TD12_TERM5
      else
C-----------------------------------------------------------------------------------------C      
	        T(2,1) = 1.0D4 * DELN
	        T_D(2,2) = 1.0D4
C
          T11_TERM1 = TAU_MAX * SQRT(2.0D0 * DEXP(1.0d0) )
          T11_TERM2 = DELT / DT
          T11_TERM3 = DEXP( -1.0D0 * DELT*DELT / DT / DT)
	        T(1,1) = T11_TERM1 * T11_TERM2 * T11_TERM3 * sign_dt
      
          TD11_TERM1 = TAU_MAX * SQRT(2.0D0 * DEXP(1.0d0) )
          TD11_TERM2 = 1.0D0/DT - 2.0D0*DELT*DELT/DT/DT/DT
          TD11_TERM3 = DEXP( -1.0D0 * DELT*DELT / DT / DT)
	        T_D(1,1) = TD11_TERM1 * TD11_TERM2 * TD11_TERM3
C
          if (T_D(1,1) .gt. 0.0d0) then    
              TAD0_TERM1 = TAU_MAX * SQRT(2.0D0 * DEXP(1.0d0) )
              TAD0_TERM2 = 1.0D0/DT
              TAD0_TERM3 = DEXP( -1.0D0 * DELN / DN)
              TAD0 =  TAD0_TERM1 * TAD0_TERM2 * TAD0_TERM3
              damage = 1.0d0 - abs(T_D(1,1) / TAD0)
          else
              damage = 1.0d0
          ENDIF   
C 
          FRIC_COEFF = 0.00D0
          T(1,1) = T(1,1) - FRIC_COEFF * T(2,1) * damage * sign_dt   
          T_D(1,2) = T_D(1,2) - FRIC_COEFF * T_D(2,2) * damage * sign_dt
          T_D(1,1) = T_D(1,1) - FRIC_COEFF * T_D(2,1) * damage * sign_dt
      end if 
C
      ZETA=0.01D0
c
    	T(1,1)   = T(1,1)   + ZETA*TAU_MAX*DELD(1)/DT
    	T_D(1,1) = T_D(1,1) + ZETA*TAU_MAX/DT/DTIME     
c      
    	T(2,1)   = T(2,1)   + ZETA*SIGMA_MAX*DELD(2)/DN
    	T_D(2,2) = T_D(2,2) + ZETA*SIGMA_MAX/DN/DTIME   
      
C     write(6,*) 'T(1,1)',T(1,1)     time(2)

C-----------------------------------------------------------------------------------------C
      RETURN
      END     
C
C--------------------------------------------------------------------------------
C  
      SUBROUTINE KASET2(DMATRIX,IDIMX,IDIMY)
C
      INCLUDE 'ABA_PARAM.INC'
      PARAMETER (ZERO = 0.0D0) 
      DIMENSION DMATRIX(IDIMX,IDIMY)
C      
      DO I = 1,IDIMX 
         DO J = 1,IDIMY
            DMATRIX(I,J) = ZERO
         END DO 
      END DO
C
      RETURN 
      END
C
C--------------------------------------------------------------------------------
C
      SUBROUTINE KASET1(DMATRIX, IDIMX) 
      INCLUDE 'ABA_PARAM.INC'
      PARAMETER (ZERO = 0.0D0) 
      DIMENSION DMATRIX(IDIMX)
C
      DO I=1, IDIMX
         DMATRIX(I) = ZERO 
      END DO
C      
      RETURN 
      END
C
C--------------------------------------------------------------------------------
C                        



C
C ==============================================================
C !!! NOTE: N_ELEM has to be changed according to the UEL !!!!!
C ==============================================================
C
       SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,STRAN,DSTRAN,
     2 TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,MATERL,NDI,NSHR,NTENS,
     3 NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,CELENT,
     4 DFGRD0,DFGRD1,NOEL,NPT,KSLAY,KSPT,KSTEP,KINC)
C
      use ComVarible
      INCLUDE 'ABA_PARAM.INC'
C
       CHARACTER*80 CMNAME
       DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),
     2 DDSDDT(NTENS),DRPLDE(NTENS),
     3 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     4 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
C 
       PARAMETER (ONE=1.0,TWO=2.0,THREE=3.0,SIX=6.0, HALF =0.5,
     1 N_ELEM=100000,NSDV=120)
       DATA NEWTON,TOLER/40,1.D-2/ 
C       
       
C 
C ----------------------------------------------------------- 
C          Material properties
C ----------------------------------------------------------- 
C          PROPS(1) - Young's modulus 
C          PROPS(2) - Poisson ratio 
C ----------------------------------------------------------- 
C
C	Elastic properties
C
       EMOD=PROPS(1)
       ENU=PROPS(2)
       EG=EMOD/(TWO*(ONE+ENU))
       EG2=EG*TWO
       ELAM=EG2*ENU/(ONE-TWO*ENU)
C
C	Stiffness tensor
C
       DO K1=1, NTENS
        DO K2=1, NTENS
         DDSDDE(K2, K1)=0.0
        END DO
       END DO
C
       DO K1=1, NDI
        DO K2=1, NDI
         DDSDDE(K2, K1)=ELAM
        END DO
        DDSDDE(K1, K1)=EG2+ELAM
       END DO 
C
       DO K1=NDI+1, NTENS
        DDSDDE(K1, K1)=EG
       END DO
C
C	Calculate Stresses
C
       DO K1=1, NTENS
        DO K2=1, NTENS
         STRESS(K2)=STRESS(K2)+DDSDDE(K2, K1)*DSTRAN(K1)
        END DO
       END DO 
C
       NELEMAN=NOEL-N_ELEM

       DO I=1,NSTATV
        STATEV(I)=USRVAR(NELEMAN,I,NPT)
       END DO
C       
       RETURN
       END      